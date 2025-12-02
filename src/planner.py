# src/planner.py
import math
#import pandas as pd
#from datetime import datetime, timedelta, time
#import numpy as np
#from .horizon import horizon_altitude_at

import pandas as pd
from datetime import datetime, timedelta
from src.ephemeris import get_default_ephemeris, get_body_radec
from src.altaz import altaz_from_radec

# We try to import skyfield; if not available we will raise an informative error at runtime.
try:
    from skyfield.api import Loader, Topos, load
    from skyfield.almanac import dark_twilight_day
    SKYFIELD_AVAILABLE = True
except Exception as e:
    SKYFIELD_AVAILABLE = False
    _SKYFIELD_IMPORT_ERROR = e

# Utility: compute diagonal FOV in arcminutes given focal length (mm) and sensor size (mm)
def compute_fov_arcmin(focal_length_mm, sensor_width_mm, sensor_height_mm, reducer_ratio=1.0):
    eff_focal = focal_length_mm * reducer_ratio
    # FOV in radians: 2 * arctan(sensor_dimension / (2*focal_length))
    w_rad = 2.0 * math.atan(sensor_width_mm / (2.0 * eff_focal))
    h_rad = 2.0 * math.atan(sensor_height_mm / (2.0 * eff_focal))
    # diagonal
    diag_rad = math.sqrt(w_rad**2 + h_rad**2)
    diag_deg = math.degrees(diag_rad)
    return diag_deg * 60.0  # arcminutes

def angular_separation(ra1_deg, dec1_deg, ra2_deg, dec2_deg):
    # all in degrees
    ra1 = math.radians(ra1_deg)
    ra2 = math.radians(ra2_deg)
    dec1 = math.radians(dec1_deg)
    dec2 = math.radians(dec2_deg)
    # spherical law of cosines
    return math.degrees(math.acos(
        math.sin(dec1)*math.sin(dec2) + math.cos(dec1)*math.cos(dec2)*math.cos(ra1-ra2)
    ))

def sexagesimal_to_degrees(hms_str):
    # Accepts "HH:MM:SS" or "±DD:MM:SS"
    parts = [p.strip() for p in str(hms_str).split(':')]
    if len(parts) == 3:
        a, b, c = parts
        sign = -1 if str(a).strip().startswith('-') else 1
        a = abs(float(a))
        deg = (a + float(b)/60.0 + float(c)/3600.0) * (15.0 if float(a) <= 24 else 1.0)
        # note: for RA hms -> degrees multiply hours by 15
        return sign * deg
    else:
        # try simple float
        return float(hms_str)

class Planner:
    def __init__(self, config, catalog_df, horizon_pts, tzinfo=None):
        self.config = config
        self.catalog = catalog_df.copy()
        self.horizon = horizon_pts
        self.tzinfo = tzinfo  # timezone object (optional)

        # Prepare skyfield loader (lazy)
        self._sf_load = None
        self._ts = None
        self._ephemeris = None

        if SKYFIELD_AVAILABLE:
            # use a temporary local data directory for skyfield
            self._sf_load = Loader('.skyfield_data')
            self._ts = self._sf_load.timescale()
            # eph: IERS or de421
            try:
                self._ephemeris = self._sf_load('de421.bsp')
            except Exception:
                # fallback to online load (may fail if offline)
                self._ephemeris = load('de421.bsp')

    def _ensure_skyfield(self):
        if not SKYFIELD_AVAILABLE:
            raise ImportError(f"Skyfield is required for accurate ephemeris calculations. Import error: {_SKYFIELD_IMPORT_ERROR}")

    def _site_topos(self):
        lat = self.config['location']['latitude']
        lon = self.config['location']['longitude']
        elev = self.config['location'].get('elevation_m', 0)
        return Topos(latitude_degrees=lat, longitude_degrees=lon, elevation_m=elev)

    def _compute_fov_for_optics(self, optics_entry):
        focal = optics_entry.get('focal_length_mm')
        w = optics_entry.get('sensor_width_mm')
        h = optics_entry.get('sensor_height_mm')
        reducer_ratio = float(optics_entry.get('reducer_ratio', 1.0)) if optics_entry.get('has_focal_reducer', False) else 1.0
        return compute_fov_arcmin(focal, w, h, reducer_ratio)

    def _parse_target_ra_dec(self, row):
        # convert RA H:M:S to degrees (RA) and Dec D:M:S to degrees (Dec)
        ra_hms = row.get('ra_hms')
        dec_dms = row.get('dec_dms')
        # For RA, skyfield expects hours or degrees; we'll convert to degrees
        # For dec, use degrees
        def hms_to_deg(hms):
            parts = [p for p in str(hms).split(':')]
            if len(parts) == 3:
                h = float(parts[0])
                m = float(parts[1])
                s = float(parts[2])
                return (abs(h) + m/60.0 + s/3600.0) * (15.0 if abs(h) <= 24 else 1.0) * ( -1 if str(hms).strip().startswith('-') else 1)
            else:
                return float(hms)
        def dms_to_deg(dms):
            parts = [p for p in str(dms).split(':')]
            if len(parts) == 3:
                d = float(parts[0])
                m = float(parts[1])
                s = float(parts[2])
                sign = -1 if str(dms).strip().startswith('-') else 1
                return sign * (abs(d) + m/60.0 + s/3600.0)
            else:
                return float(dms)
        ra_deg = hms_to_deg(ra_hms)
        dec_deg = dms_to_deg(dec_dms)
        return ra_deg, dec_deg

    def plan(
        self,
        date,
        optics,
        max_magnitude=12.0,
        min_altitude=20.0,
        fov_fill_range=(0.0, 1.0),
        object_list=None
    ):
        import pandas as pd
        from skyfield.api import Loader, Star, Topos

        if object_list is None:
            print("No catalog objects provided!")
            return pd.DataFrame()

        # ---------- PRELOAD SKYFIELD (only once) ----------
        load = Loader('./skyfield_data')
        ts = load.timescale()
        eph = load('de421.bsp')

        # Observing time (local midnight shifted to UTC)
        t = ts.utc(date.year, date.month, date.day, 4)

        # Observer position
        lat = self.config["location"]["latitude"]
        lon = self.config["location"]["longitude"]
        elev = self.config["location"]["elevation_m"]

        observer = eph['earth'] + Topos(
            latitude_degrees=lat,
            longitude_degrees=lon,
            elevation_m=elev
        )

        candidates = []

        counts = {
            "total": 0,
            "no_size": 0,
            "filtered_fov": 0,
            "filtered_mag": 0,
            "filtered_alt": 0,
            "passed": 0
        }

        # ---------- PROCESS CATALOG ----------
        for obj in object_list:
            counts["total"] += 1
            name = obj.get("name", "Unknown")

            # --- size -----------------------------------------------------------------
            size_deg = obj.get("size_deg")
            if size_deg is None:
                counts["no_size"] += 1
                continue

            # --- magnitude ------------------------------------------------------------
            mag = obj.get("magnitude")
            if mag is not None and mag > max_magnitude:
                counts["filtered_mag"] += 1
                continue

            # --- FOV fill -------------------------------------------------------------
            try:
                fov_fill = optics.fov_fill_fraction(size_deg)
            except:
                counts["filtered_fov"] += 1
                continue

            if not (fov_fill_range[0] <= fov_fill <= fov_fill_range[1]):
                counts["filtered_fov"] += 1
                continue

            #print(f"Object: {name}: size_deg={size_deg}, fov_fill={fov_fill:.4f}")

            # --- Skyfield alt/az ------------------------------------------------------
            ra_hours = obj["ra_deg"] / 15.0
            dec_deg = obj["dec_deg"]

            star = Star(ra_hours=ra_hours, dec_degrees=dec_deg)

            apparent = observer.at(t).observe(star).apparent()
            alt, az, dist = apparent.altaz()  # CORRECT — no observer passed

            alt_deg = alt.degrees
            az_deg = az.degrees

            if alt_deg < min_altitude:
                counts["filtered_alt"] += 1
                continue

            # --- Passed all filters ---------------------------------------------------
            counts["passed"] += 1

            obj_out = obj.copy()
            obj_out["alt_deg"] = alt_deg
            obj_out["az_deg"] = az_deg
            obj_out["fov_fill"] = fov_fill
            obj_out["visible_hours"] = None
            # rename magnitude to mag for UI consistency
            obj_out["mag"] = obj_out["magnitude"]

            candidates.append(obj_out)

        # ---------- DEBUG SUMMARY ----------
        print("Catalog summary:")
        for k, v in counts.items():
            print(f"  {k}: {v}")

        return pd.DataFrame(candidates)
