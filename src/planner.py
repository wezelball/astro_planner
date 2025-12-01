# src/planner.py
import math
import pandas as pd
from datetime import datetime, timedelta, time
import numpy as np

from .horizon import horizon_altitude_at

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

    def plan(self, date_obj, optics_name, max_magnitude=12.0, min_altitude=20.0, fov_fill_range=(0.15,0.85), sampling_min=5):
        """
        Main planner function.
        - date_obj: a datetime.date object (local date)
        - optics_name: name of optics entry from config
        - sampling_min: sampling resolution in minutes (default 5)
        Returns pandas DataFrame of candidate targets with visibility windows and score.
        """
        self._ensure_skyfield()

        # find optics entry
        optics_list = self.config.get('optics', [])
        optics_entry = None
        for o in optics_list:
            if o.get('name') == optics_name:
                optics_entry = o
                break
        if optics_entry is None:
            raise ValueError(f"Optics named '{optics_name}' not found in config.")

        # compute diag FOV (arcmin)
        diag_fov_arcmin = self._compute_fov_for_optics(optics_entry)

        # Prepare skyfield objects
        ts = self._ts
        eph = self._ephemeris
        site = self._site_topos()
        earth = eph['earth']
        location = earth + site
        sun = eph['sun']
        moon = eph['moon']

        # define time window: find astronomical night around the date
        # compute times for 0:00 through 23:59 local time in UTC-aware TimeScale
        # We'll sample from 00:00 to 23:59 local and then filter by sun altitude < -18 deg
        tz = None
        try:
            import zoneinfo
            # try to get timezone from config if provided
            tzname = self.config.get('location', {}).get('timezone', None)
            if tzname:
                tz = zoneinfo.ZoneInfo(tzname)
        except Exception:
            tz = None

        # produce local datetime for start and end of the date
        local_midnight = datetime(date_obj.year, date_obj.month, date_obj.day, 0, 0, 0)
        # create list of times every sampling_min minutes over 24 hours (local)
        times_local = [local_midnight + timedelta(minutes=i*sampling_min) for i in range(0, int(24*60/sampling_min))]
        # convert to UTC-aware times for skyfield by using ts.utc with naive dt (assumed UTC) is not correct.
        # Skyfield's ts.utc can accept (year, month, day, hour, minute, second) in UTC. To be robust without timezone info,
        # we'll treat provided local times as UTC if no timezone is available. If timezone set in config, convert.
        times_ts = []
        for dt in times_local:
            if tz:
                utc_dt = dt.replace(tzinfo=tz).astimezone(tz=__import__('datetime').timezone.utc)
            else:
                # assume system local is UTC (best-effort); skyfield will interpret as UTC
                utc_dt = dt
            times_ts.append(ts.utc(utc_dt.year, utc_dt.month, utc_dt.day, utc_dt.hour, utc_dt.minute, utc_dt.second))
        t_array = ts.utc([t.utc_datetime().year for t in times_ts],
                         [t.utc_datetime().month for t in times_ts],
                         [t.utc_datetime().day for t in times_ts],
                         [t.utc_datetime().hour for t in times_ts],
                         [t.utc_datetime().minute for t in times_ts],
                         [t.utc_datetime().second for t in times_ts])

        # compute sun altitudes to get astronomical night mask
        astromask = []
        sun_apparent = (location.at(t_array).observe(sun)).apparent()
        alt_sun, az_sun, dist = sun_apparent.altaz()
        alt_sun_deg = alt_sun.degrees
        for a in alt_sun_deg:
            astromask.append(a < -18.0)

        candidates = []
        # Precompute moon positions and illumination across times
        moonpos = location.at(t_array).observe(moon).apparent()
        moon_altaz = moonpos.altaz()
        moon_alt_deg = moon_altaz[0].degrees
        # moon separation and illumination: use skyfield's phase function if available
        try:
            from skyfield import almanac
            eph_for_phase = eph
            moon_illum = [almanac.fraction_of_moon_illuminated(eph_for_phase, t) for t in t_array]
            # fraction is 0..1; convert to percent
            moon_illum_pct = [f*100.0 for f in moon_illum]
        except Exception:
            moon_illum_pct = [0.0 for _ in t_array]

        # process targets
        for _, row in self.catalog.iterrows():
            mag = float(row.get('magnitude', 99.0))
            if mag > max_magnitude:
                continue
            ra_deg, dec_deg = self._parse_target_ra_dec(row)
            # compute apparent alt/az over times
            # Skyfield: create a fixed RA/Dec object by using skyfield's Star object via ra/dec
            from skyfield.api import Star
            star = Star(ra_hours=ra_deg/15.0, dec_degrees=dec_deg)
            astrom = location.at(t_array).observe(star).apparent()
            alt, az, d = astrom.altaz()
            alt_deg = alt.degrees
            az_deg = az.degrees

            # Build mask where object is above horizon (taking into account user's local horizon)
            ok_mask = []
            for i,(a_deg, azd) in enumerate(zip(alt_deg, az_deg)):
                # within astronomical night
                if not astromask[i]:
                    ok_mask.append(False)
                    continue
                # above local horizon plus min altitude
                horizon_alt = horizon_altitude_at(self.horizon, azd)
                min_alt_buffer = float(self.config.get('minima', {}).get('min_altitude_deg', min_altitude))
                if a_deg < (horizon_alt + min_alt_buffer):
                    ok_mask.append(False)
                    continue
                # Moon constraints for this time step
                moon_pct = moon_illum_pct[i]
                moon_alt = moon_alt_deg[i]
                # compute moon separation
                # need moon RA/DEC at this time: approximate by observing moon from earth at t_array[i]
                m = location.at(t_array[i]).observe(moon).apparent()
                m_ra, m_dec, _ = m.radec()
                # convert to degrees
                m_ra_deg = m_ra.hours * 15.0
                m_dec_deg = m_dec.degrees
                sep = angular_separation(ra_deg, dec_deg, m_ra_deg, m_dec_deg)
                # moon rules: pass if (illumination <= max) OR (moon below horizon) OR (sep >= min_sep)
                cfg_moon = self.config.get('moon', {})
                max_illum = float(cfg_moon.get('max_illumination_pct', 30.0))
                min_sep = float(cfg_moon.get('min_angular_separation_deg', 30.0))
                if (moon_pct <= max_illum) or (moon_alt < 0) or (sep >= min_sep):
                    ok_mask.append(True)
                else:
                    ok_mask.append(False)

            ok_mask = np.array(ok_mask, dtype=bool)

            # apply FOV fill filter: use major axis (worst-case)
            major_arcmin = float(row.get('major_arcmin', 0.0))
            fov_fill = major_arcmin / diag_fov_arcmin if diag_fov_arcmin > 0 else 0.0
            if fov_fill < fov_fill_range[0] or fov_fill > fov_fill_range[1]:
                continue

            # Now find contiguous windows in ok_mask where value True
            windows = []
            i = 0
            N = len(ok_mask)
            while i < N:
                if ok_mask[i]:
                    start = i
                    while i < N and ok_mask[i]:
                        i += 1
                    end = i-1
                    # convert indices to times
                    t_start = t_array[start].utc_datetime()
                    t_end = t_array[end].utc_datetime()
                    windows.append( (t_start, t_end, (end-start+1)*sampling_min/60.0) )  # duration in hours
                else:
                    i += 1

            if not windows:
                continue

            # choose the longest window for summary
            best_window = max(windows, key=lambda w: w[2])
            duration_h = best_window[2]
            window_str = f"{best_window[0].isoformat()} to {best_window[1].isoformat()}"

            # compute a basic score (lower is better) using altitude, mag, fov_fill, moon proximity (at mid-window)
            mid_idx = int(( ( (t_array[0].utc_datetime() - t_array[0].utc_datetime()).total_seconds() ) ))
            # better: pick index of middle time between start and end
            mid_time = best_window[0] + timedelta(hours=duration_h/2.0)
            # find nearest index
            diffs = [abs((tt.utc_datetime() - mid_time).total_seconds()) for tt in t_array]
            mid_idx = int(np.argmin(diffs))
            mid_alt = alt_deg[mid_idx]
            # normalize components
            norm_alt = 1.0 - (mid_alt / 90.0)    # lower is better (higher altitude => lower penalty)
            norm_mag = min(max((mag + 2.0) / 18.0, 0.0), 1.0)  # map mag roughly -2..16 to 0..1
            norm_fov = abs(0.5 - min(max(fov_fill,0.0),1.0)) * 2.0  # 0 if ideal 0.5 fill, 1 if extreme
            # moon proximity penalty: find moon separation at mid_idx
            # compute sep at mid
            m = location.at(t_array[mid_idx]).observe(moon).apparent()
            m_ra, m_dec, _ = m.radec()
            sep_mid = angular_separation(ra_deg, dec_deg, m_ra.hours*15.0, m_dec.degrees)
            norm_moon = 1.0 - min(sep_mid / 180.0, 1.0)  # closer => higher penalty
            # weights
            w_alt=0.4; w_mag=0.35; w_fov=0.15; w_moon=0.1
            score = w_alt*norm_alt + w_mag*norm_mag + w_fov*norm_fov + w_moon*norm_moon

            candidates.append({
                'Name': row.get('name'),
                'Type': row.get('type'),
                'Mag': mag,
                'Size_arcmin': major_arcmin,
                'FOV_fill': round(fov_fill,3),
                'Visible_window': window_str,
                'Duration_h': round(duration_h,2),
                'Score': round(score,3),
            })

        if not candidates:
            return pd.DataFrame([])
        df = pd.DataFrame(candidates)
        df = df.sort_values('Score')
        return df
