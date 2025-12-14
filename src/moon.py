# src/moon.py
"""
moon.py - helper functions for Moon altitude and Moon-object separation.

Functions:
- moon_alt_az(eph, ts, observer_geocentric, t): returns (alt_deg, az_deg) for moon at time t
- moon_alt_az_vector(eph, ts, observer_geocentric, times): same but for array of times (returns arrays)
- moon_separation_deg(eph, ts, observer_geocentric, t, star): angular separation (deg) between Moon and star at time t
- moon_separation_deg_vector(eph, ts, observer_geocentric, times, star): vectorized form
Notes:
- observer_geocentric must be created as: observer = eph['earth'] + wgs84.latlon(lat, lon, elevation_m=...)
- star should be a skyfield.api.Star object (or any Skyfield position-like object)
"""

from skyfield.api import Star
import numpy as np
from skyfield.api import Star, wgs84
#from skyfield.api import Star
from skyfield.api import Loader
from skyfield import almanac
from datetime import timedelta

def moon_position_topocentric(eph, observer, t_snapshot):
    """Return (moon_alt_deg, moon_az_deg, moon_apparent)"""
    moon = eph['moon']
    astrometric_moon = observer.at(t_snapshot).observe(moon)
    moon_apparent = astrometric_moon.apparent()
    alt, az, _ = moon_apparent.altaz()
    return float(alt.degrees), float(az.degrees), moon_apparent

def moon_separation_deg_from_apparent(obj_apparent, moon_apparent):
    """Return angular separation (degrees) between an object apparent and moon_apparent.
       obj_apparent should be a Skyfield apparent position (ICRF/ICRF-like) or a Skyfield
       star 'apparent()' result; moon_apparent is the same type for the Moon."""
    # This assumes both are topocentric/apparent at the same time and observer.
    return float(obj_apparent.separation_from(moon_apparent).degrees)

def moon_position(eph, observer, t_snapshot):
    """
    Compute topocentric apparent Moon position at a given snapshot.

    Parameters
    ----------
    eph : Skyfield ephemeris (Loader(...)('de421.bsp') or similar)
    observer : geocentric observer (eph['earth'] + wgs84.latlon(...))
    t_snapshot : Skyfield Time object

    Returns
    -------
    (moon_alt_deg, moon_az_deg, moon_apparent)
    - moon_alt_deg, moon_az_deg are floats in degrees
    - moon_apparent is a Skyfield Apparent object for separation calculations
    """
    moon = eph['moon']
    astrometric_moon = observer.at(t_snapshot).observe(moon)
    moon_apparent = astrometric_moon.apparent()
    alt, az, _ = moon_apparent.altaz()
    return float(alt.degrees), float(az.degrees), moon_apparent

def moon_alt_az(eph, ts, observer_geocentric, t):
    """
    Return (alt_deg, az_deg) for the Moon at Skyfield time t as seen from observer_geocentric.
    t may be a single Skyfield Time (not array).
    """
    moon = eph['moon']
    astrometric = observer_geocentric.at(t).observe(moon)
    apparent = astrometric.apparent()
    alt, az, _ = apparent.altaz()
    return float(alt.degrees), float(az.degrees)

def moon_alt_az_vector(eph, ts, observer_geocentric, times):
    """
    times can be a Skyfield Time array (vectorized).
    Returns (alt_deg_array, az_deg_array) as numpy arrays.
    """
    moon = eph['moon']
    astrometric = observer_geocentric.at(times).observe(moon)
    apparent = astrometric.apparent()
    alt, az, _ = apparent.altaz()
    return np.array(alt.degrees), np.array(az.degrees)

def moon_separation_deg(eph, observer, t_snapshot, star):
    """
    Compute angular separation (degrees) between a catalog object and the Moon,
    evaluated topocentrically at the same observer/time.

    Parameters
    ----------
    eph : Skyfield ephemeris
    observer : geocentric observer (eph['earth'] + wgs84.latlon(...))
    t_snapshot : Skyfield Time object
    star : Skyfield Star object (constructed from RA/Dec in J2000)

    Returns
    -------
    sep_deg : float (degrees)
    """
    # Compute topocentric apparent position of the Moon once
    moon = eph['moon']
    astrometric_moon = observer.at(t_snapshot).observe(moon)
    moon_apparent = astrometric_moon.apparent()

    # Compute topocentric apparent position of the star at the same time
    astrometric_star = observer.at(t_snapshot).observe(star)
    star_apparent = astrometric_star.apparent()

    # Separation (degrees)
    sep_deg = float(star_apparent.separation_from(moon_apparent).degrees)
    return sep_deg

def moon_separation_deg_vector(eph, ts, observer_geocentric, times, star):
    """
    Vectorized angular separations for an array of Skyfield times.
    Returns numpy array of separations (degrees).
    """
    moon_astrom = observer_geocentric.at(times).observe(eph['moon']).apparent()
    star_astrom = observer_geocentric.at(times).observe(star).apparent()
    # separation_from supports arrays and returns an Angle array
    sep = star_astrom.separation_from(moon_astrom)
    return np.array(sep.degrees)


def moon_phase_info(eph, ts, t_snapshot):
    """
    Return a tuple (illuminated_pct, waxing_bool, phase_name, moon_age_days).

    - eph: loaded skyfield ephemeris (Loader(...)(...))
    - ts: skyfield timescale (loader.timescale())
    - t_snapshot: skyfield Time object (ts.utc(...))

    Uses:
      - almanac.fraction_illuminated(eph, 'moon', t) for illuminated fraction.
      - almanac.moon_phases(eph) + almanac.find_discrete to find nearest new moon
        and compute moon age (days).
    """
    # Illuminated fraction: 0..1
    try:
        frac = almanac.fraction_illuminated(eph, 'moon', t_snapshot)
        illuminated_pct = float(frac * 100.0)
    except Exception as e:
        illuminated_pct = None

    # Waxing / waning and phase name: use discrete moon phase events
    phase_func = almanac.moon_phases(eph)

    # Search a 40-day window before and after t_snapshot to find events.
    # 30 days is a lunar month ~29.53d; 40 is safe.
    t0 = ts.utc(t_snapshot.utc_datetime() - timedelta(days=40))
    t1 = ts.utc(t_snapshot.utc_datetime() + timedelta(days=40))

    times, events = almanac.find_discrete(t0, t1, phase_func)

    # events is array of integers: 0=new, 1=first quarter, 2=full, 3=last quarter
    # Find index of the last event at or before t_snapshot
    idx = None
    for i, tt in enumerate(times):
        if tt.tt > t_snapshot.tt - 1e-9:  # tt.tt is TT; a small tolerance
            idx = i - 1
            break
    if idx is None:
        idx = len(times) - 1

    if idx < 0:
        # fallback: use first event
        idx = 0

    last_event_time = times[idx]
    last_event_type = int(events[idx])

    # Phase name for the last event and next event
    phase_names = {0: "New Moon", 1: "First Quarter", 2: "Full Moon", 3: "Last Quarter"}
    last_phase_name = phase_names.get(last_event_type, "Unknown")

    # Moon age = days since last new moon (if last_event_type==0, simple).
    # If the last event isn't a New Moon, step backward to find the most recent New Moon in the list.
    # Find last new-moon time:
    last_new_idx = None
    for j in range(idx, -1, -1):
        if int(events[j]) == 0:
            last_new_idx = j
            break
    if last_new_idx is None:
        # find forward until you find a new moon; fallback to computing difference to last_event_time
        # but this is unlikely because our 40-day window should include a new moon.
        moon_age_days = None
    else:
        new_time = times[last_new_idx]
        # compute age in days using TT seconds difference converted to days (safer via .utc_datetime)
        dt = t_snapshot.utc_datetime() - new_time.utc_datetime()
        moon_age_days = dt.total_seconds() / 86400.0

    # waxing/waning: find next major phase after last_new_idx (or compare relative phases)
    # Simpler approach: find the next phase event after t_snapshot (if it exists),
    # and if that next event is Full (2) and last was New (0), it's waxing; if next is New, waning.
    waxing = None
    # find next event index after t_snapshot
    next_idx = None
    for k, tt in enumerate(times):
        if tt.tt > t_snapshot.tt + 1e-9:
            next_idx = k
            break
    if next_idx is not None:
        prev_ev = events[idx] if idx >= 0 else None
        next_ev = events[next_idx]
        # If we're between New (0) and Full (2) it's waxing; between Full and New it's waning.
        # More robust: if next_ev in (1,2) and last_event_type in (0,1) -> waxing
        if int(next_ev) in (1, 2) and int(last_event_type) in (0, 1):
            waxing = True
        else:
            # otherwise consider waning
            waxing = False
    else:
        # fallback: set waxing based on phase angle sign using almanac.moon_phase (radians)
        try:
            phase_rad = float(almanac.moon_phase(eph, t_snapshot).radians)
            # if phase in (0, pi) -> waxing; else waning (same logic as earlier)
            waxing = (0 < phase_rad < 3.141592653589793)
        except Exception:
            waxing = None

    return illuminated_pct, waxing, last_phase_name, moon_age_days