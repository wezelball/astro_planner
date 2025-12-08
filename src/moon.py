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


def moon_phase_info(t_snapshot, eph):
    """
    Returns:
        illuminated_pct  (float 0–100)
        phase_label      ("Waxing" or "Waning")
        moon_age_days    (float)
    """
    ts = t_snapshot.ts

    # Skyfield moon phase illumination (0–1)
    illuminated = almanac.fraction_illuminated(eph, 'moon', t_snapshot)
    illuminated_pct = float(illuminated * 100.0)

    # Determine waxing vs waning:
    # 0=new → 0.5=full → 1=new
    phase = almanac.moon_phase(eph, t_snapshot).radians
    # Phase angle grows from 0 → 2π each cycle
    # Waxing if 0 < phase < π
    waxing = (0 < phase < 3.14159265)
    phase_label = "Waxing" if waxing else "Waning"

    # Moon age in days
    cycle_len_days = 29.530588
    moon_age_days = (phase / (2 * 3.14159265)) * cycle_len_days

    return illuminated_pct, phase_label, moon_age_days



