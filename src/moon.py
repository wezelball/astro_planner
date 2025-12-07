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

def moon_separation_deg(eph, ts, observer_geocentric, t, star):
    """
    Angular separation (deg) between Moon and star at time t as seen from observer_geocentric.
    star should be a Skyfield Star (or similar) object.
    """
    moon_astrom = observer_geocentric.at(t).observe(eph['moon']).apparent()
    star_astrom = observer_geocentric.at(t).observe(star).apparent()
    # Use Skyfield's separation_from
    sep = star_astrom.separation_from(moon_astrom)
    return float(sep.degrees)

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
