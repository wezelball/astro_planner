"""
altaz.py
Alt/Az coordinate utilities.

This module provides low-level astronomical coordinate transformations,
including conversion from equatorial coordinates (RA/Dec) to horizontal
coordinates (Altitude/Azimuth), along with related quantities such as
airmass and parallactic angle.

All angles are in degrees unless otherwise stated.

Conventions:
- Longitude: positive East of Greenwich
- Latitude: positive North
- Azimuth: 0° = North, increasing eastward (astronomical convention)
- Time: UTC
"""

import numpy as np
from astropy.time import Time
from math import sin, cos, asin, atan2, radians, degrees, sqrt, tan

# -----------------------------------------------------------------------------
# UTILITY FUNCTIONS
# -----------------------------------------------------------------------------

def wrap_angle_deg(angle):
    """
    Wrap an angle to the range [0, 360).

    Args:
        angle (float): Angle in degrees.

    Returns:
        float: Wrapped angle in degrees.
    """
    return angle % 360.0

def wrap_angle_pm180(angle):
    """
    Wrap an angle to the range [-180, +180).

    Args:
        angle (float): Angle in degrees.

    Returns:
        float: Wrapped angle in degrees.
    """
    a = angle % 360.0
    return a - 360.0 if a >= 180.0 else a

# -----------------------------------------------------------------------------
# ALT/AZ CORE
# -----------------------------------------------------------------------------

def compute_local_sidereal_time(jd, longitude_deg):
    """
    Compute Local Sidereal Time (LST) using the IAU 2006 expression.

    Args:
        jd (float): Julian Date (UTC).
        longitude_deg (float): Observer longitude in degrees
            (positive East of Greenwich).

    Returns:
        float: Local Sidereal Time in degrees [0, 360).
    """
    T = (jd - 2451545.0) / 36525.0

    # GMST in seconds
    gmst_sec = (
        67310.54841 +
        (876600.0 * 3600 + 8640184.812866) * T +
        0.093104 * T**2 -
        6.2e-6 * T**3
    )

    gmst_deg = (gmst_sec / 240.0) % 360.0   # convert to degrees

    lst = wrap_angle_deg(gmst_deg + longitude_deg)
    return lst

def radec_to_altaz(ra_deg, dec_deg, time_utc, lat_deg, lon_deg):
    """
    Convert equatorial coordinates (RA/Dec) to horizontal coordinates (Alt/Az).

    Args:
        ra_deg (float): Right Ascension in degrees.
        dec_deg (float): Declination in degrees.
        time_utc (str or astropy.time.Time): Observation time (UTC).
        lat_deg (float): Observer latitude in degrees (+North).
        lon_deg (float): Observer longitude in degrees (+East).

    Returns:
        tuple:
            alt_deg (float): Altitude in degrees.
            az_deg (float): Azimuth in degrees (0° = North, increasing East).
            ha_deg (float): Hour angle in degrees (-180 to +180).
    """

    # Ensure astropy Time
    t = Time(time_utc, scale="utc")

    jd = t.jd
    lst = compute_local_sidereal_time(jd, lon_deg)

    # Hour Angle (degrees)
    ha = wrap_angle_pm180(lst - ra_deg)

    # Convert to radians
    ha_r  = radians(ha)
    dec_r = radians(dec_deg)
    lat_r = radians(lat_deg)

    # Altitude
    alt_r = asin(
        sin(dec_r) * sin(lat_r) +
        cos(dec_r) * cos(lat_r) * cos(ha_r)
    )

    # Azimuth (measured from North, increasing Eastward)
    y = -sin(ha_r) * cos(dec_r)
    x = (sin(dec_r) - sin(alt_r) * sin(lat_r)) / (cos(alt_r) * cos(lat_r))
    az_r = atan2(y, x)

    az_deg = wrap_angle_deg(degrees(az_r))
    alt_deg = degrees(alt_r)

    return alt_deg, az_deg, ha

# -----------------------------------------------------------------------------
# AIR MASS & PARALLACTIC ANGLE
# -----------------------------------------------------------------------------

def airmass_kasten_young(alt_deg):
    """
    Compute airmass using the Kasten & Young (1989) model.

    Valid for altitudes above the horizon.

    Args:
        alt_deg (float): Altitude in degrees.

    Returns:
        float: Airmass value, or infinity if altitude ≤ 0°.
    """
    if alt_deg <= 0:
        return np.inf

    z = 90.0 - alt_deg
    return 1.0 / (cos(radians(z)) + 0.50572 * (6.07995 + z)**-1.6364)

def parallactic_angle(ha_deg, dec_deg, lat_deg):
    """
    Compute the parallactic angle.

    The parallactic angle is the angle between the great circle
    through the object and the zenith, and the hour circle of the object.

    Args:
        ha_deg (float): Hour angle in degrees.
        dec_deg (float): Declination in degrees.
        lat_deg (float): Observer latitude in degrees.

    Returns:
        float: Parallactic angle in degrees.
               Positive values indicate rotation from North toward East.
    """
    ha = radians(ha_deg)
    dec = radians(dec_deg)
    lat = radians(lat_deg)

    y = sin(ha)
    x = tan(lat) * cos(dec) - sin(dec) * cos(ha)
    return degrees(atan2(y, x))

# -----------------------------------------------------------------------------
# MASTER WRAPPER
# -----------------------------------------------------------------------------

def altaz_from_radec(ra_deg, dec_deg, time_utc, lat_deg, lon_deg):
    """
    Compute a full set of horizontal-coordinate quantities from RA/Dec.

    This is a convenience wrapper that returns:
      - Altitude
      - Azimuth
      - Hour angle
      - Airmass
      - Parallactic angle

    Args:
        ra_deg (float): Right Ascension in degrees.
        dec_deg (float): Declination in degrees.
        time_utc (str or astropy.time.Time): Observation time (UTC).
        lat_deg (float): Observer latitude in degrees (+North).
        lon_deg (float): Observer longitude in degrees (+East).

    Returns:
        dict: Dictionary containing:
            - alt_deg
            - az_deg
            - ha_deg
            - airmass
            - parallactic_angle_deg
    """
    alt, az, ha = radec_to_altaz(ra_deg, dec_deg, time_utc, lat_deg, lon_deg)

    am = airmass_kasten_young(alt)
    pa = parallactic_angle(ha, dec_deg, lat_deg)

    return {
        "alt_deg": alt,
        "az_deg": az,
        "ha_deg": ha,
        "airmass": am,
        "parallactic_angle_deg": pa
    }
