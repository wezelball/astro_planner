import numpy as np
from astropy.time import Time
from math import sin, cos, asin, atan2, radians, degrees, sqrt, tan

# -----------------------------------------------------------------------------
# UTILITY FUNCTIONS
# -----------------------------------------------------------------------------

def wrap_angle_deg(angle):
    """Wrap angle to the range [0, 360)."""
    return angle % 360.0

def wrap_angle_pm180(angle):
    """Wrap angle to [-180, +180)."""
    a = angle % 360.0
    return a - 360.0 if a >= 180.0 else a

# -----------------------------------------------------------------------------
# ALT/AZ CORE
# -----------------------------------------------------------------------------

def compute_local_sidereal_time(jd, longitude_deg):
    """
    Compute Local Sidereal Time (LST) using IAU 2006 expressions.

    jd : float (Julian Date)
    longitude_deg : observer longitude, +East of Greenwich
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
    Convert RA/Dec to Alt/Az for a given observer location & time.

    time_utc : string or astropy Time
    lon_deg  : geographic longitude (+E)
    lat_deg  : geographic latitude (+N)
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
    Kasten & Young 1989 airmass model.
    Input altitude in degrees.
    """
    if alt_deg <= 0:
        return np.inf

    z = 90.0 - alt_deg
    return 1.0 / (cos(radians(z)) + 0.50572 * (6.07995 + z)**-1.6364)

def parallactic_angle(ha_deg, dec_deg, lat_deg):
    """
    Parallactic angle (degrees).
    Positive = North to East.
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
    Full Alt/Az solution including:
      - Alt
      - Az
      - Hour angle
      - Airmass
      - Parallactic angle
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
