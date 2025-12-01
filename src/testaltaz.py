from ephemeris import load_ephemeris, get_body_radec
from altaz import altaz_from_radec

# ----------------------------------------------------------------------
# Load ephemeris and timescale
# ----------------------------------------------------------------------
ts, eph = load_ephemeris(download_if_missing=True)

# ----------------------------------------------------------------------
# Observer location
# ----------------------------------------------------------------------
lat = 35.0   # degrees North
lon = -80.0  # degrees East

# ----------------------------------------------------------------------
# UTC time for calculation
# ----------------------------------------------------------------------
time_utc = "2025-01-01T03:00:00"

# ----------------------------------------------------------------------
# Example object: M31 (RA/Dec)
# ----------------------------------------------------------------------
obj = {"ra": 10.6847083, "dec": 41.26875}

# Get RA/Dec from ephemeris (for a dictionary, just returns the RA/Dec)
ra, dec = get_body_radec(eph, obj, ts, ts.utc(2025, 1, 1, 3, 0, 0))

# ----------------------------------------------------------------------
# Compute Alt/Az
# ----------------------------------------------------------------------
altaz = altaz_from_radec(ra, dec, time_utc, lat, lon)

print("RA  :", ra)
print("Dec :", dec)
print("ALT :", altaz["alt_deg"])
print("AZ  :", altaz["az_deg"])
print("Airmass:", altaz["airmass"])
print("Parallactic angle:", altaz["parallactic_angle_deg"])
