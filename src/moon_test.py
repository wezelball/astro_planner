# src/moon_test.py
from skyfield.api import Loader, Topos, Star, wgs84
from datetime import datetime
from moon import moon_position, moon_separation_deg

# --- Configure this section for your site / test ---
LAT = 37.8      # example: your site latitude
LON = -77.9     # example: your site longitude
ELEV_M = 116    # your site elevation in meters

# Snapshot to test (UTC). Adjust if you want local -> convert to UTC yourself.
# Example: 2025-12-07 22:00 local (if local is EST you may need to shift to UTC)
SNAPSHOT = datetime(2025, 12, 7, 22, 0)  # year, month, day, hour, minute (UTC assumed here)

# Test objects (use RA/Dec in degrees). Replace with exact values you used for comparison.
test_objects = [
    {"name": "NGC0598", "ra_deg": 23.4625, "dec_deg": 30.6603},
    {"name": "NGC2403", "ra_deg": 114.1242, "dec_deg": 65.6945},
    {"name": "NGC7331", "ra_deg": 339.0146, "dec_deg": 34.4208},
    {"name": "NGC1023", "ra_deg": 40.1125, "dec_deg": 39.6181},
]

# --- Load Skyfield data ---
load = Loader('./skyfield_data')
ts = load.timescale()
eph = load('de421.bsp')

# Build geocentric observer
earth = eph['earth']
observer = earth + wgs84.latlon(latitude_degrees=LAT,
                                longitude_degrees=LON,
                                elevation_m=ELEV_M)

# Build Skyfield time for snapshot (UTC)
t_snapshot = ts.utc(SNAPSHOT.year, SNAPSHOT.month, SNAPSHOT.day,
                    SNAPSHOT.hour, SNAPSHOT.minute)

# Moon position
moon_alt, moon_az, moon_apparent = moon_position(eph, observer, t_snapshot)
print(f"Snapshot (UTC): {SNAPSHOT.isoformat()}  Observer: lat={LAT}, lon={LON}, elev={ELEV_M}m")
print(f"Moon altitude: {moon_alt:.6f} deg, azimuth: {moon_az:.6f} deg\n")

# Compute separations
for obj in test_objects:
    star = Star(ra_hours=obj["ra_deg"] / 15.0, dec_degrees=obj["dec_deg"])
    sep = moon_separation_deg(eph, observer, t_snapshot, star)
    print(f"{obj['name']:8s}  RA={obj['ra_deg']:8.4f}  Dec={obj['dec_deg']:8.4f}  Sep={sep:.6f} deg  (rounded -> {int(sep)}°)")
