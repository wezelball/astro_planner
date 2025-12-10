# test_moon_info.py
from skyfield.api import Loader
from datetime import datetime
from zoneinfo import ZoneInfo
from moon import moon_phase_info

LOCAL_TZ = ZoneInfo("America/New_York")

load = Loader('./skyfield_data')
ts = load.timescale()
eph = load('de421.bsp')
#eph = load('de430.bsp')   # instead of de421.bsp
#eph = load('de440s.bsp')   # instead of de421.bsp

# Your test: 2025-12-09 22:00 local
local_dt = datetime(2025, 12, 9, 22, 0, tzinfo=LOCAL_TZ)
utc_dt = local_dt.astimezone(ZoneInfo("UTC"))
t = ts.utc(utc_dt.year, utc_dt.month, utc_dt.day, utc_dt.hour, utc_dt.minute)

illum_pct, waxing, phase_name, age_days = moon_phase_info(eph, ts, t)
print("illum_pct:", illum_pct)
print("waxing:", waxing)
print("phase_name:", phase_name)
print("age_days:", age_days)
print("utc snapshot:", utc_dt.isoformat())
