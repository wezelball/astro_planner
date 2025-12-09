import streamlit as st
from datetime import datetime, date, time
from zoneinfo import ZoneInfo   # Python 3.9+
from src.config_loader import load_config_example
#from src.catalog import load_catalog_sample
from src.catalog import load_openngc_catalog

from src.horizon import parse_horizon_file
from src.planner import Planner
from src.optics import Optics
from src.moon import moon_phase_info, moon_position_topocentric
import pandas as pd
from math import isnan
from src.plotting import plot_sky_polar
from skyfield.api import Loader, wgs84

LOCAL_TZ = ZoneInfo("America/New_York")

@st.cache_data
def load_opengc_catalog(path="data/NGC.csv"):
    df = pd.read_csv(path, sep=";")

    objects = []
    for _, row in df.iterrows():
        name = row["Name"]

        # RA (convert from HH:MM:SS.s)
        try:
            ra_h, ra_m, ra_s = row["RA"].split(":")
            ra_deg = 15 * (float(ra_h) + float(ra_m)/60 + float(ra_s)/3600)
        except:
            continue

        # Dec (convert from ±DD:MM:SS)
        try:
            sign = -1 if row["Dec"].strip().startswith("-") else 1
            d = row["Dec"].replace("+", "").replace("-", "")
            dec_d, dec_m, dec_s = d.split(":")
            dec_deg = sign * (float(dec_d) + float(dec_m)/60 + float(dec_s)/3600)
        except:
            continue

        # Use V-Mag if present
        vmag = row.get("V-Mag", None)
        if isinstance(vmag, str) and vmag.strip() == "":
            vmag = None
        try:
            vmag = float(vmag)
        except:
            vmag = None

        # Size (convert arcmin → degrees)
        maj = float(row["MajAx"]) if not isnan(row["MajAx"]) else None
        minr = float(row["MinAx"]) if not isnan(row["MinAx"]) else None
        size_deg = (maj / 60.0) if maj else None

        objects.append({
            "name": name,
            "ra_deg": ra_deg,
            "dec_deg": dec_deg,
            "magnitude": vmag,
            "size_deg": size_deg,
            "type": row["Type"],
        })

    return objects

# Load once
object_list = load_opengc_catalog()

# -------------------------------------------------------
# Debug: check how many OpenNGC objects have size data
# -------------------------------------------------------
num_objects = len(object_list)
#num_with_size = sum(1 for o in object_list if o["size_deg"] is not None)
#num_with_mag  = sum(1 for o in object_list if o["magnitude"] is not None)

#st.sidebar.write("Catalog summary:")
st.sidebar.write(f"Total objects loaded: {num_objects}")
#st.sidebar.write(f"Objects with size: {num_with_size}")
#st.sidebar.write(f"Objects with V-magnitude: {num_with_mag}")
# -------------------------------------------------------
# Debug: check how many OpenNGC objects have size data
# -------------------------------------------------------

st.set_page_config(page_title="Astro Planner", layout="wide")
st.title("Astro Planner — Night Target Selector (MVP)")

# Load example config and catalogs
config_path = "config/config.toml"
config = load_config_example(config_path)

st.sidebar.header("Configuration")
st.sidebar.write(f"Location: {config['location']['latitude']}, {config['location']['longitude']} (elev {config['location']['elevation_m']} m)")

# Select optics/camera
optics_names = [o['name'] for o in config['optics']]
optics_choice = st.sidebar.selectbox("Select optics/camera:", optics_names)

# Sidebar filters
st.sidebar.header("Filters")

max_mag = st.sidebar.slider(
    "Max Magnitude",
    min_value=-2.0,
    max_value=20.0,
    value=12.0,
    step=0.1
)

min_alt = st.sidebar.slider(
    "Min Altitude (deg)",
    min_value=0.0,
    max_value=90.0,
    value=20.0,
    step=0.5
)

fov_min = st.sidebar.slider(
    "FOV fill fraction (min, max)",
    0.0, 1.0,
    (config['minima']['min_fov_fill'], config['minima']['max_fov_fill'])
)

# Sidebar date & moon
st.sidebar.header("Date & Moon")
date = st.sidebar.date_input("Date to plan (local)")

# Moon separation slider
moon_sep_min = st.sidebar.slider(
    "Minimum Moon Separation (°)",
    min_value=0,
    max_value=90,
    value=30,
    step=5
)

# Select snapshot time (UTC)
snapshot_time = st.sidebar.time_input(
    "Snapshot Time (local civil time)",
    value=time(hour=4, minute=0),  # default: 04:00 UTC ≈ local midnight
    key="snapshot_time_picker"
)

# Build LOCAL datetime from selected date + time
#local_dt = datetime.combine(date, snapshot_time)
#local_dt = local_dt.replace(tzinfo=LOCAL_TZ)
local_dt = datetime.combine(date, snapshot_time).replace(tzinfo=LOCAL_TZ)

# Convert to UTC
utc_dt = local_dt.astimezone(ZoneInfo("UTC"))

# Pass UTC date/time components to planner
date_utc = utc_dt.date()
hour_utc = utc_dt.hour
minute_utc = utc_dt.minute

# ------------------------------------------------------
# Moon Phase Sidebar Widget
# ------------------------------------------------------
load = Loader('./skyfield_data')
ts = load.timescale()
t_snapshot = ts.utc(
    date.year, date.month, date.day,
    utc_dt.hour, utc_dt.minute
)
eph = load('de421.bsp')

illum_pct, phase_label, age_days = moon_phase_info(t_snapshot, eph)

st.sidebar.subheader("Moon Phase")
st.sidebar.write(f"Illumination: **{illum_pct:.1f}%**")
st.sidebar.write(f"Phase: **{phase_label}**")
st.sidebar.write(f"Moon age: **{age_days:.1f} days**")

# Moon alt/az at snapshot time
load = Loader('./skyfield_data')
ts = load.timescale()
eph = load('de421.bsp')

t_snapshot = ts.utc(utc_dt.year, utc_dt.month, utc_dt.day,
                    utc_dt.hour, utc_dt.minute)

lat = config["location"]["latitude"]
lon = config["location"]["longitude"]
elev = config["location"]["elevation_m"]
observer = eph['earth'] + wgs84.latlon(lat, lon, elevation_m=elev)

moon_alt_deg, moon_az_deg, _ = moon_position_topocentric(eph, observer, t_snapshot)

st.sidebar.write(f"**Moon Alt:** {moon_alt_deg:.1f}°")
st.sidebar.write(f"**Moon Az:** {moon_az_deg:.1f}°")

# Circular progress indicator (illumination fraction)
st.sidebar.progress(illum_pct / 100.0)

# Define selectable types (matches Planner.filter logic)
#selectable_types = ["G", "RfN", "HII", "OCl", "PN", "Neb", "Cl+N", "SNR", "EmN", "GCl"]
selectable_types = ["G", "OCl", "PN", "Nebula", "SNR", "GCl"]
selected_type = st.sidebar.selectbox("Select object type", selectable_types)

# Sidebar sorting
st.sidebar.header("Sorting")
sort_by_visibility = st.sidebar.checkbox("Sort by visibility duration (hours)", value=False)

# Load catalogs and horizon
catalog_df = load_openngc_catalog("data/NGC.csv")

horizon = parse_horizon_file("horizon/GSpring_horizon.txt")

# Optics instance
optics_dict = next(o for o in config['optics'] if o['name'] == optics_choice)
optics = Optics(
    telescope_aperture_mm=optics_dict['telescope_aperture_mm'],
    telescope_focal_length_mm=optics_dict['telescope_focal_length_mm'],
    camera_sensor_width_mm=optics_dict['camera_sensor_width_mm'],
    camera_sensor_height_mm=optics_dict['camera_sensor_height_mm'],
    camera_pixel_width=optics_dict['camera_pixel_width'],
    camera_pixel_height=optics_dict['camera_pixel_height'],
    focal_reducer=optics_dict.get('focal_reducer', 1.0)
)

# Instantiate planner
planner = Planner(config, catalog_df, horizon)

results = planner.plan(
    date=date_utc,
    optics=optics,          # Optics object we built earlier
    max_magnitude=max_mag,
    min_altitude=min_alt,
    fov_fill_range=fov_min,
    object_list=object_list,
    selected_type=selected_type,
    hour_utc = hour_utc,
    minute_utc = minute_utc,
    horizon = horizon,
    moon_sep_min=moon_sep_min
)

# Plot results
st.subheader("Sky Plot")

fig = plot_sky_polar(results, horizon)
st.plotly_chart(fig, width = "stretch")

# Sort results if needed
if not results.empty:
    if sort_by_visibility:
        results = results.sort_values(by='visible_hours', ascending=False)
    else:
        results = results.sort_values(by='mag', ascending=True)

# Display results
st.subheader("Candidate Targets")
if not results.empty:
    st.write(f"Found {len(results)} candidates for {date}")
    display_columns = ['name', 'type', 'mag', 'size_arcmin', 'alt_deg', 'az_deg', 'visible_hours', 'moon_sep_deg']

    # keep numeric rounding, then select only the columns we want to show,
    # and finally convert to strings for left-justified display
    df_left_justified = results.round(2)
    df_left_justified = df_left_justified.loc[:, display_columns].astype(str)

    st.dataframe(df_left_justified, width="stretch")
else:
    st.write("No candidates found with current filters.")


#st.info(
#    "This is an MVP skeleton. The src modules contain the structure and placeholder functions.\n"
#    "Next steps: implement precise ephemeris calculations using Skyfield, visibility sampling, scoring, and Streamlit UI refinements."
#)
