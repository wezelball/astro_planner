import streamlit as st
from datetime import time
from src.config_loader import load_config_example
#from src.catalog import load_catalog_sample
from src.catalog import load_openngc_catalog

from src.horizon import parse_horizon_file
from src.planner import Planner
from src.optics import Optics
import pandas as pd
from math import isnan

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
max_mag = st.sidebar.number_input("Max magnitude", value=12.0, step=0.5)
min_alt = st.sidebar.number_input(
    "Min altitude (deg)",
    value=float(config['minima']['min_altitude_deg']),
    step=1.0
)
fov_min = st.sidebar.slider(
    "FOV fill fraction (min, max)",
    0.0, 1.0,
    (config['minima']['min_fov_fill'], config['minima']['max_fov_fill'])
)

# Sidebar date & moon
st.sidebar.header("Date & Moon")
date = st.sidebar.date_input("Date to plan (local)")

# Select snapshot time (UTC)
snapshot_time = st.sidebar.time_input(
    "Snapshot Time (UTC)",
    value=time(hour=4, minute=0),  # default: 04:00 UTC ≈ local midnight
    key="snapshot_time_picker"
)

# Extract hour and minute
hour_utc = snapshot_time.hour
minute_utc = snapshot_time.minute

st.sidebar.checkbox("Show Clear Sky Clock (if within 48h)", value=False)

# Define selectable types (matches Planner.filter logic)
#selectable_types = ["G", "RfN", "HII", "OCl", "PN", "Neb", "Cl+N", "SNR", "EmN", "GCl"]
selectable_types = ["G", "OCl", "PN", "Nebula", "SNR", "GCl"]
selected_type = st.sidebar.selectbox("Select object type", selectable_types)

# Sidebar sorting
st.sidebar.header("Sorting")
sort_by_visibility = st.sidebar.checkbox("Sort by visibility duration (hours)", value=False)

# Load catalogs and horizon
#catalog = load_catalog_sample("data/messier_sample.csv")
catalog_df = load_openngc_catalog("data/NGC.csv")

horizon = parse_horizon_file("horizon/horizon_sample.txt")

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
    date=date,
    optics=optics,          # Optics object we built earlier
    max_magnitude=max_mag,
    min_altitude=min_alt,
    fov_fill_range=fov_min,
    object_list=object_list,
    selected_type=selected_type,
    hour_utc = hour_utc
)

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
    display_columns = ['name', 'type', 'mag', 'size_deg', 'alt_deg', 'az_deg', 'visible_hours']
    st.table(results[display_columns])
else:
    st.write("No candidates found with current filters.")

st.info(
    "This is an MVP skeleton. The src modules contain the structure and placeholder functions.\n"
    "Next steps: implement precise ephemeris calculations using Skyfield, visibility sampling, scoring, and Streamlit UI refinements."
)
