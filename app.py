import streamlit as st
from src.config_loader import load_config_example
from src.catalog import load_catalog_sample
from src.horizon import parse_horizon_file
from src.planner import Planner

st.set_page_config(page_title="Astro Planner", layout="wide")
st.title("Astro Planner — Night Target Selector (MVP)")

# Load example config and catalogs
config_path = "config/config.toml"
config = load_config_example(config_path)

st.sidebar.header("Configuration")
st.sidebar.write(
    f"Location: {config['location']['latitude']}, "
    f"{config['location']['longitude']} "
    f"(elev {config['location']['elevation_m']} m)"
)

optics_names = [o['name'] for o in config['optics']]
optics_choice = st.sidebar.selectbox("Select optics/camera:", optics_names)

# ---------------------------
# Filters
# ---------------------------
st.sidebar.header("Filters")

max_mag = st.sidebar.number_input(
    "Max magnitude",
    value=float(12.0),   # make explicitly float
    step=0.5             # float step
)

min_alt = st.sidebar.number_input(
    "Min altitude (deg)",
    value=float(config['minima']['min_altitude_deg']),  # MUST be float
    step=1.0                                              # float step
)

# Streamlit Sliders must receive floats too if min/max/current are mixed types
fov_min = st.sidebar.slider(
    "Min FOV fill",
    float(0.0),                      # ensure float
    float(1.0),                      # ensure float
    (
        float(config['minima']['min_fov_fill']),
        float(config['minima']['max_fov_fill'])
    )
)

# ---------------------------
# Date / Moon / Weather
# ---------------------------
st.sidebar.header("Date & Moon")
date = st.sidebar.date_input("Date to plan (local)")
st.sidebar.checkbox("Show Clear Sky Clock (if within 48h)", value=False)

# ---------------------------
# Load catalogs & horizon
# ---------------------------
catalog = load_catalog_sample("data/messier_sample.csv")
horizon = parse_horizon_file("horizon/horizon_sample.txt")

# ---------------------------
# Planning logic
# ---------------------------
planner = Planner(config, catalog, horizon)

results = planner.plan(
    date,
    optics_choice,
    max_magnitude=max_mag,
    min_altitude=min_alt,
    fov_fill_range=fov_min
)

# ---------------------------
# Output
# ---------------------------
st.subheader("Candidate Targets")

if results is not None and not results.empty:
    st.write(f"Found {len(results)} candidates for {date}")
    st.table(results)
else:
    st.write("No candidates found with current filters.")

st.info(
    "This is an MVP skeleton. The src modules contain the structure and "
    "placeholder functions.\nNext steps: implement precise ephemeris "
    "calculations using Skyfield, visibility sampling, scoring, and "
    "Streamlit UI refinements."
)

