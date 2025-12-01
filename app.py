import streamlit as st
from src.config_loader import load_config_example
from src.catalog import load_catalog_sample
from src.horizon import parse_horizon_file
from src.planner import Planner
from src.equipment import Telescope, Camera, Optics

st.set_page_config(page_title="Astro Planner", layout="wide")
st.title("Astro Planner — Night Target Selector (MVP)")

# ----------------------------
# Load config and catalogs
# ----------------------------
config_path = "config/config.toml"
config = load_config_example(config_path)

catalog = load_catalog_sample("data/messier_sample.csv")
horizon = parse_horizon_file("horizon/horizon_sample.txt")

# ----------------------------
# Sidebar: Configuration
# ----------------------------
st.sidebar.header("Configuration")
st.sidebar.write(
    f"Location: {config['location']['latitude']}, "
    f"{config['location']['longitude']} (elev {config['location']['elevation_m']} m)"
)

# ----------------------------
# Optics selection
# ----------------------------
optics_names = [o['name'] for o in config['optics']]
optics_choice_name = st.sidebar.selectbox("Select optics/camera:", optics_names)
optics_dict = next(o for o in config['optics'] if o['name'] == optics_choice_name)

# Build Telescope and Camera objects
telescope = Telescope(
    name=optics_dict['name'],
    aperture_mm=optics_dict['telescope_aperture_mm'],
    focal_length_mm=optics_dict['telescope_focal_length_mm']
)

camera = Camera(
    name=optics_dict['name'],
    sensor_width_mm=optics_dict['camera_sensor_width_mm'],
    sensor_height_mm=optics_dict['camera_sensor_height_mm'],
    pixel_width=optics_dict['camera_pixel_width'],
    pixel_height=optics_dict['camera_pixel_height']
)

# Build Optics object with optional focal reducer
optics = Optics(telescope, camera, focal_reducer=optics_dict.get('focal_reducer', 1.0))

# Display optics info
st.sidebar.write(f"Effective focal length (mm): {optics.effective_focal_length_mm:.1f}")
fov_x, fov_y = optics.fov_deg
st.sidebar.write(f"Field of View (deg): {fov_x:.2f} x {fov_y:.2f}")
scale_x, scale_y = optics.pixel_scale_arcsec
st.sidebar.write(f"Pixel scale (arcsec/pixel): {scale_x:.2f} x {scale_y:.2f}")

# ----------------------------
# Filters
# ----------------------------
st.sidebar.header("Filters")
max_mag = st.sidebar.number_input("Max magnitude", value=12.0, step=0.5)
min_alt = st.sidebar.number_input(
    "Min altitude (deg)",
    value=float(config['minima']['min_altitude_deg']),
    step=1.0
)

fov_min_fill = config['minima']['min_fov_fill']
fov_max_fill = config['minima']['max_fov_fill']
st.sidebar.write(f"Automatically applying FOV fill range: {fov_min_fill:.2f} - {fov_max_fill:.2f}")

# ----------------------------
# Date & Moon
# ----------------------------
st.sidebar.header("Date & Moon")
date = st.sidebar.date_input("Date to plan (local)")
st.sidebar.checkbox("Show Clear Sky Clock (if within 48h)", value=False)

# ----------------------------
# Planner
# ----------------------------
planner = Planner(config, catalog, horizon)
results = planner.plan(
    date,
    optics_choice_name,
    max_magnitude=max_mag,
    min_altitude=min_alt,
    fov_fill_range=(fov_min_fill, fov_max_fill),
    optics=optics  # Pass Optics object for FOV fill calculation
)

# ----------------------------
# Display results
# ----------------------------
st.subheader("Candidate Targets")
if results is not None and not results.empty:
    st.write(f"Found {len(results)} candidates for {date}")
    st.table(results)
else:
    st.write("No candidates found with current filters.")

st.info(
    "This is an MVP skeleton. The src modules contain placeholder functions. "
    "Next steps: implement precise ephemeris calculations using Skyfield, "
    "visibility sampling, scoring, and Streamlit UI refinements."
)
