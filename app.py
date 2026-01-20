import sys
from streamlit.web import cli as stcli  # for recent Streamlit versions
import streamlit as st
from datetime import datetime, date, time, timedelta, timezone
from zoneinfo import ZoneInfo   # Python 3.9+
from src.config_loader import load_config
from src.catalog import load_openngc_catalog
from src.horizon import parse_horizon_file
from src.planner import Planner
from src.optics import Optics
from src.moon import (
    moon_phase_info,
    moon_position_topocentric,
)
import pandas as pd
from math import isnan
from src.plotting import plot_sky_polar
from skyfield.api import Loader, wgs84
from skyfield import almanac

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

# Format times to local timezone for display (if found)
def _format_tt_local(tt):
    if tt is None:
        return "N/A"
    # convert Skyfield Time to timezone-aware UTC datetime, then to LOCAL_TZ
    dt_utc = tt.utc_datetime().replace(tzinfo=timezone.utc)
    dt_local = dt_utc.astimezone(LOCAL_TZ)
    return dt_local.strftime("%Y-%m-%d %H:%M")    


def main():
    #import streamlit as st
    st.title("Wezelball's Astro Planner ")

    # -------------------------------------------------------
    # App setup
    # -------------------------------------------------------

    # Load once
    object_list = load_opengc_catalog()

    # -------------------------------------------------------
    # Debug: check how many OpenNGC objects have size data
    # -------------------------------------------------------
    num_objects = len(object_list)

    st.set_page_config(page_title="Astro Planner", layout="wide")
    st.title("Astro Planner — Night Target Selector (MVP)")

    # Load example config and catalogs
    config_path = "config/config.toml"
    config = load_config(config_path)

    st.sidebar.header("Configuration")

    # Select optics/camera
    optics_names = [o['name'] for o in config['optics']]
    optics_choice = st.sidebar.selectbox("Select optics/camera:", optics_names)

    # Topos (lat/lon) and geocentric observer
    lat = config["location"]["latitude"]
    lon = config["location"]["longitude"]
    elev = config["location"]["elevation_m"]
    load = Loader('./skyfield_data')
    ts = load.timescale()
    eph = load('de421.bsp')

    # Used for almanac functions
    topos = wgs84.latlon(lat, lon, elevation_m=elev)
    # used for alt/az, Moon, Sun
    observer = eph['earth'] + topos


    # -------------------------------------------------------
    # Sidebar: configuration & filters
    # -------------------------------------------------------

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

    # Define selectable types (matches Planner.filter logic)
    selectable_types = ["G", "OCl", "PN", "Nebula", "SNR", "GCl"]
    selected_type = st.sidebar.selectbox("Select object type", selectable_types)

    # Sidebar date & moon
    st.sidebar.header("Date/Time")

    date = st.sidebar.date_input("Date to plan (local)")

    # Select snapshot time (UTC)
    snapshot_time = st.sidebar.time_input(
        "Snapshot Time (local civil time)",
        value=time(hour=4, minute=0),  # default: 04:00 UTC ≈ local midnight
        key="snapshot_time_picker"
    )

    # Build LOCAL datetime from selected date + time
    local_dt = datetime.combine(date, snapshot_time).replace(tzinfo=LOCAL_TZ)

    # Convert to UTC
    utc_dt = local_dt.astimezone(ZoneInfo("UTC"))

    # Pass UTC date/time components to planner
    date_utc = utc_dt.date()
    hour_utc = utc_dt.hour
    minute_utc = utc_dt.minute

    # Use UTC date/time (utc_dt) consistently for the snapshot used for moon calculations
    t_snapshot = ts.utc(
        utc_dt.year, utc_dt.month, utc_dt.day,
        utc_dt.hour, utc_dt.minute
    )


    # Sidebar date & moon
    st.sidebar.header("Moon")

    # Moon separation slider
    moon_sep_min = st.sidebar.slider(
        "Minimum Moon Separation (°)",
        min_value=0,
        max_value=90,
        value=30,
        step=5
    )

    # ------------------------------------------------------
    # Moon Phase Sidebar Widget
    # ------------------------------------------------------

    illum_pct, waxing, phase_name, age_days = moon_phase_info(eph, ts, t_snapshot)
    waxing_text = "Waxing" if waxing else "Waning"
    moon_alt_deg, moon_az_deg, _ = moon_position_topocentric(eph, observer, t_snapshot)

    # ---------------------------
    #  Observer Group
    # ---------------------------
    with st.sidebar.expander("📍 Observer", expanded=False):
        col1, col2 = st.columns(2)
        col1.markdown(f"**Lat**<br><span style='font-size:14px'>{config['location']['latitude']:.4f}°</span>", unsafe_allow_html=True)
        col2.markdown(f"**Lon**<br><span style='font-size:14px'>{config['location']['longitude']:.4f}°</span>", unsafe_allow_html=True)

        st.markdown(f"**Elevation**<br><span style='font-size:14px'>{config['location']['elevation_m']} m</span>", unsafe_allow_html=True)

    # ---------------------------
    #  Moon Group
    # ---------------------------
    with st.sidebar.expander("🌙 Moon", expanded=True):

        # Phase %  Progress Bar
        st.write(f"**Phase:** {illum_pct:.1f}% illuminated")
        st.progress(illum_pct / 100.0)

        # Age + Illum
        col1, col2 = st.columns(2)
        col1.markdown(
            f"**Age**<br><span style='font-size:14px'>{age_days:.2f} d</span>",
            unsafe_allow_html=True,
        )
        col2.markdown(
            f"**Illum**<br><span style='font-size:14px'>{illum_pct:.1f}%</span>",
            unsafe_allow_html=True,
        )

        # Alt + Az
        col3, col4 = st.columns(2)
        col3.markdown(
            f"**Alt**<br><span style='font-size:14px'>{moon_alt_deg:.1f}°</span>",
            unsafe_allow_html=True,
        )
        col4.markdown(
            f"**Az**<br><span style='font-size:14px'>{moon_az_deg:.1f}°</span>",
            unsafe_allow_html=True,
        )

    # ---------------------------
    #  Sun Group
    # ---------------------------
    with st.sidebar.expander("☀️ Sun", expanded=True):

        # Topocentric sun position at snapshot
        sun = eph['sun']
        sun_astrom = observer.at(t_snapshot).observe(sun).apparent()
        sun_alt, sun_az, _ = sun_astrom.altaz()
        sun_alt_deg = float(sun_alt.degrees)
        sun_az_deg = float(sun_az.degrees)

        # Display compact alt/az
        col1, col2 = st.columns(2)
        col1.markdown(f"**Alt**<br><span style='font-size:14px'>{sun_alt_deg:.1f}°</span>", unsafe_allow_html=True)
        col2.markdown(f"**Az**<br><span style='font-size:14px'>{sun_az_deg:.1f}°</span>", unsafe_allow_html=True)

        # ----------------------------------------
        # Astronomical dusk / dawn (Sun at -18°)
        # ----------------------------------------
        f = almanac.dark_twilight_day(eph, topos)

        # Search window: from yesterday to +2 days (safe)
        t0 = ts.utc((utc_dt - timedelta(days=1)).year,
                    (utc_dt - timedelta(days=1)).month,
                    (utc_dt - timedelta(days=1)).day)
        t1 = ts.utc((utc_dt + timedelta(days=2)).year,
                    (utc_dt + timedelta(days=2)).month,
                    (utc_dt + timedelta(days=2)).day)

        times, states = almanac.find_discrete(t0, t1, f)

        next_astro_sunset_t = None   # entry into astronomical night
        next_astro_sunrise_t = None  # exit from astronomical night

        for tt, state in zip(times, states):
            if tt.tt <= t_snapshot.tt:
                continue

            # state == 0 means astronomical night
            if state == 0 and next_astro_sunset_t is None:
                next_astro_sunset_t = tt

            if state != 0 and next_astro_sunrise_t is None and next_astro_sunset_t is not None:
                next_astro_sunrise_t = tt

            if next_astro_sunset_t is not None and next_astro_sunrise_t is not None:
                break

        astro_sunset_str  = _format_tt_local(next_astro_sunset_t)
        astro_sunrise_str = _format_tt_local(next_astro_sunrise_t)

        col3, col4 = st.columns(2)
        col3.markdown(
            f"**Astro Dusk**<br><span style='font-size:13px'>{astro_sunset_str}</span>",
            unsafe_allow_html=True,
        )
        col4.markdown(
            f"**Astro Dawn**<br><span style='font-size:13px'>{astro_sunrise_str}</span>",
            unsafe_allow_html=True,
        )

        st.caption("Astronomical twilight (Sun at −18°), local time")
        

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

    # -------------------------------------------------------
    # Planning & results
    # -------------------------------------------------------

    # Instantiate planner
    planner = Planner(config, catalog_df, horizon)

    # Call the heart of the application, Planner.plan()
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

if __name__ == "__main__":
    # If already running under Streamlit, just call main()
    from streamlit import runtime
    if runtime.exists():
        main()
    else:
        # Reinvoke via "streamlit run app.py"
        sys.argv = ["streamlit", "run", sys.argv[0]]
        sys.exit(stcli.main())
