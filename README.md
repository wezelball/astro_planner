# Astro Planner - MVP

This is an initial skeleton for *Astro Planner* — a Streamlit-based web dashboard that helps pick astrophotography targets
based on your equipment, local horizon, catalogs (Messier/NGC/IC), and moon conditions.

This repository contains:
- `app.py` - Streamlit entrypoint (MVP UI)
- `config/config.toml` - example user configuration
- `data/messier_sample.csv` - small sample catalog entries
- `horizon/horizon_sample.txt` - sample Cartes du Ciel compatible horizon file
- `src/` - core modules (config_loader, catalog, horizon, planner)
- `requirements.txt` - suggested Python packages

Notes:
- The astro computation code uses `skyfield`/`astropy` in comments and is structured to be plug-and-play.
- You can run the app locally with Streamlit once dependencies are installed:
  ```bash
  pip install -r requirements.txt
  streamlit run app.py
  ```
