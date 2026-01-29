# src/astro_planner/paths.py

from pathlib import Path

def get_project_root():
    """Heuristic to find project root: src/astro_planner → astro_planner/"""
    current_dir = Path(__file__).parent
    src_dir = current_dir.parent
    project_root = src_dir.parent
    return project_root

ROOT_DIR = get_project_root()

# User data / config: live in your workdir
DATA_DIR = ROOT_DIR / "data"              # private NGC.csv etc.
CONFIG_DIR = ROOT_DIR / "config"          # private config.toml
HORIZON_DIR = ROOT_DIR / "horizon"        # private horizon file
SKYFIELD_DIR = ROOT_DIR / "skyfield_data" # local de421.bsp

if __name__ == "__main__":
    print("Package directory:", Path(__file__).parent)
    print("Project root     :", ROOT_DIR)
    print("Data dir         :", DATA_DIR)
    print("Config dir       :", CONFIG_DIR)
    print("SKYFIELD dir     :", SKYFIELD_DIR)
    print("Data files:", list(DATA_DIR.glob("*.csv")))
