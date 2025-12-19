"""
config_loader.py
Configuration loading utilities.

This module provides a lightweight TOML configuration loader for the
astro planner application. It performs minimal validation and applies
reasonable defaults for missing configuration sections.

This is intentionally not a strict schema validator; it is designed
to keep the application runnable even with partial configuration files.
"""
import tomli

def load_config(path='config/config.toml'):
    """
    Load the application configuration from a TOML file.

    The configuration is parsed into a plain Python dictionary.
    Only minimal validation is performed; missing sections are
    populated with default values to ensure the application
    can run with incomplete configs.

    Default sections added if missing:

    minima:
        - min_altitude_deg (float): Minimum object altitude [deg]
        - min_fov_fill (float): Minimum field-of-view fill fraction
        - max_fov_fill (float): Maximum field-of-view fill fraction

    moon:
        - max_illumination_pct (float): Maximum allowed Moon illumination [%]
        - min_angular_separation_deg (float): Minimum Moon-object separation [deg]

    Args:
        path (str): Path to the TOML configuration file.

    Returns:
        dict: Configuration dictionary with defaults applied.
    """
    with open(path, 'rb') as f:
        cfg = tomli.load(f)
    # Provide some normalization/defaults
    if 'minima' not in cfg:
        cfg['minima'] = {'min_altitude_deg':20, 'min_fov_fill':0.15, 'max_fov_fill':0.85}
    if 'moon' not in cfg:
        cfg['moon'] = {'max_illumination_pct':30, 'min_angular_separation_deg':30}
    return cfg
