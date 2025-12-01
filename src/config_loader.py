import tomli

def load_config_example(path='config/config.toml'):
    """Load the TOML config file. Returns a simple dict.
    This function does minimal validation and returns defaults if keys missing.
    """
    with open(path, 'rb') as f:
        cfg = tomli.load(f)
    # Provide some normalization/defaults
    if 'minima' not in cfg:
        cfg['minima'] = {'min_altitude_deg':20, 'min_fov_fill':0.15, 'max_fov_fill':0.85}
    if 'moon' not in cfg:
        cfg['moon'] = {'max_illumination_pct':30, 'min_angular_separation_deg':30}
    return cfg
