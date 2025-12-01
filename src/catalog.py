import pandas as pd

def load_catalog_sample(path):
    """Load a simple CSV catalog into a pandas DataFrame.
    Expected columns: name, ra_hms, dec_dms, type, major_arcmin, minor_arcmin, magnitude, notes
    """
    df = pd.read_csv(path)
    return df
