import pandas as pd
from astropy.coordinates import Angle
import numpy as np

def load_openngc_catalog(path: str) -> pd.DataFrame:
    """
    Load and convert the OpenNGC NGC.csv catalog into a format
    suitable for the Planner, including:
        - RA/Dec converted to decimal degrees
        - Size converted from arcminutes to degrees
        - Magnitude selection (V-Mag preferred)
    """
    # OpenNGC uses semicolon separator
    df = pd.read_csv(path, sep=';', dtype=str)

    # --- RA/Dec conversion ---
    def parse_ra(ra_str):
        try:
            return Angle(ra_str, unit="hourangle").degree
        except Exception:
            return np.nan

    def parse_dec(dec_str):
        try:
            return Angle(dec_str, unit="deg").degree
        except Exception:
            return np.nan

    df["ra_deg"] = df["RA"].apply(parse_ra)
    df["dec_deg"] = df["Dec"].apply(parse_dec)

    # --- Object size (MajAx & MinAx are arcminutes) ---
    def parse_size(val):
        try:
            return float(val) / 60.0  # arcmin → degrees
        except:
            return np.nan

    df["major_arcmin"] = df["MajAx"].apply(parse_size)
    df["minor_arcmin"] = df["MinAx"].apply(parse_size)

    # Use the larger axis as object "size"
    df["size_deg"] = df[["major_arcmin", "minor_arcmin"]].max(axis=1)

    # --- Magnitude selection ---
    # Prefer V-Mag > B-Mag > SurfBr
    def choose_mag(row):
        for m in ["V-Mag", "B-Mag", "SurfBr"]:
            val = row.get(m)
            try:
                return float(val)
            except:
                pass
        return np.nan

    df["magnitude"] = df.apply(choose_mag, axis=1)

    # --- Clean and restrict fields ---
    cleaned = df[[
        "Name", "Type", "ra_deg", "dec_deg", "size_deg",
        "magnitude", "Const", "Common names"
    ]].rename(columns={
        "Name": "name",
        "Type": "type",
        "Const": "constellation",
        "Common names": "common_names",
    })

    # Drop objects without coordinates
    cleaned = cleaned.dropna(subset=["ra_deg", "dec_deg"])

    return cleaned
