"""
catalog.py
Catalog loading and normalization utilities.

This module provides helpers to load astronomical object catalogs
and normalize them into a format suitable for the Planner.

Currently supported:
- OpenNGC CSV catalog (NGC.csv)

All angular quantities are converted to decimal degrees.
"""

import pandas as pd
from astropy.coordinates import Angle
import numpy as np

def load_openngc_catalog(path: str) -> pd.DataFrame:
    """
    Load and normalize the OpenNGC NGC.csv catalog.

    The catalog is converted into a Planner-friendly format, including:
      - Right Ascension converted from HH:MM:SS to decimal degrees
      - Declination converted to decimal degrees
      - Object size converted from arcminutes to degrees
      - A single representative magnitude selected by priority

    Magnitude selection priority:
      1. V-Mag
      2. B-Mag
      3. Surface brightness (SurfBr)

    Args:
        path (str): Path to the OpenNGC CSV file (semicolon-delimited).

    Returns:
        pandas.DataFrame: Cleaned catalog with the following columns:
            - name (str)
            - type (str)
            - ra_deg (float)
            - dec_deg (float)
            - size_deg (float)
            - magnitude (float)
            - constellation (str)
            - common_names (str)

        Rows lacking valid RA/Dec are removed.
    """
    # OpenNGC uses semicolon separator
    df = pd.read_csv(path, sep=';', dtype=str)

    # --- RA/Dec conversion ---
    
    def parse_ra(ra_str):
        """
        Parse Right Ascension from sexagesimal string to degrees.

        Returns NaN on parse failure.
        """
        try:
            return Angle(ra_str, unit="hourangle").degree
        except Exception:
            return np.nan

    def parse_dec(dec_str):
        """
        Parse Declination from sexagesimal string to degrees.

        Returns NaN on parse failure.
        """
        try:
            return Angle(dec_str, unit="deg").degree
        except Exception:
            return np.nan

    df["ra_deg"] = df["RA"].apply(parse_ra)
    df["dec_deg"] = df["Dec"].apply(parse_dec)

    # --- Object size (MajAx & MinAx are arcminutes) ---
    def parse_size(val):
        """
        Convert size from arcminutes to degrees.

        Returns NaN if value is missing or invalid.
        """
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
        """
        Select the best available magnitude for an object.

        Priority order:
            V-Mag → B-Mag → SurfBr

        Returns NaN if none are valid.
        """
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
