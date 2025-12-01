"""
src/ephemeris.py

Utility functions to load and cache a skyfield ephemeris (e.g., de421.bsp).
This module provides a thin wrapper around skyfield.api.Loader to:
 - use a local cache directory (".skyfield_data" by default),
 - attempt to load a local file if present,
 - fall back to downloading the ephemeris if necessary,
 - expose the timescale and ephemeris objects for downstream modules.

Design goals:
 - Fail-fast and informative if Skyfield is not installed.
 - Do not perform large downloads at import time; provide functions that the app calls explicitly.
"""

from pathlib import Path
import os

# Try to import Skyfield; if not present, set a flag and keep error for later.
try:
    from skyfield.api import Loader
    SKYFIELD_AVAILABLE = True
except Exception as _e:
    SKYFIELD_AVAILABLE = False
    _SKYFIELD_IMPORT_ERROR = _e

_DEFAULT_CACHE_DIR = ".skyfield_data"
_DEFAULT_EPHEMERIS = "de421.bsp"

from skyfield.api import Star


def _ensure_cache_dir(cache_dir):
    p = Path(cache_dir)
    p.mkdir(parents=True, exist_ok=True)
    return str(p)

def load_ephemeris(ephem_name=_DEFAULT_EPHEMERIS, cache_dir=_DEFAULT_CACHE_DIR, download_if_missing=True):
    """
    Return (timescale, ephemeris) tuple using Skyfield Loader.
    - ephem_name: filename of ephemeris (e.g., 'de421.bsp' or 'de440s.bsp')
    - cache_dir: local directory to cache ephemeris & related data
    - download_if_missing: if True, attempt to download the file if not present locally

    Raises:
      ImportError if Skyfield is not installed.
      FileNotFoundError if the ephemeris is not found and download_if_missing is False.
      Any exception raised by skyfield Loader when downloading/loading.
    """
    if not SKYFIELD_AVAILABLE:
        raise ImportError("Skyfield is required for ephemeris loading but is not installed. "
                          f"Original import error: {_SKYFIELD_IMPORT_ERROR!r}")

    cache_dir = _ensure_cache_dir(cache_dir)
    loader = Loader(cache_dir)

    # First try to load the ephemeris from the local cache (Loader handles caching).
    try:
        ts = loader.timescale()
        # Loader.__call__ will return a path or download as needed.
        eph = loader(ephem_name)
        return ts, eph
    except FileNotFoundError:
        if not download_if_missing:
            raise
        # Try to download using loader (Loader will attempt to fetch if not present).
        eph = loader(ephem_name)
        ts = loader.timescale()
        return ts, eph

def get_default_ephemeris(cache_dir=_DEFAULT_CACHE_DIR, download_if_missing=True):
    """
    Convenience function that returns (ts, eph) for the default ephemeris file.
    """
    return load_ephemeris(_DEFAULT_EPHEMERIS, cache_dir=cache_dir, download_if_missing=download_if_missing)

def ephemeris_available(cache_dir=_DEFAULT_CACHE_DIR):
    """
    Return True if Skyfield is installed and the default ephemeris file exists in cache_dir.
    Does NOT attempt to download the file.
    """
    if not SKYFIELD_AVAILABLE:
        return False
    p = Path(cache_dir) / _DEFAULT_EPHEMERIS
    return p.exists()

def get_body_radec(ephem, name, ts, t):
    """
    Returns (RA(deg), Dec(deg)) for a catalog object.
    `name` can be:
      - a string matching a solar-system body (e.g., 'sun', 'moon')
      - a dict with 'ra' and 'dec' fields (fixed coordinates)
    """
    # Handle custom RA/Dec objects first
    if isinstance(name, dict) and 'ra' in name and 'dec' in name:
        return float(name['ra']), float(name['dec'])

    # Handle known solar system bodies
    if isinstance(name, str):
        if name.lower() in ephem:
            astrometric = ephem[name.lower()].at(t)
            ra, dec, _ = astrometric.radec()
            return float(ra.hours) * 15.0, float(dec.degrees)

    raise ValueError(f"Unknown target name or missing RA/Dec: {name}")

if __name__ == "__main__":
    # Quick smoke test when run as a script: do not attempt to download automatically.
    print("Skyfield available:", SKYFIELD_AVAILABLE)
    if SKYFIELD_AVAILABLE:
        from pprint import pprint
        cache = _ensure_cache_dir(_DEFAULT_CACHE_DIR)
        print("Cache dir:", cache)
        print("Default ephemeris present:", ephemeris_available(cache))
        print("To load the ephemeris, call get_default_ephemeris(cache_dir, download_if_missing=True)")
    else:
        print("Skyfield import error:", _SKYFIELD_IMPORT_ERROR)
