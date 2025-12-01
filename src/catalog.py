# src/catalog.py

import csv

def load_catalog_sample(csv_path):
    """
    Load a sample catalog CSV file and return a list of dicts.

    Expected CSV columns:
      - name: object name (string)
      - mag: magnitude (float)
      - size_deg: angular size in degrees (float)
      - type: optional, e.g., galaxy, nebula

    Returns:
      list of dicts, e.g.
      [
          {"name": "M1", "mag": 8.4, "size_deg": 0.33, "type": "SNR", "altitude_deg": 0.0},
          ...
      ]
    """
    catalog = []
    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            obj = {
                "name": row.get("name", "Unknown"),
                "mag": float(row.get("mag", 99.0)),
                "size_deg": float(row.get("size_deg", 1.0)),
                "type": row.get("type", "Unknown"),
                "altitude_deg": 0.0  # placeholder; will be computed later
            }
            catalog.append(obj)
    return catalog
