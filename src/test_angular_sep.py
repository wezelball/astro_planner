#src/test_angular_sep.py

import numpy as np

def ra_to_degrees(ra_hms):
    """
    Convert RA from hours:minutes:seconds (or decimal hours) to degrees.
    
    Parameters:
    ra_hms: str 'H:M:S', tuple/list (H,M,S), or float decimal hours
    
    Returns:
    RA in degrees (0-360) [web:58][web:59].
    """
    if isinstance(ra_hms, str):
        h, m, s = map(float, ra_hms.split(':'))
    elif isinstance(ra_hms, (list, tuple)):
        h, m, s = ra_hms
    else:  # decimal hours
        return np.radians(ra_hms * 15)
    
    # Convert to decimal hours
    decimal_hours = h + m/60 + s/3600
    return decimal_hours * 15  # 15 deg/hour

def ra_degrees_vectorized(ra_hms_array):
    """
    Vectorized conversion for array of RA inputs.
    """
    degrees = np.zeros(len(ra_hms_array))
    for i, ra in enumerate(ra_hms_array):
        degrees[i] = ra_to_degrees(ra)
    return degrees

def angular_separation_vectorized(ra1, dec1, ra2, dec2):
    """
    Vectorized angular separation for arrays of RA/Dec coordinates.
    
    Parameters:
    ra1, dec1: Arrays of first object's coordinates (degrees), shape (N,)
    ra2, dec2: Arrays of second object's coordinates (degrees), shape (M,)
    
    Returns:
    2D array of separations (degrees), shape (N, M) [web:38][web:48].
    """
    # Convert to radians with broadcasting
    ra1_rad = np.radians(ra1[:, np.newaxis])   # Shape (N, 1)
    dec1_rad = np.radians(dec1[:, np.newaxis]) # Shape (N, 1)
    ra2_rad = np.radians(ra2)                  # Shape (M,)
    dec2_rad = np.radians(dec2)                # Shape (M,)
    
    # Broadcasting creates (N, M) arrays automatically
    cos_theta = (np.sin(dec1_rad) * np.sin(dec2_rad) + 
                 np.cos(dec1_rad) * np.cos(dec2_rad) * 
                 np.cos(ra1_rad - ra2_rad))
    
    theta = np.degrees(np.arccos(np.clip(cos_theta, -1.0, 1.0)))
    return theta


# Vectorized for Messier objects
messier_ras = ['00:42:44', '08:40:22', '16:41:41', '13:42:11', '01:35:20', '07:13:12.55']  # M31, M44, M13, M3, M33, Moon
messier_deg = ra_degrees_vectorized(messier_ras)
print("\nMessier RA in degrees:", messier_deg)


# Example: Moon separation from Messier objects
messier_ras = np.array([10.53, 130.1, 250.4, 205.6, 23.8])  # M31, M44, M13, M3, M33 (degrees)
messier_decs = np.array([41.3, 19.7, 28.4, 26.3, 30.8])
moon_ra, moon_dec = np.array([108.3]), np.array([26.3])  # Single Moon position

separations = angular_separation_vectorized(messier_ras, messier_decs, 
                                          moon_ra, moon_dec)
print("Moon separations from Messier objects (degrees):")
for name, sep in zip(["M31 (NGC224)", "M44 (NGC2632)", "M13 (NGC6205)", "M3 (NGC5272)", "M33 (NGC598)"], separations.flatten()):
    print(f"{name}: {sep:.2f}°")