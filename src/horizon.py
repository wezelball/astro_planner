import math

def parse_horizon_file(path):
    """Parse a horizon file with lines 'az alt' OR 'az,alt' (degrees).
    Returns a sorted list of (az, alt) tuples.
    """
    pts = []
    with open(path, 'r') as f:
        for line in f:
            line = line.split('#', 1)[0].strip()
            if not line:
                continue

            # Accept either comma OR whitespace separation
            if ',' in line:
                parts = line.split(',', 1)
            else:
                parts = line.split()

            if len(parts) != 2:
                continue

            try:
                az = float(parts[0])
                alt = float(parts[1])
                pts.append((az % 360.0, alt))
            except ValueError:
                continue

    pts.sort()

    # ensure wrap-around
    if pts and pts[0][0] != 0.0:
        pts = [(0.0, pts[0][1])] + pts

    return pts


def horizon_altitude_at(pts, az):
    """Linear interpolation of horizon altitude at given az (degrees)."""
    az = az % 360.0
    if not pts:
        return 0.0
    for i in range(len(pts)-1):
        a0, h0 = pts[i]
        a1, h1 = pts[i+1]
        if a0 <= az <= a1:
            # interpolate
            if a1 == a0:
                return h0
            t = (az - a0) / (a1 - a0)
            return h0 + t*(h1 - h0)
    # wrap-around between last and first
    a0, h0 = pts[-1]
    a1, h1 = pts[0][0]+360.0, pts[0][1]
    t = (az - a0) / (a1 - a0)
    return h0 + t*(h1 - h0)
