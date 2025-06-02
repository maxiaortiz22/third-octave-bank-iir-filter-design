def compute_band_edges(center_freqs):
    """Compute band edges for third-octave bands"""
    bands = []
    for f_c in center_freqs:
        f_lower = f_c * (2 ** (-1/6))  # Lower bandedge
        f_upper = f_c * (2 ** (1/6))   # Upper bandedge
        bands.append((f_c, f_lower, f_upper))
    return bands