import numpy as np
from scipy import signal

def design_cascaded_filters(low_freq, high_freq, order, fs):
    """Design cascaded biquad filters for better numerical stability"""
    # Normalized frequencies with safety bounds
    low = max(low_freq / (fs / 2), 1e-6)
    high = min(high_freq / (fs / 2), 0.99)
    
    if order == 2:
        sos = signal.butter(2, [low, high], btype='bandpass', output='sos')
    else:  # order == 4
        sos = signal.butter(4, [low, high], btype='bandpass', output='sos')
    
    return sos

def normalize_sos_gain(sos, center_freq, fs):
    """Normalize SOS filter to have unity gain at center frequency"""
    w = 2 * np.pi * center_freq / fs
    
    # Calculate total gain through all sections
    total_gain = 1.0
    for section in sos:
        b = section[0:3]
        a = section[3:6]
        
        # Evaluate transfer function at center frequency
        z = np.exp(1j * w)
        H_num = b[0] + b[1] * z**(-1) + b[2] * z**(-2)
        H_den = a[0] + a[1] * z**(-1) + a[2] * z**(-2)
        H = H_num / H_den
        total_gain *= abs(H)
    
    # Normalize first section's b coefficients
    if total_gain > 0:
        sos[0, 0:3] /= total_gain
    
    return sos

def sos_to_cascaded_biquads(sos):
    """Convert SOS to format suitable for C++ implementation"""
    biquads = []
    for section in sos:
        biquad = {
            'b': section[0:3].copy(),  # b0, b1, b2
            'a': section[3:6].copy()   # a0, a1, a2 (a0 should be 1.0)
        }
        # Ensure a0 = 1.0
        if biquad['a'][0] != 1.0 and biquad['a'][0] != 0.0:
            biquad['b'] /= biquad['a'][0]
            biquad['a'] /= biquad['a'][0]
        biquads.append(biquad)
    
    return biquads