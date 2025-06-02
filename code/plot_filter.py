import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

def plot_filter_responses(filters, bands, sr, title):
    """Plot frequency responses of cascaded biquad filters for third-octave bands."""
    
    plt.figure(figsize=(12, 6))
    
    for band_idx, biquads in enumerate(filters):
        if not biquads:  # Skip if filter design failed
            continue
            
        # Compute overall frequency response by cascading biquads
        w = np.logspace(np.log10(15), np.log10(22000), 2048)
        H_total = np.ones(len(w), dtype=complex)
        
        for biquad in biquads:
            b, a = biquad['b'], biquad['a']
            _, h = signal.freqz(b, a, worN=w, fs=sr)
            H_total *= h
        
        # Plot with selective labeling
        label = f"{bands[band_idx][0]:.1f} Hz" if band_idx % 4 == 0 else ""
        plt.semilogx(w, 20 * np.log10(np.abs(H_total)), label=label)
    
    plt.title(title)
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Magnitude (dB)")
    plt.grid(True, which='both', ls='--', alpha=0.7)
    plt.xticks([20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000],
               [r'$20$', r'$50$', r'$100$', r'$200$', r'$500$', r'$1k$', r'$2k$', r'$5k$', r'$10k$', r'$20k$'],
               rotation=45) 
    plt.xlim(15, 22000)
    plt.ylim(-60, 5)
    if any(band_idx % 4 == 0 for band_idx in range(len(filters))):
        plt.legend(loc='lower left', fontsize='small', ncol=2)
    plt.tight_layout()