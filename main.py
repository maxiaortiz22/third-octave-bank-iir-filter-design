import matplotlib.pyplot as plt
from code.frequency_cut import compute_band_edges
from code.filter_design import design_cascaded_filters, normalize_sos_gain, sos_to_cascaded_biquads
from code.plot_filter import plot_filter_responses
from code.export_to_cpp import generate_cpp_initialization

if __name__ == "__main__":
    # Sample rate (Hz)
    SAMPLE_RATE = 48000.0

    # Nominal center frequencies (Hz)
    NOMINAL_CENTER_FREQS = [
        20.0, 25.0, 31.5, 40.0, 50.0, 63.0, 80.0, 100.0, 125.0, 160.0,
        200.0, 250.0, 315.0, 400.0, 500.0, 630.0, 800.0, 1000.0, 1250.0,
        1600.0, 2000.0, 2500.0, 3150.0, 4000.0, 5000.0, 6300.0, 8000.0,
        10000.0, 12500.0, 16000.0, 20000.0
    ]

    # Compute band edges
    bands = compute_band_edges(NOMINAL_CENTER_FREQS)
    NUM_BANDS = len(bands)

    # Calculate coefficients
    filters_order_2 = []
    filters_order_4 = []

    print(f"Designing {NUM_BANDS} third-octave filters...")
    for band_idx, (center_freq, low_freq, high_freq) in enumerate(bands):
        print(f"Band {band_idx}: {center_freq} Hz ({low_freq:.1f} - {high_freq:.1f} Hz)")
        
        # Order 2
        try:
            sos_2 = design_cascaded_filters(low_freq, high_freq, 2, SAMPLE_RATE)
            sos_2_norm = normalize_sos_gain(sos_2.copy(), center_freq, SAMPLE_RATE)
            biquads_2 = sos_to_cascaded_biquads(sos_2_norm)
            filters_order_2.append(biquads_2)
        except Exception as e:
            print(f"Error designing order 2 filter for {center_freq} Hz: {e}")
            filters_order_2.append([])
        
        # Order 4
        try:
            sos_4 = design_cascaded_filters(low_freq, high_freq, 4, SAMPLE_RATE)
            sos_4_norm = normalize_sos_gain(sos_4.copy(), center_freq, SAMPLE_RATE)
            biquads_4 = sos_to_cascaded_biquads(sos_4_norm)
            filters_order_4.append(biquads_4)
        except Exception as e:
            print(f"Error designing order 4 filter for {center_freq} Hz: {e}")
            filters_order_4.append([])

    # Plot responses
    plot_filter_responses(filters_order_2, bands, SAMPLE_RATE, "Order 2 Third-Octave Filters (Cascaded Biquads)")
    plt.show()

    plot_filter_responses(filters_order_4, bands, SAMPLE_RATE, "Order 4 Third-Octave Filters (Cascaded Biquads)")
    plt.show()

    generate_cpp_initialization(bands, filters_order_2, filters_order_4)
    print(f"\nC++ coefficient initialization code generated in 'cpp_coefficient_initialization.txt'")
    print("Copy and paste this into your C++ file to replace the initializeCoefficients() function.")

    # Verification: Print some sample coefficients
    print("\nSample coefficients for verification:")
    print("Order 2, Band 0 (20 Hz):")
    if filters_order_2[0]:
        biquad = filters_order_2[0][0]
        print(f"  b: [{biquad['b'][0]:.6e}, {biquad['b'][1]:.6e}, {biquad['b'][2]:.6e}]")
        print(f"  a: [{biquad['a'][0]:.6e}, {biquad['a'][1]:.6e}, {biquad['a'][2]:.6e}]")

    print("Order 4, Band 15 (1000 Hz):")
    if filters_order_4[15]:
        for i, biquad in enumerate(filters_order_4[15]):
            print(f"  Section {i+1}:")
            print(f"    b: [{biquad['b'][0]:.6e}, {biquad['b'][1]:.6e}, {biquad['b'][2]:.6e}]")
            print(f"    a: [{biquad['a'][0]:.6e}, {biquad['a'][1]:.6e}, {biquad['a'][2]:.6e}]")