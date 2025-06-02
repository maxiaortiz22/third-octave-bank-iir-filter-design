def generate_cpp_initialization(bands, filters_order_2, filters_order_4):
    """Generate C++ coefficient initialization code"""

    with open("data/cpp_coefficient_initialization.txt", "w") as f:
        print("// C++ Coefficient Initialization for Third-Octave Filters", file=f)
        print("void initializeCoefficients() {", file=f)
        print("    // Initialize all bands with their center frequencies", file=f)
        print("    for (int band = 0; band < NUM_BANDS; band++) {", file=f)
        print("        pFilterBank->bands[band].center_freq = center_frequencies[band];", file=f)
        print("    }", file=f)
        print("", file=f)
        
        # Order 2 coefficients
        print("    if (pFilterBank->filter_order == 2) {", file=f)
        for band_idx, biquads in enumerate(filters_order_2):
            if not biquads:
                continue
            center_freq = bands[band_idx][0]
            print(f"        // Band {band_idx}: {center_freq} Hz", file=f)
            print(f"        pFilterBank->bands[{band_idx}].num_sections = {len(biquads)};", file=f)
            
            for sec_idx, biquad in enumerate(biquads):
                print(f"        pFilterBank->bands[{band_idx}].sections[{sec_idx}].b[0] = {biquad['b'][0]:.12e};", file=f)
                print(f"        pFilterBank->bands[{band_idx}].sections[{sec_idx}].b[1] = {biquad['b'][1]:.12e};", file=f)
                print(f"        pFilterBank->bands[{band_idx}].sections[{sec_idx}].b[2] = {biquad['b'][2]:.12e};", file=f)
                print(f"        pFilterBank->bands[{band_idx}].sections[{sec_idx}].a[0] = {biquad['a'][0]:.12e};", file=f)
                print(f"        pFilterBank->bands[{band_idx}].sections[{sec_idx}].a[1] = {biquad['a'][1]:.12e};", file=f)
                print(f"        pFilterBank->bands[{band_idx}].sections[{sec_idx}].a[2] = {biquad['a'][2]:.12e};", file=f)
                print(f"        pFilterBank->bands[{band_idx}].sections[{sec_idx}].z[0] = 0.0;", file=f)
                print(f"        pFilterBank->bands[{band_idx}].sections[{sec_idx}].z[1] = 0.0;", file=f)
            print("", file=f)
        
        # Order 4 coefficients
        print("    } else { // filter_order == 4", file=f)
        for band_idx, biquads in enumerate(filters_order_4):
            if not biquads:
                continue
            center_freq = bands[band_idx][0]
            print(f"        // Band {band_idx}: {center_freq} Hz", file=f)
            print(f"        pFilterBank->bands[{band_idx}].num_sections = {len(biquads)};", file=f)
            
            for sec_idx, biquad in enumerate(biquads):
                print(f"        pFilterBank->bands[{band_idx}].sections[{sec_idx}].b[0] = {biquad['b'][0]:.12e};", file=f)
                print(f"        pFilterBank->bands[{band_idx}].sections[{sec_idx}].b[1] = {biquad['b'][1]:.12e};", file=f)
                print(f"        pFilterBank->bands[{band_idx}].sections[{sec_idx}].b[2] = {biquad['b'][2]:.12e};", file=f)
                print(f"        pFilterBank->bands[{band_idx}].sections[{sec_idx}].a[0] = {biquad['a'][0]:.12e};", file=f)
                print(f"        pFilterBank->bands[{band_idx}].sections[{sec_idx}].a[1] = {biquad['a'][1]:.12e};", file=f)
                print(f"        pFilterBank->bands[{band_idx}].sections[{sec_idx}].a[2] = {biquad['a'][2]:.12e};", file=f)
                print(f"        pFilterBank->bands[{band_idx}].sections[{sec_idx}].z[0] = 0.0;", file=f)
                print(f"        pFilterBank->bands[{band_idx}].sections[{sec_idx}].z[1] = 0.0;", file=f)
            print("", file=f)
        
        print("    }", file=f)
        print("}", file=f)