//
// Third-Octave Bandpass Filter Bank Implementation
//

#pragma once

#include <stdio.h>

static enum ChannelType {
    LEFT_CHANNEL = 0,
    RIGHT_CHANNEL = 1,
    STEREO_CHANNEL = 2
};

// Function declarations
void third_octave_filter_alloc(float sample_rate, int filter_order);
void third_octave_filter_free(void);
void third_octave_filter_process(float *data, int buffer_size);
float third_octave_filter_getValue(int param);
void third_octave_filter_setValue(int param, float val);