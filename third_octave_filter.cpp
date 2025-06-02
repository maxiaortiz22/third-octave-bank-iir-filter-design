#include "third_octave_filter.h"
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Filter bank constants
#define NUM_BANDS 31 // All center frequencies for third-octave bands
#define MAX_BIQUAD_SECTIONS 2 // Max 2 sections for 4th order (2 biquads)

// Third-octave filter bank using cascaded biquad sections
// This approach provides better numerical stability than high-order transfer functions
// Sample rate: 48000.0 Hz
// Number of bands: 31

struct BiquadSection {
    double b[3];  // Numerator coefficients
    double a[3];  // Denominator coefficients (a[0] = 1.0)
    double z[2];  // Delay line (initialize to 0)
};

struct FilterBand {
    BiquadSection sections[2];  // Max 2 sections for 4th order
    int num_sections;
    double center_freq;
};

struct third_octave_filter {
    float bypass;
    float sample_rate;
    int filter_order; // 2 or 4
    int initialized;
    int channelType;
    struct FilterBand bands[NUM_BANDS]; // One filter per band
    float mic_constant; // Mic constant for calibration
    const float alpha = 0.99f; // Smoothing factor (adjustable, 0 < alpha < 1)
    int integrationTime;
    int samplesCount;
    float temporal_sum[NUM_BANDS]; // Temporal sum for each band
    float volume_level[NUM_BANDS]; // Volume level for each band
    float smoothed_level[NUM_BANDS]; // Smoothed level for each band
    uint64_t maxNumberOfSamples;
};

// Center frequencies (Hz) for third-octave bands
static const double center_frequencies[NUM_BANDS] = {
        20.0, 25.0, 31.5, 40.0, 50.0, 63.0, 80.0, 100.0, 125.0, 160.0,
        200.0, 250.0, 315.0, 400.0, 500.0, 630.0, 800.0, 1000.0, 1250.0,
        1600.0, 2000.0, 2500.0, 3150.0, 4000.0, 5000.0, 6300.0, 8000.0,
        10000.0, 12500.0, 16000.0, 20000.0
};

static struct third_octave_filter *pFilterBank;

// Process a single sample through one biquad section (Direct Form I)
double processBiquad(struct BiquadSection* section, double input) {
    // Direct Form I: y[n] = b0*x[n] + b1*x[n-1] + b2*x[n-2] - a1*y[n-1] - a2*y[n-2]
    double output = section->b[0] * input + section->z[0];
    section->z[0] = section->b[1] * input - section->a[1] * output + section->z[1];
    section->z[1] = section->b[2] * input - section->a[2] * output;
    return output;
}

// Process a sample through an entire filter band (cascade of biquads)
double processFilterBand(struct FilterBand* band, double input) {
    double output = input;
    for (int i = 0; i < band->num_sections; i++) {
        output = processBiquad(&band->sections[i], output);
    }
    return output;
}

// Process filtered samples and obtain average for each band
void initializeFilterBank(FilterBand* bands, int filter_order) {
    if (filter_order == 2) {
        // Band 0: 20.0 Hz
        bands[0].center_freq = 20.0;
        bands[0].num_sections = 2;
        bands[0].sections[0].b[0] = 9.1839749966e-08;
        bands[0].sections[0].b[1] = 1.8367949993e-07;
        bands[0].sections[0].b[2] = 9.1839749966e-08;
        bands[0].sections[0].a[0] = 1.0000000000e+00;
        bands[0].sections[0].a[1] = -1.9995282685e+00;
        bands[0].sections[0].a[2] = 9.9953634283e-01;
        bands[0].sections[0].z[0] = 0.0;
        bands[0].sections[0].z[1] = 0.0;
        bands[0].sections[1].b[0] = 1.0000000000e+00;
        bands[0].sections[1].b[1] = -2.0000000000e+00;
        bands[0].sections[1].b[2] = 1.0000000000e+00;
        bands[0].sections[1].a[0] = 1.0000000000e+00;
        bands[0].sections[1].a[1] = -1.9996006861e+00;
        bands[0].sections[1].a[2] = 9.9960650149e-01;
        bands[0].sections[1].z[0] = 0.0;
        bands[0].sections[1].z[1] = 0.0;

        // Band 1: 25.0 Hz
        bands[1].center_freq = 25.0;
        bands[1].num_sections = 2;
        bands[1].sections[0].b[0] = 1.4348423658e-07;
        bands[1].sections[0].b[1] = 2.8696847317e-07;
        bands[1].sections[0].b[2] = 1.4348423658e-07;
        bands[1].sections[0].a[0] = 1.0000000000e+00;
        bands[1].sections[0].a[1] = -1.9994078468e+00;
        bands[1].sections[0].a[2] = 9.9942046219e-01;
        bands[1].sections[0].z[0] = 0.0;
        bands[1].sections[0].z[1] = 0.0;
        bands[1].sections[1].b[0] = 1.0000000000e+00;
        bands[1].sections[1].b[1] = -2.0000000000e+00;
        bands[1].sections[1].b[2] = 1.0000000000e+00;
        bands[1].sections[1].a[0] = 1.0000000000e+00;
        bands[1].sections[1].a[1] = -1.9994990648e+00;
        bands[1].sections[1].a[2] = 9.9950815101e-01;
        bands[1].sections[1].z[0] = 0.0;
        bands[1].sections[1].z[1] = 0.0;

        // Band 2: 31.5 Hz
        bands[2].center_freq = 31.5;
        bands[2].num_sections = 2;
        bands[2].sections[0].b[0] = 2.2776385304e-07;
        bands[2].sections[0].b[1] = 4.5552770607e-07;
        bands[2].sections[0].b[2] = 2.2776385304e-07;
        bands[2].sections[0].a[0] = 1.0000000000e+00;
        bands[2].sections[0].a[1] = -1.9992498108e+00;
        bands[2].sections[0].a[2] = 9.9926983749e-01;
        bands[2].sections[0].z[0] = 0.0;
        bands[2].sections[0].z[1] = 0.0;
        bands[2].sections[1].b[0] = 1.0000000000e+00;
        bands[2].sections[1].b[1] = -2.0000000000e+00;
        bands[2].sections[1].b[2] = 1.0000000000e+00;
        bands[2].sections[1].a[0] = 1.0000000000e+00;
        bands[2].sections[1].a[1] = -1.9993658855e+00;
        bands[2].sections[1].a[2] = 9.9938030979e-01;
        bands[2].sections[1].z[0] = 0.0;
        bands[2].sections[1].z[1] = 0.0;

        // Band 3: 40.0 Hz
        bands[3].center_freq = 40.0;
        bands[3].num_sections = 2;
        bands[3].sections[0].b[0] = 3.6720162518e-07;
        bands[3].sections[0].b[1] = 7.3440325036e-07;
        bands[3].sections[0].b[2] = 3.6720162518e-07;
        bands[3].sections[0].a[0] = 1.0000000000e+00;
        bands[3].sections[0].a[1] = -1.9990406114e+00;
        bands[3].sections[0].a[2] = 9.9907290111e-01;
        bands[3].sections[0].z[0] = 0.0;
        bands[3].sections[0].z[1] = 0.0;
        bands[3].sections[1].b[0] = 1.0000000000e+00;
        bands[3].sections[1].b[1] = -2.0000000000e+00;
        bands[3].sections[1].b[2] = 1.0000000000e+00;
        bands[3].sections[1].a[0] = 1.0000000000e+00;
        bands[3].sections[1].a[1] = -1.9991899002e+00;
        bands[3].sections[1].a[2] = 9.9921315736e-01;
        bands[3].sections[1].z[0] = 0.0;
        bands[3].sections[1].z[1] = 0.0;

        // Band 4: 50.0 Hz
        bands[4].center_freq = 50.0;
        bands[4].num_sections = 2;
        bands[4].sections[0].b[0] = 5.7362965624e-07;
        bands[4].sections[0].b[1] = 1.1472593125e-06;
        bands[4].sections[0].b[2] = 5.7362965624e-07;
        bands[4].sections[0].a[0] = 1.0000000000e+00;
        bands[4].sections[0].a[1] = -1.9987908144e+00;
        bands[4].sections[0].a[2] = 9.9884126116e-01;
        bands[4].sections[0].z[0] = 0.0;
        bands[4].sections[0].z[1] = 0.0;
        bands[4].sections[1].b[0] = 1.0000000000e+00;
        bands[4].sections[1].b[1] = -2.0000000000e+00;
        bands[4].sections[1].b[2] = 1.0000000000e+00;
        bands[4].sections[1].a[0] = 1.0000000000e+00;
        bands[4].sections[1].a[1] = -1.9989802073e+00;
        bands[4].sections[1].a[2] = 9.9901654301e-01;
        bands[4].sections[1].z[0] = 0.0;
        bands[4].sections[1].z[1] = 0.0;

        // Band 5: 63.0 Hz
        bands[5].center_freq = 63.0;
        bands[5].num_sections = 2;
        bands[5].sections[0].b[0] = 9.1044093041e-07;
        bands[5].sections[0].b[1] = 1.8208818608e-06;
        bands[5].sections[0].b[2] = 9.1044093041e-07;
        bands[5].sections[0].a[0] = 1.0000000000e+00;
        bands[5].sections[0].a[1] = -1.9984601330e+00;
        bands[5].sections[0].a[2] = 9.9854020996e-01;
        bands[5].sections[0].z[0] = 0.0;
        bands[5].sections[0].z[1] = 0.0;
        bands[5].sections[1].b[0] = 1.0000000000e+00;
        bands[5].sections[1].b[1] = -2.0000000000e+00;
        bands[5].sections[1].b[2] = 1.0000000000e+00;
        bands[5].sections[1].a[0] = 1.0000000000e+00;
        bands[5].sections[1].a[1] = -1.9987033227e+00;
        bands[5].sections[1].a[2] = 9.9876100175e-01;
        bands[5].sections[1].z[0] = 0.0;
        bands[5].sections[1].z[1] = 0.0;

        // Band 6: 80.0 Hz
        bands[6].center_freq = 80.0;
        bands[6].num_sections = 2;
        bands[6].sections[0].b[0] = 1.4675488516e-06;
        bands[6].sections[0].b[1] = 2.9350977032e-06;
        bands[6].sections[0].b[2] = 1.4675488516e-06;
        bands[6].sections[0].a[0] = 1.0000000000e+00;
        bands[6].sections[0].a[1] = -1.9980175675e+00;
        bands[6].sections[0].a[2] = 9.9814666549e-01;
        bands[6].sections[0].z[0] = 0.0;
        bands[6].sections[0].z[1] = 0.0;
        bands[6].sections[1].b[0] = 1.0000000000e+00;
        bands[6].sections[1].b[1] = -2.0000000000e+00;
        bands[6].sections[1].b[2] = 1.0000000000e+00;
        bands[6].sections[1].a[0] = 1.0000000000e+00;
        bands[6].sections[1].a[1] = -1.9983339386e+00;
        bands[6].sections[1].a[2] = 9.9842693007e-01;
        bands[6].sections[1].z[0] = 0.0;
        bands[6].sections[1].z[1] = 0.0;

        // Band 7: 100.0 Hz
        bands[7].center_freq = 100.0;
        bands[7].num_sections = 2;
        bands[7].sections[0].b[0] = 2.2920635944e-06;
        bands[7].sections[0].b[1] = 4.5841271889e-06;
        bands[7].sections[0].b[2] = 2.2920635944e-06;
        bands[7].sections[0].a[0] = 1.0000000000e+00;
        bands[7].sections[0].a[1] = -1.9974822047e+00;
        bands[7].sections[0].a[2] = 9.9768387234e-01;
        bands[7].sections[0].z[0] = 0.0;
        bands[7].sections[0].z[1] = 0.0;
        bands[7].sections[1].b[0] = 1.0000000000e+00;
        bands[7].sections[1].b[1] = -2.0000000000e+00;
        bands[7].sections[1].b[2] = 1.0000000000e+00;
        bands[7].sections[1].a[0] = 1.0000000000e+00;
        bands[7].sections[1].a[1] = -1.9978887759e+00;
        bands[7].sections[1].a[2] = 9.9803404586e-01;
        bands[7].sections[1].z[0] = 0.0;
        bands[7].sections[1].z[1] = 0.0;

        // Band 8: 125.0 Hz
        bands[8].center_freq = 125.0;
        bands[8].num_sections = 2;
        bands[8].sections[0].b[0] = 3.5794339408e-06;
        bands[8].sections[0].b[1] = 7.1588678816e-06;
        bands[8].sections[0].b[2] = 3.5794339408e-06;
        bands[8].sections[0].a[0] = 1.0000000000e+00;
        bands[8].sections[0].a[1] = -1.9967906744e+00;
        bands[8].sections[0].a[2] = 9.9710568599e-01;
        bands[8].sections[0].z[0] = 0.0;
        bands[8].sections[0].z[1] = 0.0;
        bands[8].sections[1].b[0] = 1.0000000000e+00;
        bands[8].sections[1].b[1] = -2.0000000000e+00;
        bands[8].sections[1].b[2] = 1.0000000000e+00;
        bands[8].sections[1].a[0] = 1.0000000000e+00;
        bands[8].sections[1].a[1] = -1.9973162276e+00;
        bands[8].sections[1].a[2] = 9.9754315463e-01;
        bands[8].sections[1].z[0] = 0.0;
        bands[8].sections[1].z[1] = 0.0;

        // Band 9: 160.0 Hz
        bands[9].center_freq = 160.0;
        bands[9].num_sections = 2;
        bands[9].sections[0].b[0] = 5.8601557463e-06;
        bands[9].sections[0].b[1] = 1.1720311493e-05;
        bands[9].sections[0].b[2] = 5.8601557463e-06;
        bands[9].sections[0].a[0] = 1.0000000000e+00;
        bands[9].sections[0].a[1] = -1.9957808984e+00;
        bands[9].sections[0].a[2] = 9.9629679590e-01;
        bands[9].sections[0].z[0] = 0.0;
        bands[9].sections[0].z[1] = 0.0;
        bands[9].sections[1].b[0] = 1.0000000000e+00;
        bands[9].sections[1].b[1] = -2.0000000000e+00;
        bands[9].sections[1].b[2] = 1.0000000000e+00;
        bands[9].sections[1].a[0] = 1.0000000000e+00;
        bands[9].sections[1].a[1] = -1.9964846396e+00;
        bands[9].sections[1].a[2] = 9.9685630460e-01;
        bands[9].sections[1].z[0] = 0.0;
        bands[9].sections[1].z[1] = 0.0;

        // Band 10: 200.0 Hz
        bands[10].center_freq = 200.0;
        bands[10].num_sections = 2;
        bands[10].sections[0].b[0] = 9.1486666575e-06;
        bands[10].sections[0].b[1] = 1.8297333315e-05;
        bands[10].sections[0].b[2] = 9.1486666575e-06;
        bands[10].sections[0].a[0] = 1.0000000000e+00;
        bands[10].sections[0].a[1] = -1.9945674699e+00;
        bands[10].sections[0].a[2] = 9.9537316779e-01;
        bands[10].sections[0].z[0] = 0.0;
        bands[10].sections[0].z[1] = 0.0;
        bands[10].sections[1].b[0] = 1.0000000000e+00;
        bands[10].sections[1].b[1] = -2.0000000000e+00;
        bands[10].sections[1].b[2] = 1.0000000000e+00;
        bands[10].sections[1].a[0] = 1.0000000000e+00;
        bands[10].sections[1].a[1] = -1.9954914094e+00;
        bands[10].sections[1].a[2] = 9.9607189798e-01;
        bands[10].sections[1].z[0] = 0.0;
        bands[10].sections[1].z[1] = 0.0;

        // Band 11: 250.0 Hz
        bands[11].center_freq = 250.0;
        bands[11].num_sections = 2;
        bands[11].sections[0].b[0] = 1.4279529675e-05;
        bands[11].sections[0].b[1] = 2.8559059350e-05;
        bands[11].sections[0].b[2] = 1.4279529675e-05;
        bands[11].sections[0].a[0] = 1.0000000000e+00;
        bands[11].sections[0].a[1] = -1.9929617346e+00;
        bands[11].sections[0].a[2] = 9.9421986348e-01;
        bands[11].sections[0].z[0] = 0.0;
        bands[11].sections[0].z[1] = 0.0;
        bands[11].sections[1].b[0] = 1.0000000000e+00;
        bands[11].sections[1].b[1] = -2.0000000000e+00;
        bands[11].sections[1].b[2] = 1.0000000000e+00;
        bands[11].sections[1].a[0] = 1.0000000000e+00;
        bands[11].sections[1].a[1] = -1.9941856867e+00;
        bands[11].sections[1].a[2] = 9.9509223080e-01;
        bands[11].sections[1].z[0] = 0.0;
        bands[11].sections[1].z[1] = 0.0;

        // Band 12: 315.0 Hz
        bands[12].center_freq = 315.0;
        bands[12].num_sections = 2;
        bands[12].sections[0].b[0] = 2.2638747004e-05;
        bands[12].sections[0].b[1] = 4.5277494009e-05;
        bands[12].sections[0].b[2] = 2.2638747004e-05;
        bands[12].sections[0].a[0] = 1.0000000000e+00;
        bands[12].sections[0].a[1] = -1.9907268354e+00;
        bands[12].sections[0].a[2] = 9.9272262119e-01;
        bands[12].sections[0].z[0] = 0.0;
        bands[12].sections[0].z[1] = 0.0;
        bands[12].sections[1].b[0] = 1.0000000000e+00;
        bands[12].sections[1].b[1] = -2.0000000000e+00;
        bands[12].sections[1].b[2] = 1.0000000000e+00;
        bands[12].sections[1].a[0] = 1.0000000000e+00;
        bands[12].sections[1].a[1] = -1.9923817986e+00;
        bands[12].sections[1].a[2] = 9.9382004797e-01;
        bands[12].sections[1].z[0] = 0.0;
        bands[12].sections[1].z[1] = 0.0;

        // Band 13: 400.0 Hz
        bands[13].center_freq = 400.0;
        bands[13].num_sections = 2;
        bands[13].sections[0].b[0] = 3.6438801475e-05;
        bands[13].sections[0].b[1] = 7.2877602950e-05;
        bands[13].sections[0].b[2] = 3.6438801475e-05;
        bands[13].sections[0].a[0] = 1.0000000000e+00;
        bands[13].sections[0].a[1] = -1.9875534886e+00;
        bands[13].sections[0].a[2] = 9.9076821039e-01;
        bands[13].sections[0].z[0] = 0.0;
        bands[13].sections[0].z[1] = 0.0;
        bands[13].sections[1].b[0] = 1.0000000000e+00;
        bands[13].sections[1].b[1] = -2.0000000000e+00;
        bands[13].sections[1].b[2] = 1.0000000000e+00;
        bands[13].sections[1].a[0] = 1.0000000000e+00;
        bands[13].sections[1].a[1] = -1.9898416841e+00;
        bands[13].sections[1].a[2] = 9.9215875807e-01;
        bands[13].sections[1].z[0] = 0.0;
        bands[13].sections[1].z[1] = 0.0;

        // Band 14: 500.0 Hz
        bands[14].center_freq = 500.0;
        bands[14].num_sections = 2;
        bands[14].sections[0].b[0] = 5.6814507890e-05;
        bands[14].sections[0].b[1] = 1.1362901578e-04;
        bands[14].sections[0].b[2] = 5.6814507890e-05;
        bands[14].sections[0].a[0] = 1.0000000000e+00;
        bands[14].sections[0].a[1] = -1.9834575730e+00;
        bands[14].sections[0].a[2] = 9.8847404755e-01;
        bands[14].sections[0].z[0] = 0.0;
        bands[14].sections[0].z[1] = 0.0;
        bands[14].sections[1].b[0] = 1.0000000000e+00;
        bands[14].sections[1].b[1] = -2.0000000000e+00;
        bands[14].sections[1].b[2] = 1.0000000000e+00;
        bands[14].sections[1].a[0] = 1.0000000000e+00;
        bands[14].sections[1].a[1] = -1.9865911375e+00;
        bands[14].sections[1].a[2] = 9.9020763562e-01;
        bands[14].sections[1].z[0] = 0.0;
        bands[14].sections[1].z[1] = 0.0;

        // Band 15: 630.0 Hz
        bands[15].center_freq = 630.0;
        bands[15].num_sections = 2;
        bands[15].sections[0].b[0] = 8.9949759973e-05;
        bands[15].sections[0].b[1] = 1.7989951995e-04;
        bands[15].sections[0].b[2] = 8.9949759973e-05;
        bands[15].sections[0].a[0] = 1.0000000000e+00;
        bands[15].sections[0].a[1] = -1.9775496894e+00;
        bands[15].sections[0].a[2] = 9.8550001917e-01;
        bands[15].sections[0].z[0] = 0.0;
        bands[15].sections[0].z[1] = 0.0;
        bands[15].sections[1].b[0] = 1.0000000000e+00;
        bands[15].sections[1].b[1] = -2.0000000000e+00;
        bands[15].sections[1].b[2] = 1.0000000000e+00;
        bands[15].sections[1].a[0] = 1.0000000000e+00;
        bands[15].sections[1].a[1] = -1.9819432110e+00;
        bands[15].sections[1].a[2] = 9.8767646726e-01;
        bands[15].sections[1].z[0] = 0.0;
        bands[15].sections[1].z[1] = 0.0;

        // Band 16: 800.0 Hz
        bands[16].center_freq = 800.0;
        bands[16].num_sections = 2;
        bands[16].sections[0].b[0] = 1.4452155018e-04;
        bands[16].sections[0].b[1] = 2.8904310037e-04;
        bands[16].sections[0].b[2] = 1.4452155018e-04;
        bands[16].sections[0].a[0] = 1.0000000000e+00;
        bands[16].sections[0].a[1] = -1.9688355849e+00;
        bands[16].sections[0].a[2] = 9.8162535328e-01;
        bands[16].sections[0].z[0] = 0.0;
        bands[16].sections[0].z[1] = 0.0;
        bands[16].sections[1].b[0] = 1.0000000000e+00;
        bands[16].sections[1].b[1] = -2.0000000000e+00;
        bands[16].sections[1].b[2] = 1.0000000000e+00;
        bands[16].sections[1].a[0] = 1.0000000000e+00;
        bands[16].sections[1].a[1] = -1.9751484285e+00;
        bands[16].sections[1].a[2] = 9.8437528452e-01;
        bands[16].sections[1].z[0] = 0.0;
        bands[16].sections[1].z[1] = 0.0;

        // Band 17: 1000.0 Hz
        bands[17].center_freq = 1000.0;
        bands[17].num_sections = 2;
        bands[17].sections[0].b[0] = 2.2486138577e-04;
        bands[17].sections[0].b[1] = 4.4972277155e-04;
        bands[17].sections[0].b[2] = 2.2486138577e-04;
        bands[17].sections[0].a[0] = 1.0000000000e+00;
        bands[17].sections[0].a[1] = -1.9571616553e+00;
        bands[17].sections[0].a[2] = 9.7708815318e-01;
        bands[17].sections[0].z[0] = 0.0;
        bands[17].sections[0].z[1] = 0.0;
        bands[17].sections[1].b[0] = 1.0000000000e+00;
        bands[17].sections[1].b[1] = -2.0000000000e+00;
        bands[17].sections[1].b[2] = 1.0000000000e+00;
        bands[17].sections[1].a[0] = 1.0000000000e+00;
        bands[17].sections[1].a[1] = -1.9661212031e+00;
        bands[17].sections[1].a[2] = 9.8050392671e-01;
        bands[17].sections[1].z[0] = 0.0;
        bands[17].sections[1].z[1] = 0.0;

        // Band 18: 1250.0 Hz
        bands[18].center_freq = 1250.0;
        bands[18].num_sections = 2;
        bands[18].sections[0].b[0] = 3.4949840016e-04;
        bands[18].sections[0].b[1] = 6.9899680031e-04;
        bands[18].sections[0].b[2] = 3.4949840016e-04;
        bands[18].sections[0].a[0] = 1.0000000000e+00;
        bands[18].sections[0].a[1] = -1.9404316792e+00;
        bands[18].sections[0].a[2] = 9.7144942329e-01;
        bands[18].sections[0].z[0] = 0.0;
        bands[18].sections[0].z[1] = 0.0;
        bands[18].sections[1].b[0] = 1.0000000000e+00;
        bands[18].sections[1].b[1] = -2.0000000000e+00;
        bands[18].sections[1].b[2] = 1.0000000000e+00;
        bands[18].sections[1].a[0] = 1.0000000000e+00;
        bands[18].sections[1].a[1] = -1.9532794096e+00;
        bands[18].sections[1].a[2] = 9.7568282934e-01;
        bands[18].sections[1].z[0] = 0.0;
        bands[18].sections[1].z[1] = 0.0;

        // Band 19: 1600.0 Hz
        bands[19].center_freq = 1600.0;
        bands[19].num_sections = 2;
        bands[19].sections[0].b[0] = 5.6842540279e-04;
        bands[19].sections[0].b[1] = 1.1368508056e-03;
        bands[19].sections[0].b[2] = 5.6842540279e-04;
        bands[19].sections[0].a[0] = 1.0000000000e+00;
        bands[19].sections[0].a[1] = -1.9130826710e+00;
        bands[19].sections[0].a[2] = 9.6361754164e-01;
        bands[19].sections[0].z[0] = 0.0;
        bands[19].sections[0].z[1] = 0.0;
        bands[19].sections[1].b[0] = 1.0000000000e+00;
        bands[19].sections[1].b[1] = -2.0000000000e+00;
        bands[19].sections[1].b[2] = 1.0000000000e+00;
        bands[19].sections[1].a[0] = 1.0000000000e+00;
        bands[19].sections[1].a[1] = -1.9324274349e+00;
        bands[19].sections[1].a[2] = 9.6896533897e-01;
        bands[19].sections[1].z[0] = 0.0;
        bands[19].sections[1].z[1] = 0.0;

        // Band 20: 2000.0 Hz
        bands[20].center_freq = 2000.0;
        bands[20].num_sections = 2;
        bands[20].sections[0].b[0] = 8.8077669909e-04;
        bands[20].sections[0].b[1] = 1.7615533982e-03;
        bands[20].sections[0].b[2] = 8.8077669909e-04;
        bands[20].sections[0].a[0] = 1.0000000000e+00;
        bands[20].sections[0].a[1] = -1.8763388948e+00;
        bands[20].sections[0].a[2] = 9.5475792455e-01;
        bands[20].sections[0].z[0] = 0.0;
        bands[20].sections[0].z[1] = 0.0;
        bands[20].sections[1].b[0] = 1.0000000000e+00;
        bands[20].sections[1].b[1] = -2.0000000000e+00;
        bands[20].sections[1].b[2] = 1.0000000000e+00;
        bands[20].sections[1].a[0] = 1.0000000000e+00;
        bands[20].sections[1].a[1] = -1.9045584369e+00;
        bands[20].sections[1].a[2] = 9.6133091745e-01;
        bands[20].sections[1].z[0] = 0.0;
        bands[20].sections[1].z[1] = 0.0;

        // Band 21: 2500.0 Hz
        bands[21].center_freq = 2500.0;
        bands[21].num_sections = 2;
        bands[21].sections[0].b[0] = 1.3620124551e-03;
        bands[21].sections[0].b[1] = 2.7240249103e-03;
        bands[21].sections[0].b[2] = 1.3620124551e-03;
        bands[21].sections[0].a[0] = 1.0000000000e+00;
        bands[21].sections[0].a[1] = -1.8224286352e+00;
        bands[21].sections[0].a[2] = 9.4382382636e-01;
        bands[21].sections[0].z[0] = 0.0;
        bands[21].sections[0].z[1] = 0.0;
        bands[21].sections[1].b[0] = 1.0000000000e+00;
        bands[21].sections[1].b[1] = -2.0000000000e+00;
        bands[21].sections[1].b[2] = 1.0000000000e+00;
        bands[21].sections[1].a[0] = 1.0000000000e+00;
        bands[21].sections[1].a[1] = -1.8638026219e+00;
        bands[21].sections[1].a[2] = 9.5184627152e-01;
        bands[21].sections[1].z[0] = 0.0;
        bands[21].sections[1].z[1] = 0.0;

        // Band 22: 3150.0 Hz
        bands[22].center_freq = 3150.0;
        bands[22].num_sections = 2;
        bands[22].sections[0].b[0] = 2.1336121687e-03;
        bands[22].sections[0].b[1] = 4.2672243374e-03;
        bands[22].sections[0].b[2] = 2.1336121687e-03;
        bands[22].sections[0].a[0] = 1.0000000000e+00;
        bands[22].sections[0].a[1] = -1.7396600701e+00;
        bands[22].sections[0].a[2] = 9.2985039510e-01;
        bands[22].sections[0].z[0] = 0.0;
        bands[22].sections[0].z[1] = 0.0;
        bands[22].sections[1].b[0] = 1.0000000000e+00;
        bands[22].sections[1].b[1] = -2.0000000000e+00;
        bands[22].sections[1].b[2] = 1.0000000000e+00;
        bands[22].sections[1].a[0] = 1.0000000000e+00;
        bands[22].sections[1].a[1] = -1.8013022953e+00;
        bands[22].sections[1].a[2] = 9.3960162339e-01;
        bands[22].sections[1].z[0] = 0.0;
        bands[22].sections[1].z[1] = 0.0;

        // Band 23: 4000.0 Hz
        bands[23].center_freq = 4000.0;
        bands[23].num_sections = 2;
        bands[23].sections[0].b[0] = 3.3814676189e-03;
        bands[23].sections[0].b[1] = 6.7629352377e-03;
        bands[23].sections[0].b[2] = 3.3814676189e-03;
        bands[23].sections[0].a[0] = 1.0000000000e+00;
        bands[23].sections[0].a[1] = -1.6111845677e+00;
        bands[23].sections[0].a[2] = 9.1200507593e-01;
        bands[23].sections[0].z[0] = 0.0;
        bands[23].sections[0].z[1] = 0.0;
        bands[23].sections[1].b[0] = 1.0000000000e+00;
        bands[23].sections[1].b[1] = -2.0000000000e+00;
        bands[23].sections[1].b[2] = 1.0000000000e+00;
        bands[23].sections[1].a[0] = 1.0000000000e+00;
        bands[23].sections[1].a[1] = -1.7041154148e+00;
        bands[23].sections[1].a[2] = 9.2370966205e-01;
        bands[23].sections[1].z[0] = 0.0;
        bands[23].sections[1].z[1] = 0.0;

        // Band 24: 5000.0 Hz
        bands[24].center_freq = 5000.0;
        bands[24].num_sections = 2;
        bands[24].sections[0].b[0] = 5.1786052377e-03;
        bands[24].sections[0].b[1] = 1.0357210475e-02;
        bands[24].sections[0].b[2] = 5.1786052377e-03;
        bands[24].sections[0].a[0] = 1.0000000000e+00;
        bands[24].sections[0].a[1] = -1.4334898271e+00;
        bands[24].sections[0].a[2] = 8.9166200446e-01;
        bands[24].sections[0].z[0] = 0.0;
        bands[24].sections[0].z[1] = 0.0;
        bands[24].sections[1].b[0] = 1.0000000000e+00;
        bands[24].sections[1].b[1] = -2.0000000000e+00;
        bands[24].sections[1].b[2] = 1.0000000000e+00;
        bands[24].sections[1].a[0] = 1.0000000000e+00;
        bands[24].sections[1].a[1] = -1.5689457900e+00;
        bands[24].sections[1].a[2] = 9.0514133105e-01;
        bands[24].sections[1].z[0] = 0.0;
        bands[24].sections[1].z[1] = 0.0;

        // Band 25: 6300.0 Hz
        bands[25].center_freq = 6300.0;
        bands[25].num_sections = 2;
        bands[25].sections[0].b[0] = 8.0135241686e-03;
        bands[25].sections[0].b[1] = 1.6027048337e-02;
        bands[25].sections[0].b[2] = 8.0135241686e-03;
        bands[25].sections[0].a[0] = 1.0000000000e+00;
        bands[25].sections[0].a[1] = -1.1659539968e+00;
        bands[25].sections[0].a[2] = 8.6633608642e-01;
        bands[25].sections[0].z[0] = 0.0;
        bands[25].sections[0].z[1] = 0.0;
        bands[25].sections[1].b[0] = 1.0000000000e+00;
        bands[25].sections[1].b[1] = -2.0000000000e+00;
        bands[25].sections[1].b[2] = 1.0000000000e+00;
        bands[25].sections[1].a[0] = 1.0000000000e+00;
        bands[25].sections[1].a[1] = -1.3632319368e+00;
        bands[25].sections[1].a[2] = 8.8111168516e-01;
        bands[25].sections[1].z[0] = 0.0;
        bands[25].sections[1].z[1] = 0.0;

        // Band 26: 8000.0 Hz
        bands[26].center_freq = 8000.0;
        bands[26].num_sections = 2;
        bands[26].sections[0].b[0] = 1.2505435661e-02;
        bands[26].sections[0].b[1] = 2.5010871322e-02;
        bands[26].sections[0].b[2] = 1.2505435661e-02;
        bands[26].sections[0].a[0] = 1.0000000000e+00;
        bands[26].sections[0].a[1] = -7.6892892667e-01;
        bands[26].sections[0].a[2] = 8.3529708847e-01;
        bands[26].sections[0].z[0] = 0.0;
        bands[26].sections[0].z[1] = 0.0;
        bands[26].sections[1].b[0] = 1.0000000000e+00;
        bands[26].sections[1].b[1] = -2.0000000000e+00;
        bands[26].sections[1].b[2] = 1.0000000000e+00;
        bands[26].sections[1].a[0] = 1.0000000000e+00;
        bands[26].sections[1].a[1] = -1.0520645541e+00;
        bands[26].sections[1].a[2] = 8.4964354776e-01;
        bands[26].sections[1].z[0] = 0.0;
        bands[26].sections[1].z[1] = 0.0;

        // Band 27: 10000.0 Hz
        bands[27].center_freq = 10000.0;
        bands[27].num_sections = 2;
        bands[27].sections[0].b[0] = 1.8821791870e-02;
        bands[27].sections[0].b[1] = 3.7643583740e-02;
        bands[27].sections[0].b[2] = 1.8821791870e-02;
        bands[27].sections[0].a[0] = 1.0000000000e+00;
        bands[27].sections[0].a[1] = -2.6268885378e-01;
        bands[27].sections[0].a[2] = 8.0220300414e-01;
        bands[27].sections[0].z[0] = 0.0;
        bands[27].sections[0].z[1] = 0.0;
        bands[27].sections[1].b[0] = 1.0000000000e+00;
        bands[27].sections[1].b[1] = -2.0000000000e+00;
        bands[27].sections[1].b[2] = 1.0000000000e+00;
        bands[27].sections[1].a[0] = 1.0000000000e+00;
        bands[27].sections[1].a[1] = -6.4302702672e-01;
        bands[27].sections[1].a[2] = 8.1205737082e-01;
        bands[27].sections[1].z[0] = 0.0;
        bands[27].sections[1].z[1] = 0.0;

        // Band 28: 12500.0 Hz
        bands[28].center_freq = 12500.0;
        bands[28].num_sections = 2;
        bands[28].sections[0].b[0] = 2.8110768870e-02;
        bands[28].sections[0].b[1] = -5.6221537740e-02;
        bands[28].sections[0].b[2] = 2.8110768870e-02;
        bands[28].sections[0].a[0] = 1.0000000000e+00;
        bands[28].sections[0].a[1] = -1.0383271516e-01;
        bands[28].sections[0].a[2] = 7.6289180572e-01;
        bands[28].sections[0].z[0] = 0.0;
        bands[28].sections[0].z[1] = 0.0;
        bands[28].sections[1].b[0] = 1.0000000000e+00;
        bands[28].sections[1].b[1] = 2.0000000000e+00;
        bands[28].sections[1].b[2] = 1.0000000000e+00;
        bands[28].sections[1].a[0] = 1.0000000000e+00;
        bands[28].sections[1].a[1] = 3.7366504960e-01;
        bands[28].sections[1].a[2] = 7.6727247890e-01;
        bands[28].sections[1].z[0] = 0.0;
        bands[28].sections[1].z[1] = 0.0;

        // Band 29: 16000.0 Hz
        bands[29].center_freq = 16000.0;
        bands[29].num_sections = 2;
        bands[29].sections[0].b[0] = 4.3369979301e-02;
        bands[29].sections[0].b[1] = -8.6739958603e-02;
        bands[29].sections[0].b[2] = 4.3369979301e-02;
        bands[29].sections[0].a[0] = 1.0000000000e+00;
        bands[29].sections[0].a[1] = 6.0774354142e-01;
        bands[29].sections[0].a[2] = 6.8364720419e-01;
        bands[29].sections[0].z[0] = 0.0;
        bands[29].sections[0].z[1] = 0.0;
        bands[29].sections[1].b[0] = 1.0000000000e+00;
        bands[29].sections[1].b[1] = 2.0000000000e+00;
        bands[29].sections[1].b[2] = 1.0000000000e+00;
        bands[29].sections[1].a[0] = 1.0000000000e+00;
        bands[29].sections[1].a[1] = 1.1492032881e+00;
        bands[29].sections[1].a[2] = 7.3744978539e-01;
        bands[29].sections[1].z[0] = 0.0;
        bands[29].sections[1].z[1] = 0.0;

        // Band 30: 20000.0 Hz
        bands[30].center_freq = 20000.0;
        bands[30].num_sections = 2;
        bands[30].sections[0].b[0] = 6.3852757781e-02;
        bands[30].sections[0].b[1] = -1.2770551556e-01;
        bands[30].sections[0].b[2] = 6.3852757781e-02;
        bands[30].sections[0].a[0] = 1.0000000000e+00;
        bands[30].sections[0].a[1] = 1.1475077614e+00;
        bands[30].sections[0].a[2] = 5.2760390543e-01;
        bands[30].sections[0].z[0] = 0.0;
        bands[30].sections[0].z[1] = 0.0;
        bands[30].sections[1].b[0] = 1.0000000000e+00;
        bands[30].sections[1].b[1] = 2.0000000000e+00;
        bands[30].sections[1].b[2] = 1.0000000000e+00;
        bands[30].sections[1].a[0] = 1.0000000000e+00;
        bands[30].sections[1].a[1] = 1.7589790538e+00;
        bands[30].sections[1].a[2] = 8.0671714926e-01;
        bands[30].sections[1].z[0] = 0.0;
        bands[30].sections[1].z[1] = 0.0;

    } else { // filter_order == 4
        // Band 0: 20.0 Hz
        bands[0].center_freq = 20.0;
        bands[0].num_sections = 4;
        bands[0].sections[0].b[0] = 8.4350901384e-15;
        bands[0].sections[0].b[1] = 1.6870180277e-14;
        bands[0].sections[0].b[2] = 8.4350901384e-15;
        bands[0].sections[0].a[0] = 1.0000000000e+00;
        bands[0].sections[0].a[1] = -1.9994076617e+00;
        bands[0].sections[0].a[2] = 9.9941515210e-01;
        bands[0].sections[0].z[0] = 0.0;
        bands[0].sections[0].z[1] = 0.0;
        bands[0].sections[1].b[0] = 1.0000000000e+00;
        bands[0].sections[1].b[1] = 2.0000000000e+00;
        bands[0].sections[1].b[2] = 1.0000000000e+00;
        bands[0].sections[1].a[0] = 1.0000000000e+00;
        bands[0].sections[1].a[1] = -1.9994587246e+00;
        bands[0].sections[1].a[2] = 9.9946499258e-01;
        bands[0].sections[1].z[0] = 0.0;
        bands[0].sections[1].z[1] = 0.0;
        bands[0].sections[2].b[0] = 1.0000000000e+00;
        bands[0].sections[2].b[1] = -2.0000000000e+00;
        bands[0].sections[2].b[2] = 1.0000000000e+00;
        bands[0].sections[2].a[0] = 1.0000000000e+00;
        bands[0].sections[2].a[1] = -1.9997348532e+00;
        bands[0].sections[2].a[2] = 9.9974333928e-01;
        bands[0].sections[2].z[0] = 0.0;
        bands[0].sections[2].z[1] = 0.0;
        bands[0].sections[3].b[0] = 1.0000000000e+00;
        bands[0].sections[3].b[1] = -2.0000000000e+00;
        bands[0].sections[3].b[2] = 1.0000000000e+00;
        bands[0].sections[3].a[0] = 1.0000000000e+00;
        bands[0].sections[3].a[1] = -1.9997871914e+00;
        bands[0].sections[3].a[2] = 9.9979272575e-01;
        bands[0].sections[3].z[0] = 0.0;
        bands[0].sections[3].z[1] = 0.0;

        // Band 1: 25.0 Hz
        bands[1].center_freq = 25.0;
        bands[1].num_sections = 4;
        bands[1].sections[0].b[0] = 2.0589405689e-14;
        bands[1].sections[0].b[1] = 4.1178811379e-14;
        bands[1].sections[0].b[2] = 2.0589405689e-14;
        bands[1].sections[0].a[0] = 1.0000000000e+00;
        bands[1].sections[0].a[1] = -1.9992572907e+00;
        bands[1].sections[0].a[2] = 9.9926899360e-01;
        bands[1].sections[0].z[0] = 0.0;
        bands[1].sections[0].z[1] = 0.0;
        bands[1].sections[1].b[0] = 1.0000000000e+00;
        bands[1].sections[1].b[1] = 2.0000000000e+00;
        bands[1].sections[1].b[2] = 1.0000000000e+00;
        bands[1].sections[1].a[0] = 1.0000000000e+00;
        bands[1].sections[1].a[1] = -1.9993214924e+00;
        bands[1].sections[1].a[2] = 9.9933128540e-01;
        bands[1].sections[1].z[0] = 0.0;
        bands[1].sections[1].z[1] = 0.0;
        bands[1].sections[2].b[0] = 1.0000000000e+00;
        bands[1].sections[2].b[1] = -2.0000000000e+00;
        bands[1].sections[2].b[2] = 1.0000000000e+00;
        bands[1].sections[2].a[0] = 1.0000000000e+00;
        bands[1].sections[2].a[1] = -1.9996659254e+00;
        bands[1].sections[2].a[2] = 9.9967918444e-01;
        bands[1].sections[2].z[0] = 0.0;
        bands[1].sections[2].z[1] = 0.0;
        bands[1].sections[3].b[0] = 1.0000000000e+00;
        bands[1].sections[3].b[1] = -2.0000000000e+00;
        bands[1].sections[3].b[2] = 1.0000000000e+00;
        bands[1].sections[3].a[0] = 1.0000000000e+00;
        bands[1].sections[3].a[1] = -1.9997322667e+00;
        bands[1].sections[3].a[2] = 9.9974091387e-01;
        bands[1].sections[3].z[0] = 0.0;
        bands[1].sections[3].z[1] = 0.0;

        // Band 2: 31.5 Hz
        bands[2].center_freq = 31.5;
        bands[2].num_sections = 4;
        bands[2].sections[0].b[0] = 5.1881705198e-14;
        bands[2].sections[0].b[1] = 1.0376341040e-13;
        bands[2].sections[0].b[2] = 5.1881705198e-14;
        bands[2].sections[0].a[0] = 1.0000000000e+00;
        bands[2].sections[0].a[1] = -1.9990604418e+00;
        bands[2].sections[0].a[2] = 9.9907901952e-01;
        bands[2].sections[0].z[0] = 0.0;
        bands[2].sections[0].z[1] = 0.0;
        bands[2].sections[1].b[0] = 1.0000000000e+00;
        bands[2].sections[1].b[1] = 2.0000000000e+00;
        bands[2].sections[1].b[2] = 1.0000000000e+00;
        bands[2].sections[1].a[0] = 1.0000000000e+00;
        bands[2].sections[1].a[1] = -1.9991419467e+00;
        bands[2].sections[1].a[2] = 9.9915749275e-01;
        bands[2].sections[1].z[0] = 0.0;
        bands[2].sections[1].z[1] = 0.0;
        bands[2].sections[2].b[0] = 1.0000000000e+00;
        bands[2].sections[2].b[1] = -2.0000000000e+00;
        bands[2].sections[2].b[2] = 1.0000000000e+00;
        bands[2].sections[2].a[0] = 1.0000000000e+00;
        bands[2].sections[2].a[1] = -1.9995747402e+00;
        bands[2].sections[2].a[2] = 9.9959578936e-01;
        bands[2].sections[2].z[0] = 0.0;
        bands[2].sections[2].z[1] = 0.0;
        bands[2].sections[3].b[0] = 1.0000000000e+00;
        bands[2].sections[3].b[1] = -2.0000000000e+00;
        bands[2].sections[3].b[2] = 1.0000000000e+00;
        bands[2].sections[3].a[0] = 1.0000000000e+00;
        bands[2].sections[3].a[1] = -1.9996598346e+00;
        bands[2].sections[3].a[2] = 9.9967356241e-01;
        bands[2].sections[3].z[0] = 0.0;
        bands[2].sections[3].z[1] = 0.0;

        // Band 3: 40.0 Hz
        bands[3].center_freq = 40.0;
        bands[3].num_sections = 4;
        bands[3].sections[0].b[0] = 1.3485463392e-13;
        bands[3].sections[0].b[1] = 2.6970926783e-13;
        bands[3].sections[0].b[2] = 1.3485463392e-13;
        bands[3].sections[0].a[0] = 1.0000000000e+00;
        bands[3].sections[0].a[1] = -1.9988006937e+00;
        bands[3].sections[0].a[2] = 9.9883064642e-01;
        bands[3].sections[0].z[0] = 0.0;
        bands[3].sections[0].z[1] = 0.0;
        bands[3].sections[1].b[0] = 1.0000000000e+00;
        bands[3].sections[1].b[1] = 2.0000000000e+00;
        bands[3].sections[1].b[2] = 1.0000000000e+00;
        bands[3].sections[1].a[0] = 1.0000000000e+00;
        bands[3].sections[1].a[1] = -1.9989052058e+00;
        bands[3].sections[1].a[2] = 9.9893027092e-01;
        bands[3].sections[1].z[0] = 0.0;
        bands[3].sections[1].z[1] = 0.0;
        bands[3].sections[2].b[0] = 1.0000000000e+00;
        bands[3].sections[2].b[1] = -2.0000000000e+00;
        bands[3].sections[2].b[2] = 1.0000000000e+00;
        bands[3].sections[2].a[0] = 1.0000000000e+00;
        bands[3].sections[2].a[1] = -1.9994528050e+00;
        bands[3].sections[2].a[2] = 9.9948674484e-01;
        bands[3].sections[2].z[0] = 0.0;
        bands[3].sections[2].z[1] = 0.0;
        bands[3].sections[3].b[0] = 1.0000000000e+00;
        bands[3].sections[3].b[1] = -2.0000000000e+00;
        bands[3].sections[3].b[2] = 1.0000000000e+00;
        bands[3].sections[3].a[0] = 1.0000000000e+00;
        bands[3].sections[3].a[1] = -1.9995633591e+00;
        bands[3].sections[3].a[2] = 9.9958549419e-01;
        bands[3].sections[3].z[0] = 0.0;
        bands[3].sections[3].z[1] = 0.0;

        // Band 4: 50.0 Hz
        bands[4].center_freq = 50.0;
        bands[4].num_sections = 4;
        bands[4].sections[0].b[0] = 3.2910467246e-13;
        bands[4].sections[0].b[1] = 6.5820934491e-13;
        bands[4].sections[0].b[2] = 3.2910467246e-13;
        bands[4].sections[0].a[0] = 1.0000000000e+00;
        bands[4].sections[0].a[1] = -1.9984917277e+00;
        bands[4].sections[0].a[2] = 9.9853852191e-01;
        bands[4].sections[0].z[0] = 0.0;
        bands[4].sections[0].z[1] = 0.0;
        bands[4].sections[1].b[0] = 1.0000000000e+00;
        bands[4].sections[1].b[1] = 2.0000000000e+00;
        bands[4].sections[1].b[2] = 1.0000000000e+00;
        bands[4].sections[1].a[0] = 1.0000000000e+00;
        bands[4].sections[1].a[1] = -1.9986238581e+00;
        bands[4].sections[1].a[2] = 9.9866301706e-01;
        bands[4].sections[1].z[0] = 0.0;
        bands[4].sections[1].z[1] = 0.0;
        bands[4].sections[2].b[0] = 1.0000000000e+00;
        bands[4].sections[2].b[1] = -2.0000000000e+00;
        bands[4].sections[2].b[2] = 1.0000000000e+00;
        bands[4].sections[2].a[0] = 1.0000000000e+00;
        bands[4].sections[2].a[1] = -1.9993054451e+00;
        bands[4].sections[2].a[2] = 9.9935847260e-01;
        bands[4].sections[2].z[0] = 0.0;
        bands[4].sections[2].z[1] = 0.0;
        bands[4].sections[3].b[0] = 1.0000000000e+00;
        bands[4].sections[3].b[1] = -2.0000000000e+00;
        bands[4].sections[3].b[2] = 1.0000000000e+00;
        bands[4].sections[3].a[0] = 1.0000000000e+00;
        bands[4].sections[3].a[1] = -1.9994473101e+00;
        bands[4].sections[3].a[2] = 9.9948189432e-01;
        bands[4].sections[3].z[0] = 0.0;
        bands[4].sections[3].z[1] = 0.0;

        // Band 5: 63.0 Hz
        bands[5].center_freq = 63.0;
        bands[5].num_sections = 4;
        bands[5].sections[0].b[0] = 8.2907310466e-13;
        bands[5].sections[0].b[1] = 1.6581462093e-12;
        bands[5].sections[0].b[2] = 8.2907310466e-13;
        bands[5].sections[0].a[0] = 1.0000000000e+00;
        bands[5].sections[0].a[1] = -1.9980846117e+00;
        bands[5].sections[0].a[2] = 9.9815888794e-01;
        bands[5].sections[0].z[0] = 0.0;
        bands[5].sections[0].z[1] = 0.0;
        bands[5].sections[1].b[0] = 1.0000000000e+00;
        bands[5].sections[1].b[1] = 2.0000000000e+00;
        bands[5].sections[1].b[2] = 1.0000000000e+00;
        bands[5].sections[1].a[0] = 1.0000000000e+00;
        bands[5].sections[1].a[1] = -1.9982535356e+00;
        bands[5].sections[1].a[2] = 9.9831569349e-01;
        bands[5].sections[1].z[0] = 0.0;
        bands[5].sections[1].z[1] = 0.0;
        bands[5].sections[2].b[0] = 1.0000000000e+00;
        bands[5].sections[2].b[1] = -2.0000000000e+00;
        bands[5].sections[2].b[2] = 1.0000000000e+00;
        bands[5].sections[2].a[0] = 1.0000000000e+00;
        bands[5].sections[2].a[1] = -1.9991075645e+00;
        bands[5].sections[2].a[2] = 9.9919174367e-01;
        bands[5].sections[2].z[0] = 0.0;
        bands[5].sections[2].z[1] = 0.0;
        bands[5].sections[3].b[0] = 1.0000000000e+00;
        bands[5].sections[3].b[1] = -2.0000000000e+00;
        bands[5].sections[3].b[2] = 1.0000000000e+00;
        bands[5].sections[3].a[0] = 1.0000000000e+00;
        bands[5].sections[3].a[1] = -1.9992923281e+00;
        bands[5].sections[3].a[2] = 9.9934723027e-01;
        bands[5].sections[3].z[0] = 0.0;
        bands[5].sections[3].z[1] = 0.0;

        // Band 6: 80.0 Hz
        bands[6].center_freq = 80.0;
        bands[6].num_sections = 4;
        bands[6].sections[0].b[0] = 2.1542619151e-12;
        bands[6].sections[0].b[1] = 4.3085238301e-12;
        bands[6].sections[0].b[2] = 2.1542619151e-12;
        bands[6].sections[0].a[0] = 1.0000000000e+00;
        bands[6].sections[0].a[1] = -1.9975429216e+00;
        bands[6].sections[0].a[2] = 9.9766266166e-01;
        bands[6].sections[0].z[0] = 0.0;
        bands[6].sections[0].z[1] = 0.0;
        bands[6].sections[1].b[0] = 1.0000000000e+00;
        bands[6].sections[1].b[1] = 2.0000000000e+00;
        bands[6].sections[1].b[2] = 1.0000000000e+00;
        bands[6].sections[1].a[0] = 1.0000000000e+00;
        bands[6].sections[1].a[1] = -1.9977614761e+00;
        bands[6].sections[1].a[2] = 9.9786168242e-01;
        bands[6].sections[1].z[0] = 0.0;
        bands[6].sections[1].z[1] = 0.0;
        bands[6].sections[2].b[0] = 1.0000000000e+00;
        bands[6].sections[2].b[1] = -2.0000000000e+00;
        bands[6].sections[2].b[2] = 1.0000000000e+00;
        bands[6].sections[2].a[0] = 1.0000000000e+00;
        bands[6].sections[2].a[1] = -1.9988380330e+00;
        bands[6].sections[2].a[2] = 9.9897375634e-01;
        bands[6].sections[2].z[0] = 0.0;
        bands[6].sections[2].z[1] = 0.0;
        bands[6].sections[3].b[0] = 1.0000000000e+00;
        bands[6].sections[3].b[1] = -2.0000000000e+00;
        bands[6].sections[3].b[2] = 1.0000000000e+00;
        bands[6].sections[3].a[0] = 1.0000000000e+00;
        bands[6].sections[3].a[1] = -1.9990826364e+00;
        bands[6].sections[3].a[2] = 9.9917115793e-01;
        bands[6].sections[3].z[0] = 0.0;
        bands[6].sections[3].z[1] = 0.0;

        // Band 7: 100.0 Hz
        bands[7].center_freq = 100.0;
        bands[7].num_sections = 4;
        bands[7].sections[0].b[0] = 5.2552700572e-12;
        bands[7].sections[0].b[1] = 1.0510540114e-11;
        bands[7].sections[0].b[2] = 5.2552700572e-12;
        bands[7].sections[0].a[0] = 1.0000000000e+00;
        bands[7].sections[0].a[1] = -1.9968921444e+00;
        bands[7].sections[0].a[2] = 9.9707918252e-01;
        bands[7].sections[0].z[0] = 0.0;
        bands[7].sections[0].z[1] = 0.0;
        bands[7].sections[1].b[0] = 1.0000000000e+00;
        bands[7].sections[1].b[1] = 2.0000000000e+00;
        bands[7].sections[1].b[2] = 1.0000000000e+00;
        bands[7].sections[1].a[0] = 1.0000000000e+00;
        bands[7].sections[1].a[1] = -1.9971712845e+00;
        bands[7].sections[1].a[2] = 9.9732781434e-01;
        bands[7].sections[1].z[0] = 0.0;
        bands[7].sections[1].z[1] = 0.0;
        bands[7].sections[2].b[0] = 1.0000000000e+00;
        bands[7].sections[2].b[1] = -2.0000000000e+00;
        bands[7].sections[2].b[2] = 1.0000000000e+00;
        bands[7].sections[2].a[0] = 1.0000000000e+00;
        bands[7].sections[2].a[1] = -1.9985053239e+00;
        bands[7].sections[2].a[2] = 9.9871736304e-01;
        bands[7].sections[2].z[0] = 0.0;
        bands[7].sections[2].z[1] = 0.0;
        bands[7].sections[3].b[0] = 1.0000000000e+00;
        bands[7].sections[3].b[1] = -2.0000000000e+00;
        bands[7].sections[3].b[2] = 1.0000000000e+00;
        bands[7].sections[3].a[0] = 1.0000000000e+00;
        bands[7].sections[3].a[1] = -1.9988257526e+00;
        bands[7].sections[3].a[2] = 9.9896405266e-01;
        bands[7].sections[3].z[0] = 0.0;
        bands[7].sections[3].z[1] = 0.0;

        // Band 8: 125.0 Hz
        bands[8].center_freq = 125.0;
        bands[8].num_sections = 4;
        bands[8].sections[0].b[0] = 1.2817574288e-11;
        bands[8].sections[0].b[1] = 2.5635148576e-11;
        bands[8].sections[0].b[2] = 1.2817574288e-11;
        bands[8].sections[0].a[0] = 1.0000000000e+00;
        bands[8].sections[0].a[1] = -1.9960581768e+00;
        bands[8].sections[0].a[2] = 9.9635031472e-01;
        bands[8].sections[0].z[0] = 0.0;
        bands[8].sections[0].z[1] = 0.0;
        bands[8].sections[1].b[0] = 1.0000000000e+00;
        bands[8].sections[1].b[1] = 2.0000000000e+00;
        bands[8].sections[1].b[2] = 1.0000000000e+00;
        bands[8].sections[1].a[0] = 1.0000000000e+00;
        bands[8].sections[1].a[1] = -1.9964163830e+00;
        bands[8].sections[1].a[2] = 9.9666087754e-01;
        bands[8].sections[1].z[0] = 0.0;
        bands[8].sections[1].z[1] = 0.0;
        bands[8].sections[2].b[0] = 1.0000000000e+00;
        bands[8].sections[2].b[1] = -2.0000000000e+00;
        bands[8].sections[2].b[2] = 1.0000000000e+00;
        bands[8].sections[2].a[0] = 1.0000000000e+00;
        bands[8].sections[2].a[1] = -1.9980657121e+00;
        bands[8].sections[2].a[2] = 9.9839696684e-01;
        bands[8].sections[2].z[0] = 0.0;
        bands[8].sections[2].z[1] = 0.0;
        bands[8].sections[3].b[0] = 1.0000000000e+00;
        bands[8].sections[3].b[1] = -2.0000000000e+00;
        bands[8].sections[3].b[2] = 1.0000000000e+00;
        bands[8].sections[3].a[0] = 1.0000000000e+00;
        bands[8].sections[3].a[1] = -1.9984891650e+00;
        bands[8].sections[3].a[2] = 9.9870522942e-01;
        bands[8].sections[3].z[0] = 0.0;
        bands[8].sections[3].z[1] = 0.0;

        // Band 9: 160.0 Hz
        bands[9].center_freq = 160.0;
        bands[9].num_sections = 4;
        bands[9].sections[0].b[0] = 3.4359359096e-11;
        bands[9].sections[0].b[1] = 6.8718718193e-11;
        bands[9].sections[0].b[2] = 3.4359359096e-11;
        bands[9].sections[0].a[0] = 1.0000000000e+00;
        bands[9].sections[0].a[1] = -1.9948524105e+00;
        bands[9].sections[0].a[2] = 9.9533079781e-01;
        bands[9].sections[0].z[0] = 0.0;
        bands[9].sections[0].z[1] = 0.0;
        bands[9].sections[1].b[0] = 1.0000000000e+00;
        bands[9].sections[1].b[1] = 2.0000000000e+00;
        bands[9].sections[1].b[2] = 1.0000000000e+00;
        bands[9].sections[1].a[0] = 1.0000000000e+00;
        bands[9].sections[1].a[1] = -1.9953275197e+00;
        bands[9].sections[1].a[2] = 9.9572790736e-01;
        bands[9].sections[1].z[0] = 0.0;
        bands[9].sections[1].z[1] = 0.0;
        bands[9].sections[2].b[0] = 1.0000000000e+00;
        bands[9].sections[2].b[1] = -2.0000000000e+00;
        bands[9].sections[2].b[2] = 1.0000000000e+00;
        bands[9].sections[2].a[0] = 1.0000000000e+00;
        bands[9].sections[2].a[1] = -1.9974059951e+00;
        bands[9].sections[2].a[2] = 9.9794859161e-01;
        bands[9].sections[2].z[0] = 0.0;
        bands[9].sections[2].z[1] = 0.0;
        bands[9].sections[3].b[0] = 1.0000000000e+00;
        bands[9].sections[3].b[1] = -2.0000000000e+00;
        bands[9].sections[3].b[2] = 1.0000000000e+00;
        bands[9].sections[3].a[0] = 1.0000000000e+00;
        bands[9].sections[3].a[1] = -1.9979890530e+00;
        bands[9].sections[3].a[2] = 9.9834298478e-01;
        bands[9].sections[3].z[0] = 0.0;
        bands[9].sections[3].z[1] = 0.0;

        // Band 10: 200.0 Hz
        bands[10].center_freq = 200.0;
        bands[10].num_sections = 4;
        bands[10].sections[0].b[0] = 8.3752740615e-11;
        bands[10].sections[0].b[1] = 1.6750548123e-10;
        bands[10].sections[0].b[2] = 8.3752740615e-11;
        bands[10].sections[0].a[0] = 1.0000000000e+00;
        bands[10].sections[0].a[1] = -1.9934198903e+00;
        bands[10].sections[0].a[2] = 9.9416691834e-01;
        bands[10].sections[0].z[0] = 0.0;
        bands[10].sections[0].z[1] = 0.0;
        bands[10].sections[1].b[0] = 1.0000000000e+00;
        bands[10].sections[1].b[1] = 2.0000000000e+00;
        bands[10].sections[1].b[2] = 1.0000000000e+00;
        bands[10].sections[1].a[0] = 1.0000000000e+00;
        bands[10].sections[1].a[1] = -1.9940374504e+00;
        bands[10].sections[1].a[2] = 9.9466271093e-01;
        bands[10].sections[1].z[0] = 0.0;
        bands[10].sections[1].z[1] = 0.0;
        bands[10].sections[2].b[0] = 1.0000000000e+00;
        bands[10].sections[2].b[1] = -2.0000000000e+00;
        bands[10].sections[2].b[2] = 1.0000000000e+00;
        bands[10].sections[2].a[0] = 1.0000000000e+00;
        bands[10].sections[2].a[1] = -1.9965888532e+00;
        bands[10].sections[2].a[2] = 9.9743642153e-01;
        bands[10].sections[2].z[0] = 0.0;
        bands[10].sections[2].z[1] = 0.0;
        bands[10].sections[3].b[0] = 1.0000000000e+00;
        bands[10].sections[3].b[1] = -2.0000000000e+00;
        bands[10].sections[3].b[2] = 1.0000000000e+00;
        bands[10].sections[3].a[0] = 1.0000000000e+00;
        bands[10].sections[3].a[1] = -1.9973762485e+00;
        bands[10].sections[3].a[2] = 9.9792914324e-01;
        bands[10].sections[3].z[0] = 0.0;
        bands[10].sections[3].z[1] = 0.0;

        // Band 11: 250.0 Hz
        bands[11].center_freq = 250.0;
        bands[11].num_sections = 4;
        bands[11].sections[0].b[0] = 2.0407136844e-10;
        bands[11].sections[0].b[1] = 4.0814273687e-10;
        bands[11].sections[0].b[2] = 2.0407136844e-10;
        bands[11].sections[0].a[0] = 1.0000000000e+00;
        bands[11].sections[0].a[1] = -1.9915476509e+00;
        bands[11].sections[0].a[2] = 9.9271399279e-01;
        bands[11].sections[0].z[0] = 0.0;
        bands[11].sections[0].z[1] = 0.0;
        bands[11].sections[1].b[0] = 1.0000000000e+00;
        bands[11].sections[1].b[1] = 2.0000000000e+00;
        bands[11].sections[1].b[2] = 1.0000000000e+00;
        bands[11].sections[1].a[0] = 1.0000000000e+00;
        bands[11].sections[1].a[1] = -1.9923565002e+00;
        bands[11].sections[1].a[2] = 9.9333279104e-01;
        bands[11].sections[1].z[0] = 0.0;
        bands[11].sections[1].z[1] = 0.0;
        bands[11].sections[2].b[0] = 1.0000000000e+00;
        bands[11].sections[2].b[1] = -2.0000000000e+00;
        bands[11].sections[2].b[2] = 1.0000000000e+00;
        bands[11].sections[2].a[0] = 1.0000000000e+00;
        bands[11].sections[2].a[1] = -1.9954727525e+00;
        bands[11].sections[2].a[2] = 9.9679660155e-01;
        bands[11].sections[2].z[0] = 0.0;
        bands[11].sections[2].z[1] = 0.0;
        bands[11].sections[3].b[0] = 1.0000000000e+00;
        bands[11].sections[3].b[1] = -2.0000000000e+00;
        bands[11].sections[3].b[2] = 1.0000000000e+00;
        bands[11].sections[3].a[0] = 1.0000000000e+00;
        bands[11].sections[3].a[1] = -1.9965484142e+00;
        bands[11].sections[3].a[2] = 9.9741206641e-01;
        bands[11].sections[3].z[0] = 0.0;
        bands[11].sections[3].z[1] = 0.0;

        // Band 12: 315.0 Hz
        bands[12].center_freq = 315.0;
        bands[12].num_sections = 4;
        bands[12].sections[0].b[0] = 5.1303990177e-10;
        bands[12].sections[0].b[1] = 1.0260798035e-09;
        bands[12].sections[0].b[2] = 5.1303990177e-10;
        bands[12].sections[0].a[0] = 1.0000000000e+00;
        bands[12].sections[0].a[1] = -1.9889785537e+00;
        bands[12].sections[0].a[2] = 9.9082838460e-01;
        bands[12].sections[0].z[0] = 0.0;
        bands[12].sections[0].z[1] = 0.0;
        bands[12].sections[1].b[0] = 1.0000000000e+00;
        bands[12].sections[1].b[1] = 2.0000000000e+00;
        bands[12].sections[1].b[2] = 1.0000000000e+00;
        bands[12].sections[1].a[0] = 1.0000000000e+00;
        bands[12].sections[1].a[1] = -1.9900579517e+00;
        bands[12].sections[1].a[2] = 9.9160649777e-01;
        bands[12].sections[1].z[0] = 0.0;
        bands[12].sections[1].z[1] = 0.0;
        bands[12].sections[2].b[0] = 1.0000000000e+00;
        bands[12].sections[2].b[1] = -2.0000000000e+00;
        bands[12].sections[2].b[2] = 1.0000000000e+00;
        bands[12].sections[2].a[0] = 1.0000000000e+00;
        bands[12].sections[2].a[1] = -1.9938647643e+00;
        bands[12].sections[2].a[2] = 9.9596549699e-01;
        bands[12].sections[2].z[0] = 0.0;
        bands[12].sections[2].z[1] = 0.0;
        bands[12].sections[3].b[0] = 1.0000000000e+00;
        bands[12].sections[3].b[1] = -2.0000000000e+00;
        bands[12].sections[3].b[2] = 1.0000000000e+00;
        bands[12].sections[3].a[0] = 1.0000000000e+00;
        bands[12].sections[3].a[1] = -1.9953696180e+00;
        bands[12].sections[3].a[2] = 9.9674023348e-01;
        bands[12].sections[3].z[0] = 0.0;
        bands[12].sections[3].z[1] = 0.0;

        // Band 13: 400.0 Hz
        bands[13].center_freq = 400.0;
        bands[13].num_sections = 4;
        bands[13].sections[0].b[0] = 1.3295202928e-09;
        bands[13].sections[0].b[1] = 2.6590405856e-09;
        bands[13].sections[0].b[2] = 1.3295202928e-09;
        bands[13].sections[0].a[0] = 1.0000000000e+00;
        bands[13].sections[0].a[1] = -1.9853891471e+00;
        bands[13].sections[0].a[2] = 9.8836803756e-01;
        bands[13].sections[0].z[0] = 0.0;
        bands[13].sections[0].z[1] = 0.0;
        bands[13].sections[1].b[0] = 1.0000000000e+00;
        bands[13].sections[1].b[1] = 2.0000000000e+00;
        bands[13].sections[1].b[2] = 1.0000000000e+00;
        bands[13].sections[1].a[0] = 1.0000000000e+00;
        bands[13].sections[1].a[1] = -1.9868594324e+00;
        bands[13].sections[1].a[2] = 9.8935344429e-01;
        bands[13].sections[1].z[0] = 0.0;
        bands[13].sections[1].z[1] = 0.0;
        bands[13].sections[2].b[0] = 1.0000000000e+00;
        bands[13].sections[2].b[1] = -2.0000000000e+00;
        bands[13].sections[2].b[2] = 1.0000000000e+00;
        bands[13].sections[2].a[0] = 1.0000000000e+00;
        bands[13].sections[2].a[1] = -1.9914945951e+00;
        bands[13].sections[2].a[2] = 9.9487981632e-01;
        bands[13].sections[2].z[0] = 0.0;
        bands[13].sections[2].z[1] = 0.0;
        bands[13].sections[3].b[0] = 1.0000000000e+00;
        bands[13].sections[3].b[1] = -2.0000000000e+00;
        bands[13].sections[3].b[2] = 1.0000000000e+00;
        bands[13].sections[3].a[0] = 1.0000000000e+00;
        bands[13].sections[3].a[1] = -1.9936533056e+00;
        bands[13].sections[3].a[2] = 9.9586229328e-01;
        bands[13].sections[3].z[0] = 0.0;
        bands[13].sections[3].z[1] = 0.0;

        // Band 14: 500.0 Hz
        bands[14].center_freq = 500.0;
        bands[14].num_sections = 4;
        bands[14].sections[0].b[0] = 3.2331582922e-09;
        bands[14].sections[0].b[1] = 6.4663165844e-09;
        bands[14].sections[0].b[2] = 3.2331582922e-09;
        bands[14].sections[0].a[0] = 1.0000000000e+00;
        bands[14].sections[0].a[1] = -1.9808342785e+00;
        bands[14].sections[0].a[2] = 9.8548141443e-01;
        bands[14].sections[0].z[0] = 0.0;
        bands[14].sections[0].z[1] = 0.0;
        bands[14].sections[1].b[0] = 1.0000000000e+00;
        bands[14].sections[1].b[1] = 2.0000000000e+00;
        bands[14].sections[1].b[2] = 1.0000000000e+00;
        bands[14].sections[1].a[0] = 1.0000000000e+00;
        bands[14].sections[1].a[1] = -1.9828178524e+00;
        bands[14].sections[1].a[2] = 9.8670912937e-01;
        bands[14].sections[1].z[0] = 0.0;
        bands[14].sections[1].z[1] = 0.0;
        bands[14].sections[2].b[0] = 1.0000000000e+00;
        bands[14].sections[2].b[1] = -2.0000000000e+00;
        bands[14].sections[2].b[2] = 1.0000000000e+00;
        bands[14].sections[2].a[0] = 1.0000000000e+00;
        bands[14].sections[2].a[1] = -1.9883190566e+00;
        bands[14].sections[2].a[2] = 9.9360424779e-01;
        bands[14].sections[2].z[0] = 0.0;
        bands[14].sections[2].z[1] = 0.0;
        bands[14].sections[3].b[0] = 1.0000000000e+00;
        bands[14].sections[3].b[1] = -2.0000000000e+00;
        bands[14].sections[3].b[2] = 1.0000000000e+00;
        bands[14].sections[3].a[0] = 1.0000000000e+00;
        bands[14].sections[3].a[1] = -1.9913808773e+00;
        bands[14].sections[3].a[2] = 9.9483028062e-01;
        bands[14].sections[3].z[0] = 0.0;
        bands[14].sections[3].z[1] = 0.0;

        // Band 15: 630.0 Hz
        bands[15].center_freq = 630.0;
        bands[15].num_sections = 4;
        bands[15].sections[0].b[0] = 8.1076056950e-09;
        bands[15].sections[0].b[1] = 1.6215211390e-08;
        bands[15].sections[0].b[2] = 8.1076056950e-09;
        bands[15].sections[0].a[0] = 1.0000000000e+00;
        bands[15].sections[0].a[1] = -1.9743792869e+00;
        bands[15].sections[0].a[2] = 9.8174157148e-01;
        bands[15].sections[0].z[0] = 0.0;
        bands[15].sections[0].z[1] = 0.0;
        bands[15].sections[1].b[0] = 1.0000000000e+00;
        bands[15].sections[1].b[1] = 2.0000000000e+00;
        bands[15].sections[1].b[2] = 1.0000000000e+00;
        bands[15].sections[1].a[0] = 1.0000000000e+00;
        bands[15].sections[1].a[1] = -1.9771156395e+00;
        bands[15].sections[1].a[2] = 9.8328164323e-01;
        bands[15].sections[1].z[0] = 0.0;
        bands[15].sections[1].z[1] = 0.0;
        bands[15].sections[2].b[0] = 1.0000000000e+00;
        bands[15].sections[2].b[1] = -2.0000000000e+00;
        bands[15].sections[2].b[2] = 1.0000000000e+00;
        bands[15].sections[2].a[0] = 1.0000000000e+00;
        bands[15].sections[2].a[1] = -1.9835671939e+00;
        bands[15].sections[2].a[2] = 9.9194883501e-01;
        bands[15].sections[2].z[0] = 0.0;
        bands[15].sections[2].z[1] = 0.0;
        bands[15].sections[3].b[0] = 1.0000000000e+00;
        bands[15].sections[3].b[1] = -2.0000000000e+00;
        bands[15].sections[3].b[2] = 1.0000000000e+00;
        bands[15].sections[3].a[0] = 1.0000000000e+00;
        bands[15].sections[3].a[1] = -1.9880183216e+00;
        bands[15].sections[3].a[2] = 9.9348999484e-01;
        bands[15].sections[3].z[0] = 0.0;
        bands[15].sections[3].z[1] = 0.0;

        // Band 16: 800.0 Hz
        bands[16].center_freq = 800.0;
        bands[16].num_sections = 4;
        bands[16].sections[0].b[0] = 2.0941054136e-08;
        bands[16].sections[0].b[1] = 4.1882108272e-08;
        bands[16].sections[0].b[2] = 2.0941054136e-08;
        bands[16].sections[0].a[0] = 1.0000000000e+00;
        bands[16].sections[0].a[1] = -1.9650345671e+00;
        bands[16].sections[0].a[2] = 9.7687277199e-01;
        bands[16].sections[0].z[0] = 0.0;
        bands[16].sections[0].z[1] = 0.0;
        bands[16].sections[1].b[0] = 1.0000000000e+00;
        bands[16].sections[1].b[1] = 2.0000000000e+00;
        bands[16].sections[1].b[2] = 1.0000000000e+00;
        bands[16].sections[1].a[0] = 1.0000000000e+00;
        bands[16].sections[1].a[1] = -1.9688992907e+00;
        bands[16].sections[1].a[2] = 9.7881656097e-01;
        bands[16].sections[1].z[0] = 0.0;
        bands[16].sections[1].z[1] = 0.0;
        bands[16].sections[2].b[0] = 1.0000000000e+00;
        bands[16].sections[2].b[1] = -2.0000000000e+00;
        bands[16].sections[2].b[2] = 1.0000000000e+00;
        bands[16].sections[2].a[0] = 1.0000000000e+00;
        bands[16].sections[2].a[1] = -1.9762940882e+00;
        bands[16].sections[2].a[2] = 9.8978904502e-01;
        bands[16].sections[2].z[0] = 0.0;
        bands[16].sections[2].z[1] = 0.0;
        bands[16].sections[3].b[0] = 1.0000000000e+00;
        bands[16].sections[3].b[1] = -2.0000000000e+00;
        bands[16].sections[3].b[2] = 1.0000000000e+00;
        bands[16].sections[3].a[0] = 1.0000000000e+00;
        bands[16].sections[3].a[1] = -1.9829266057e+00;
        bands[16].sections[3].a[2] = 9.9173946151e-01;
        bands[16].sections[3].z[0] = 0.0;
        bands[16].sections[3].z[1] = 0.0;

        // Band 17: 1000.0 Hz
        bands[17].center_freq = 1000.0;
        bands[17].num_sections = 4;
        bands[17].sections[0].b[0] = 5.0727813667e-08;
        bands[17].sections[0].b[1] = 1.0145562733e-07;
        bands[17].sections[0].b[2] = 5.0727813667e-08;
        bands[17].sections[0].a[0] = 1.0000000000e+00;
        bands[17].sections[0].a[1] = -1.9527423270e+00;
        bands[17].sections[0].a[2] = 9.7117632899e-01;
        bands[17].sections[0].z[0] = 0.0;
        bands[17].sections[0].z[1] = 0.0;
        bands[17].sections[1].b[0] = 1.0000000000e+00;
        bands[17].sections[1].b[1] = 2.0000000000e+00;
        bands[17].sections[1].b[2] = 1.0000000000e+00;
        bands[17].sections[1].a[0] = 1.0000000000e+00;
        bands[17].sections[1].a[1] = -1.9581399076e+00;
        bands[17].sections[1].a[2] = 9.7358775957e-01;
        bands[17].sections[1].z[0] = 0.0;
        bands[17].sections[1].z[1] = 0.0;
        bands[17].sections[2].b[0] = 1.0000000000e+00;
        bands[17].sections[2].b[1] = -2.0000000000e+00;
        bands[17].sections[2].b[2] = 1.0000000000e+00;
        bands[17].sections[2].a[0] = 1.0000000000e+00;
        bands[17].sections[2].a[1] = -1.9662099097e+00;
        bands[17].sections[2].a[2] = 9.8725563041e-01;
        bands[17].sections[2].z[0] = 0.0;
        bands[17].sections[2].z[1] = 0.0;
        bands[17].sections[3].b[0] = 1.0000000000e+00;
        bands[17].sections[3].b[1] = -2.0000000000e+00;
        bands[17].sections[3].b[2] = 1.0000000000e+00;
        bands[17].sections[3].a[0] = 1.0000000000e+00;
        bands[17].sections[3].a[1] = -1.9759327057e+00;
        bands[17].sections[3].a[2] = 9.8968290755e-01;
        bands[17].sections[3].z[0] = 0.0;
        bands[17].sections[3].z[1] = 0.0;

        // Band 18: 1250.0 Hz
        bands[18].center_freq = 1250.0;
        bands[18].num_sections = 4;
        bands[18].sections[0].b[0] = 1.2264796379e-07;
        bands[18].sections[0].b[1] = 2.4529592757e-07;
        bands[18].sections[0].b[2] = 1.2264796379e-07;
        bands[18].sections[0].a[0] = 1.0000000000e+00;
        bands[18].sections[0].a[1] = -1.9354279749e+00;
        bands[18].sections[0].a[2] = 9.6410370396e-01;
        bands[18].sections[0].z[0] = 0.0;
        bands[18].sections[0].z[1] = 0.0;
        bands[18].sections[1].b[0] = 1.0000000000e+00;
        bands[18].sections[1].b[1] = 2.0000000000e+00;
        bands[18].sections[1].b[2] = 1.0000000000e+00;
        bands[18].sections[1].a[0] = 1.0000000000e+00;
        bands[18].sections[1].a[1] = -1.9430468640e+00;
        bands[18].sections[1].a[2] = 9.6708776378e-01;
        bands[18].sections[1].z[0] = 0.0;
        bands[18].sections[1].z[1] = 0.0;
        bands[18].sections[2].b[0] = 1.0000000000e+00;
        bands[18].sections[2].b[1] = -2.0000000000e+00;
        bands[18].sections[2].b[2] = 1.0000000000e+00;
        bands[18].sections[2].a[0] = 1.0000000000e+00;
        bands[18].sections[2].a[1] = -1.9513014962e+00;
        bands[18].sections[2].a[2] = 9.8410082611e-01;
        bands[18].sections[2].z[0] = 0.0;
        bands[18].sections[2].z[1] = 0.0;
        bands[18].sections[3].b[0] = 1.0000000000e+00;
        bands[18].sections[3].b[1] = -2.0000000000e+00;
        bands[18].sections[3].b[2] = 1.0000000000e+00;
        bands[18].sections[3].a[0] = 1.0000000000e+00;
        bands[18].sections[3].a[1] = -1.9656730433e+00;
        bands[18].sections[3].a[2] = 9.8711620638e-01;
        bands[18].sections[3].z[0] = 0.0;
        bands[18].sections[3].z[1] = 0.0;

        // Band 19: 1600.0 Hz
        bands[19].center_freq = 1600.0;
        bands[19].num_sections = 4;
        bands[19].sections[0].b[0] = 3.2479650500e-07;
        bands[19].sections[0].b[1] = 6.4959301000e-07;
        bands[19].sections[0].b[2] = 3.2479650500e-07;
        bands[19].sections[0].a[0] = 1.0000000000e+00;
        bands[19].sections[0].a[1] = -1.9076132574e+00;
        bands[19].sections[0].a[2] = 9.5429138292e-01;
        bands[19].sections[0].z[0] = 0.0;
        bands[19].sections[0].z[1] = 0.0;
        bands[19].sections[1].b[0] = 1.0000000000e+00;
        bands[19].sections[1].b[1] = 2.0000000000e+00;
        bands[19].sections[1].b[2] = 1.0000000000e+00;
        bands[19].sections[1].a[0] = 1.0000000000e+00;
        bands[19].sections[1].a[1] = -1.9188940539e+00;
        bands[19].sections[1].a[2] = 9.5805299174e-01;
        bands[19].sections[1].z[0] = 0.0;
        bands[19].sections[1].z[1] = 0.0;
        bands[19].sections[2].b[0] = 1.0000000000e+00;
        bands[19].sections[2].b[1] = -2.0000000000e+00;
        bands[19].sections[2].b[2] = 1.0000000000e+00;
        bands[19].sections[2].a[0] = 1.0000000000e+00;
        bands[19].sections[2].a[1] = -1.9261821256e+00;
        bands[19].sections[2].a[2] = 9.7970772245e-01;
        bands[19].sections[2].z[0] = 0.0;
        bands[19].sections[2].z[1] = 0.0;
        bands[19].sections[3].b[0] = 1.0000000000e+00;
        bands[19].sections[3].b[1] = -2.0000000000e+00;
        bands[19].sections[3].b[2] = 1.0000000000e+00;
        bands[19].sections[3].a[0] = 1.0000000000e+00;
        bands[19].sections[3].a[1] = -1.9485002511e+00;
        bands[19].sections[3].a[2] = 9.8352928341e-01;
        bands[19].sections[3].z[0] = 0.0;
        bands[19].sections[3].z[1] = 0.0;

        // Band 20: 2000.0 Hz
        bands[20].center_freq = 2000.0;
        bands[20].num_sections = 4;
        bands[20].sections[0].b[0] = 7.8083643709e-07;
        bands[20].sections[0].b[1] = 1.5616728742e-06;
        bands[20].sections[0].b[2] = 7.8083643709e-07;
        bands[20].sections[0].a[0] = 1.0000000000e+00;
        bands[20].sections[0].a[1] = -1.8708395354e+00;
        bands[20].sections[0].a[2] = 9.4320474652e-01;
        bands[20].sections[0].z[0] = 0.0;
        bands[20].sections[0].z[1] = 0.0;
        bands[20].sections[1].b[0] = 1.0000000000e+00;
        bands[20].sections[1].b[1] = 2.0000000000e+00;
        bands[20].sections[1].b[2] = 1.0000000000e+00;
        bands[20].sections[1].a[0] = 1.0000000000e+00;
        bands[20].sections[1].a[1] = -1.8870605258e+00;
        bands[20].sections[1].a[2] = 9.4781723028e-01;
        bands[20].sections[1].z[0] = 0.0;
        bands[20].sections[1].z[1] = 0.0;
        bands[20].sections[2].b[0] = 1.0000000000e+00;
        bands[20].sections[2].b[1] = -2.0000000000e+00;
        bands[20].sections[2].b[2] = 1.0000000000e+00;
        bands[20].sections[2].a[0] = 1.0000000000e+00;
        bands[20].sections[2].a[1] = -1.8915106784e+00;
        bands[20].sections[2].a[2] = 9.7472306463e-01;
        bands[20].sections[2].z[0] = 0.0;
        bands[20].sections[2].z[1] = 0.0;
        bands[20].sections[3].b[0] = 1.0000000000e+00;
        bands[20].sections[3].b[1] = -2.0000000000e+00;
        bands[20].sections[3].b[2] = 1.0000000000e+00;
        bands[20].sections[3].a[0] = 1.0000000000e+00;
        bands[20].sections[3].a[1] = -1.9249074526e+00;
        bands[20].sections[3].a[2] = 9.7943740106e-01;
        bands[20].sections[3].z[0] = 0.0;
        bands[20].sections[3].z[1] = 0.0;

        // Band 21: 2500.0 Hz
        bands[21].center_freq = 2500.0;
        bands[21].num_sections = 4;
        bands[21].sections[0].b[0] = 1.8702246496e-06;
        bands[21].sections[0].b[1] = 3.7404492992e-06;
        bands[21].sections[0].b[2] = 1.8702246496e-06;
        bands[21].sections[0].a[0] = 1.0000000000e+00;
        bands[21].sections[0].a[1] = -1.8176383212e+00;
        bands[21].sections[0].a[2] = 9.2953724538e-01;
        bands[21].sections[0].z[0] = 0.0;
        bands[21].sections[0].z[1] = 0.0;
        bands[21].sections[1].b[0] = 1.0000000000e+00;
        bands[21].sections[1].b[1] = 2.0000000000e+00;
        bands[21].sections[1].b[2] = 1.0000000000e+00;
        bands[21].sections[1].a[0] = 1.0000000000e+00;
        bands[21].sections[1].a[1] = -1.8411002801e+00;
        bands[21].sections[1].a[2] = 9.3515068421e-01;
        bands[21].sections[1].z[0] = 0.0;
        bands[21].sections[1].z[1] = 0.0;
        bands[21].sections[2].b[0] = 1.0000000000e+00;
        bands[21].sections[2].b[1] = -2.0000000000e+00;
        bands[21].sections[2].b[2] = 1.0000000000e+00;
        bands[21].sections[2].a[0] = 1.0000000000e+00;
        bands[21].sections[2].a[1] = -1.8394492824e+00;
        bands[21].sections[2].a[2] = 9.6855054962e-01;
        bands[21].sections[2].z[0] = 0.0;
        bands[21].sections[2].z[1] = 0.0;
        bands[21].sections[3].b[0] = 1.0000000000e+00;
        bands[21].sections[3].b[1] = -2.0000000000e+00;
        bands[21].sections[3].b[2] = 1.0000000000e+00;
        bands[21].sections[3].a[0] = 1.0000000000e+00;
        bands[21].sections[3].a[1] = -1.8895659971e+00;
        bands[21].sections[3].a[2] = 9.7433041351e-01;
        bands[21].sections[3].z[0] = 0.0;
        bands[21].sections[3].z[1] = 0.0;

        // Band 22: 3150.0 Hz
        bands[22].center_freq = 3150.0;
        bands[22].num_sections = 4;
        bands[22].sections[0].b[0] = 4.5991007162e-06;
        bands[22].sections[0].b[1] = 9.1982014324e-06;
        bands[22].sections[0].b[2] = 4.5991007162e-06;
        bands[22].sections[0].a[0] = 1.0000000000e+00;
        bands[22].sections[0].a[1] = -1.7370051632e+00;
        bands[22].sections[0].a[2] = 9.1208572767e-01;
        bands[22].sections[0].z[0] = 0.0;
        bands[22].sections[0].z[1] = 0.0;
        bands[22].sections[1].b[0] = 1.0000000000e+00;
        bands[22].sections[1].b[1] = 2.0000000000e+00;
        bands[22].sections[1].b[2] = 1.0000000000e+00;
        bands[22].sections[1].a[0] = 1.0000000000e+00;
        bands[22].sections[1].a[1] = -1.7715029014e+00;
        bands[22].sections[1].a[2] = 9.1888423312e-01;
        bands[22].sections[1].z[0] = 0.0;
        bands[22].sections[1].z[1] = 0.0;
        bands[22].sections[2].b[0] = 1.0000000000e+00;
        bands[22].sections[2].b[1] = -2.0000000000e+00;
        bands[22].sections[2].b[2] = 1.0000000000e+00;
        bands[22].sections[2].a[0] = 1.0000000000e+00;
        bands[22].sections[2].a[1] = -1.7578103425e+00;
        bands[22].sections[2].a[2] = 9.6063188158e-01;
        bands[22].sections[2].z[0] = 0.0;
        bands[22].sections[2].z[1] = 0.0;
        bands[22].sections[3].b[0] = 1.0000000000e+00;
        bands[22].sections[3].b[1] = -2.0000000000e+00;
        bands[22].sections[3].b[2] = 1.0000000000e+00;
        bands[22].sections[3].a[0] = 1.0000000000e+00;
        bands[22].sections[3].a[1] = -1.8341411742e+00;
        bands[22].sections[3].a[2] = 9.6769749730e-01;
        bands[22].sections[3].z[0] = 0.0;
        bands[22].sections[3].z[1] = 0.0;

        // Band 23: 4000.0 Hz
        bands[23].center_freq = 4000.0;
        bands[23].num_sections = 4;
        bands[23].sections[0].b[0] = 1.1583372579e-05;
        bands[23].sections[0].b[1] = 2.3166745159e-05;
        bands[23].sections[0].b[2] = 1.1583372579e-05;
        bands[23].sections[0].a[0] = 1.0000000000e+00;
        bands[23].sections[0].a[1] = -1.6133028848e+00;
        bands[23].sections[0].a[2] = 8.8980273776e-01;
        bands[23].sections[0].z[0] = 0.0;
        bands[23].sections[0].z[1] = 0.0;
        bands[23].sections[1].b[0] = 1.0000000000e+00;
        bands[23].sections[1].b[1] = 2.0000000000e+00;
        bands[23].sections[1].b[2] = 1.0000000000e+00;
        bands[23].sections[1].a[0] = 1.0000000000e+00;
        bands[23].sections[1].a[1] = -1.6646511710e+00;
        bands[23].sections[1].a[2] = 8.9792644163e-01;
        bands[23].sections[1].z[0] = 0.0;
        bands[23].sections[1].z[1] = 0.0;
        bands[23].sections[2].b[0] = 1.0000000000e+00;
        bands[23].sections[2].b[1] = -2.0000000000e+00;
        bands[23].sections[2].b[2] = 1.0000000000e+00;
        bands[23].sections[2].a[0] = 1.0000000000e+00;
        bands[23].sections[2].a[1] = -1.6285888068e+00;
        bands[23].sections[2].a[2] = 9.5047616516e-01;
        bands[23].sections[2].z[0] = 0.0;
        bands[23].sections[2].z[1] = 0.0;
        bands[23].sections[3].b[0] = 1.0000000000e+00;
        bands[23].sections[3].b[1] = -2.0000000000e+00;
        bands[23].sections[3].b[2] = 1.0000000000e+00;
        bands[23].sections[3].a[0] = 1.0000000000e+00;
        bands[23].sections[3].a[1] = -1.7461000097e+00;
        bands[23].sections[3].a[2] = 9.5901872221e-01;
        bands[23].sections[3].z[0] = 0.0;
        bands[23].sections[3].z[1] = 0.0;

        // Band 24: 5000.0 Hz
        bands[24].center_freq = 5000.0;
        bands[24].num_sections = 4;
        bands[24].sections[0].b[0] = 2.7253846053e-05;
        bands[24].sections[0].b[1] = 5.4507692106e-05;
        bands[24].sections[0].b[2] = 2.7253846053e-05;
        bands[24].sections[0].a[0] = 1.0000000000e+00;
        bands[24].sections[0].a[1] = -1.4438706954e+00;
        bands[24].sections[0].a[2] = 8.6436759582e-01;
        bands[24].sections[0].z[0] = 0.0;
        bands[24].sections[0].z[1] = 0.0;
        bands[24].sections[1].b[0] = 1.0000000000e+00;
        bands[24].sections[1].b[1] = 2.0000000000e+00;
        bands[24].sections[1].b[2] = 1.0000000000e+00;
        bands[24].sections[1].a[0] = 1.0000000000e+00;
        bands[24].sections[1].a[1] = -1.5178801537e+00;
        bands[24].sections[1].a[2] = 8.7367654899e-01;
        bands[24].sections[1].z[0] = 0.0;
        bands[24].sections[1].z[1] = 0.0;
        bands[24].sections[2].b[0] = 1.0000000000e+00;
        bands[24].sections[2].b[1] = -2.0000000000e+00;
        bands[24].sections[2].b[2] = 1.0000000000e+00;
        bands[24].sections[2].a[0] = 1.0000000000e+00;
        bands[24].sections[2].a[1] = -1.4467940547e+00;
        bands[24].sections[2].a[2] = 9.3885347125e-01;
        bands[24].sections[2].z[0] = 0.0;
        bands[24].sections[2].z[1] = 0.0;
        bands[24].sections[3].b[0] = 1.0000000000e+00;
        bands[24].sections[3].b[1] = -2.0000000000e+00;
        bands[24].sections[3].b[2] = 1.0000000000e+00;
        bands[24].sections[3].a[0] = 1.0000000000e+00;
        bands[24].sections[3].a[1] = -1.6212340474e+00;
        bands[24].sections[3].a[2] = 9.4877199428e-01;
        bands[24].sections[3].z[0] = 0.0;
        bands[24].sections[3].z[1] = 0.0;

        // Band 25: 6300.0 Hz
        bands[25].center_freq = 6300.0;
        bands[25].num_sections = 4;
        bands[25].sections[0].b[0] = 6.5526052261e-05;
        bands[25].sections[0].b[1] = 1.3105210452e-04;
        bands[25].sections[0].b[2] = 6.5526052261e-05;
        bands[25].sections[0].a[0] = 1.0000000000e+00;
        bands[25].sections[0].a[1] = -1.1907460788e+00;
        bands[25].sections[0].a[2] = 8.3256622652e-01;
        bands[25].sections[0].z[0] = 0.0;
        bands[25].sections[0].z[1] = 0.0;
        bands[25].sections[1].b[0] = 1.0000000000e+00;
        bands[25].sections[1].b[1] = 2.0000000000e+00;
        bands[25].sections[1].b[2] = 1.0000000000e+00;
        bands[25].sections[1].a[0] = 1.0000000000e+00;
        bands[25].sections[1].a[1] = -1.2973596890e+00;
        bands[25].sections[1].a[2] = 8.4270962842e-01;
        bands[25].sections[1].z[0] = 0.0;
        bands[25].sections[1].z[1] = 0.0;
        bands[25].sections[2].b[0] = 1.0000000000e+00;
        bands[25].sections[2].b[1] = -2.0000000000e+00;
        bands[25].sections[2].b[2] = 1.0000000000e+00;
        bands[25].sections[2].a[0] = 1.0000000000e+00;
        bands[25].sections[2].a[1] = -1.1689637311e+00;
        bands[25].sections[2].a[2] = 9.2434615335e-01;
        bands[25].sections[2].z[0] = 0.0;
        bands[25].sections[2].z[1] = 0.0;
        bands[25].sections[3].b[0] = 1.0000000000e+00;
        bands[25].sections[3].b[1] = -2.0000000000e+00;
        bands[25].sections[3].b[2] = 1.0000000000e+00;
        bands[25].sections[3].a[0] = 1.0000000000e+00;
        bands[25].sections[3].a[1] = -1.4275907538e+00;
        bands[25].sections[3].a[2] = 9.3532837875e-01;
        bands[25].sections[3].z[0] = 0.0;
        bands[25].sections[3].z[1] = 0.0;

        // Band 26: 8000.0 Hz
        bands[26].center_freq = 8000.0;
        bands[26].num_sections = 4;
        bands[26].sections[0].b[0] = 1.6040574966e-04;
        bands[26].sections[0].b[1] = 3.2081149932e-04;
        bands[26].sections[0].b[2] = 1.6040574966e-04;
        bands[26].sections[0].a[0] = 1.0000000000e+00;
        bands[26].sections[0].a[1] = -8.1719051702e-01;
        bands[26].sections[0].a[2] = 7.9316533363e-01;
        bands[26].sections[0].z[0] = 0.0;
        bands[26].sections[0].z[1] = 0.0;
        bands[26].sections[1].b[0] = 1.0000000000e+00;
        bands[26].sections[1].b[1] = 2.0000000000e+00;
        bands[26].sections[1].b[2] = 1.0000000000e+00;
        bands[26].sections[1].a[0] = 1.0000000000e+00;
        bands[26].sections[1].a[1] = -9.6862161371e-01;
        bands[26].sections[1].a[2] = 8.0294634638e-01;
        bands[26].sections[1].z[0] = 0.0;
        bands[26].sections[1].z[1] = 0.0;
        bands[26].sections[2].b[0] = 1.0000000000e+00;
        bands[26].sections[2].b[1] = -2.0000000000e+00;
        bands[26].sections[2].b[2] = 1.0000000000e+00;
        bands[26].sections[2].a[0] = 1.0000000000e+00;
        bands[26].sections[2].a[1] = -7.5113810411e-01;
        bands[26].sections[2].a[2] = 9.0658661275e-01;
        bands[26].sections[2].z[0] = 0.0;
        bands[26].sections[2].z[1] = 0.0;
        bands[26].sections[3].b[0] = 1.0000000000e+00;
        bands[26].sections[3].b[1] = -2.0000000000e+00;
        bands[26].sections[3].b[2] = 1.0000000000e+00;
        bands[26].sections[3].a[0] = 1.0000000000e+00;
        bands[26].sections[3].a[1] = -1.1287745450e+00;
        bands[26].sections[3].a[2] = 9.1738009123e-01;
        bands[26].sections[3].z[0] = 0.0;
        bands[26].sections[3].z[1] = 0.0;

        // Band 27: 10000.0 Hz
        bands[27].center_freq = 10000.0;
        bands[27].num_sections = 4;
        bands[27].sections[0].b[0] = 3.6551901396e-04;
        bands[27].sections[0].b[1] = 7.3103802793e-04;
        bands[27].sections[0].b[2] = 3.6551901396e-04;
        bands[27].sections[0].a[0] = 1.0000000000e+00;
        bands[27].sections[0].a[1] = -3.4184508429e-01;
        bands[27].sections[0].a[2] = 7.5011325492e-01;
        bands[27].sections[0].z[0] = 0.0;
        bands[27].sections[0].z[1] = 0.0;
        bands[27].sections[1].b[0] = 1.0000000000e+00;
        bands[27].sections[1].b[1] = 2.0000000000e+00;
        bands[27].sections[1].b[2] = 1.0000000000e+00;
        bands[27].sections[1].a[0] = 1.0000000000e+00;
        bands[27].sections[1].a[1] = -5.4362466713e-01;
        bands[27].sections[1].a[2] = 7.5678725377e-01;
        bands[27].sections[1].z[0] = 0.0;
        bands[27].sections[1].z[1] = 0.0;
        bands[27].sections[2].b[0] = 1.0000000000e+00;
        bands[27].sections[2].b[1] = -2.0000000000e+00;
        bands[27].sections[2].b[2] = 1.0000000000e+00;
        bands[27].sections[2].a[0] = 1.0000000000e+00;
        bands[27].sections[2].a[1] = -2.1297356521e-01;
        bands[27].sections[2].a[2] = 8.8786071806e-01;
        bands[27].sections[2].z[0] = 0.0;
        bands[27].sections[2].z[1] = 0.0;
        bands[27].sections[3].b[0] = 1.0000000000e+00;
        bands[27].sections[3].b[1] = -2.0000000000e+00;
        bands[27].sections[3].b[2] = 1.0000000000e+00;
        bands[27].sections[3].a[0] = 1.0000000000e+00;
        bands[27].sections[3].a[1] = -7.2743584771e-01;
        bands[27].sections[3].a[2] = 8.9536884091e-01;
        bands[27].sections[3].z[0] = 0.0;
        bands[27].sections[3].z[1] = 0.0;

        // Band 28: 12500.0 Hz
        bands[28].center_freq = 12500.0;
        bands[28].num_sections = 4;
        bands[28].sections[0].b[0] = 8.2108238703e-04;
        bands[28].sections[0].b[1] = -1.6421647741e-03;
        bands[28].sections[0].b[2] = 8.2108238703e-04;
        bands[28].sections[0].a[0] = 1.0000000000e+00;
        bands[28].sections[0].a[1] = 5.3325488961e-03;
        bands[28].sections[0].a[2] = 6.9895765693e-01;
        bands[28].sections[0].z[0] = 0.0;
        bands[28].sections[0].z[1] = 0.0;
        bands[28].sections[1].b[0] = 1.0000000000e+00;
        bands[28].sections[1].b[1] = 2.0000000000e+00;
        bands[28].sections[1].b[2] = 1.0000000000e+00;
        bands[28].sections[1].a[0] = 1.0000000000e+00;
        bands[28].sections[1].a[1] = 2.5766083001e-01;
        bands[28].sections[1].a[2] = 7.0190923815e-01;
        bands[28].sections[1].z[0] = 0.0;
        bands[28].sections[1].z[1] = 0.0;
        bands[28].sections[2].b[0] = 1.0000000000e+00;
        bands[28].sections[2].b[1] = -2.0000000000e+00;
        bands[28].sections[2].b[2] = 1.0000000000e+00;
        bands[28].sections[2].a[0] = 1.0000000000e+00;
        bands[28].sections[2].a[1] = -1.8477308461e-01;
        bands[28].sections[2].a[2] = 8.6549145081e-01;
        bands[28].sections[2].z[0] = 0.0;
        bands[28].sections[2].z[1] = 0.0;
        bands[28].sections[3].b[0] = 1.0000000000e+00;
        bands[28].sections[3].b[1] = 2.0000000000e+00;
        bands[28].sections[3].b[2] = 1.0000000000e+00;
        bands[28].sections[3].a[0] = 1.0000000000e+00;
        bands[28].sections[3].a[1] = 4.6663550942e-01;
        bands[28].sections[3].a[2] = 8.6887107446e-01;
        bands[28].sections[3].z[0] = 0.0;
        bands[28].sections[3].z[1] = 0.0;

        // Band 29: 16000.0 Hz
        bands[29].center_freq = 16000.0;
        bands[29].num_sections = 4;
        bands[29].sections[0].b[0] = 1.9719201682e-03;
        bands[29].sections[0].b[1] = -3.9438403365e-03;
        bands[29].sections[0].b[2] = 1.9719201682e-03;
        bands[29].sections[0].a[0] = 1.0000000000e+00;
        bands[29].sections[0].a[1] = 7.0843611868e-01;
        bands[29].sections[0].a[2] = 6.1282354090e-01;
        bands[29].sections[0].z[0] = 0.0;
        bands[29].sections[0].z[1] = 0.0;
        bands[29].sections[1].b[0] = 1.0000000000e+00;
        bands[29].sections[1].b[1] = -2.0000000000e+00;
        bands[29].sections[1].b[2] = 1.0000000000e+00;
        bands[29].sections[1].a[0] = 1.0000000000e+00;
        bands[29].sections[1].a[1] = 9.9892865902e-01;
        bands[29].sections[1].a[2] = 6.4930453542e-01;
        bands[29].sections[1].z[0] = 0.0;
        bands[29].sections[1].z[1] = 0.0;
        bands[29].sections[2].b[0] = 1.0000000000e+00;
        bands[29].sections[2].b[1] = 2.0000000000e+00;
        bands[29].sections[2].b[2] = 1.0000000000e+00;
        bands[29].sections[2].a[0] = 1.0000000000e+00;
        bands[29].sections[2].a[1] = 5.6042306540e-01;
        bands[29].sections[2].a[2] = 8.1439287683e-01;
        bands[29].sections[2].z[0] = 0.0;
        bands[29].sections[2].z[1] = 0.0;
        bands[29].sections[3].b[0] = 1.0000000000e+00;
        bands[29].sections[3].b[1] = 2.0000000000e+00;
        bands[29].sections[3].b[2] = 1.0000000000e+00;
        bands[29].sections[3].a[0] = 1.0000000000e+00;
        bands[29].sections[3].a[1] = 1.2866162001e+00;
        bands[29].sections[3].a[2] = 8.5610905393e-01;
        bands[29].sections[3].z[0] = 0.0;
        bands[29].sections[3].z[1] = 0.0;

        // Band 30: 20000.0 Hz
        bands[30].center_freq = 20000.0;
        bands[30].num_sections = 4;
        bands[30].sections[0].b[0] = 4.2701518998e-03;
        bands[30].sections[0].b[1] = -8.5403037995e-03;
        bands[30].sections[0].b[2] = 4.2701518998e-03;
        bands[30].sections[0].a[0] = 1.0000000000e+00;
        bands[30].sections[0].a[1] = 1.1979947229e+00;
        bands[30].sections[0].a[2] = 4.5134001086e-01;
        bands[30].sections[0].z[0] = 0.0;
        bands[30].sections[0].z[1] = 0.0;
        bands[30].sections[1].b[0] = 1.0000000000e+00;
        bands[30].sections[1].b[1] = -2.0000000000e+00;
        bands[30].sections[1].b[2] = 1.0000000000e+00;
        bands[30].sections[1].a[0] = 1.0000000000e+00;
        bands[30].sections[1].a[1] = 1.6216153437e+00;
        bands[30].sections[1].a[2] = 6.8744181415e-01;
        bands[30].sections[1].z[0] = 0.0;
        bands[30].sections[1].z[1] = 0.0;
        bands[30].sections[2].b[0] = 1.0000000000e+00;
        bands[30].sections[2].b[1] = 2.0000000000e+00;
        bands[30].sections[2].b[2] = 1.0000000000e+00;
        bands[30].sections[2].a[0] = 1.0000000000e+00;
        bands[30].sections[2].a[1] = 1.2058440678e+00;
        bands[30].sections[2].a[2] = 7.0929749752e-01;
        bands[30].sections[2].z[0] = 0.0;
        bands[30].sections[2].z[1] = 0.0;
        bands[30].sections[3].b[0] = 1.0000000000e+00;
        bands[30].sections[3].b[1] = 2.0000000000e+00;
        bands[30].sections[3].b[2] = 1.0000000000e+00;
        bands[30].sections[3].a[0] = 1.0000000000e+00;
        bands[30].sections[3].a[1] = 1.8635789025e+00;
        bands[30].sections[3].a[2] = 9.0508965329e-01;
        bands[30].sections[3].z[0] = 0.0;
        bands[30].sections[3].z[1] = 0.0;

    }
}

void resetFilterState() {
    // Clear delay lines for all bands and sections
    for (int band = 0; band < NUM_BANDS; band++) {
        for (int section = 0; section < pFilterBank->bands[band].num_sections; section++) {
            pFilterBank->bands[band].sections[section].z[0] = 0.0;
            pFilterBank->bands[band].sections[section].z[1] = 0.0;
        }
    }
    pFilterBank->initialized = 1;
}

void calculateSamplesIntegration() {
    pFilterBank->maxNumberOfSamples = (uint64_t)((pFilterBank->sample_rate / 1000) * (float)pFilterBank->integrationTime);
}

void calculateLevel() {
    // Calculate the level for each band based on the accumulated temporal sum
    for (int band = 0; band < NUM_BANDS; band++) {
        if (pFilterBank->samplesCount > 0) {
            // Calculate RMS level
            float rms = sqrt(pFilterBank->temporal_sum[band] / (float)pFilterBank->samplesCount);
            //Calculate dB level
            float level_db = 10.0f * log10f(rms) + pFilterBank->mic_constant;
            // Smooth the instantaneous level
            if (pFilterBank->smoothed_level[band] == 0) {
                pFilterBank->smoothed_level[band] = pFilterBank->volume_level[band]; // Initialize on first pass
            } else {
                pFilterBank->smoothed_level[band] = pFilterBank->alpha * pFilterBank->volume_level[band] + (1 - pFilterBank->alpha) * pFilterBank->smoothed_level[band];
            }
            pFilterBank->volume_level[band] = pFilterBank->smoothed_level[band];
            // Reset temporal sum for the next integration period
            pFilterBank->temporal_sum[band] = 0.0f;
            // Reset sample count after processing
            pFilterBank->samplesCount = 0;
        }
    }
    pFilterBank->samplesCount = 0; // Reset sample count after processing
}

void third_octave_filter_alloc(float sample_rate, int filter_order = 2) {
    if (pFilterBank == NULL) {
        pFilterBank = (struct third_octave_filter *) malloc(sizeof(struct third_octave_filter));
    }

    pFilterBank->bypass = 1; // Filter disabled by default
    pFilterBank->sample_rate = sample_rate;
    pFilterBank->filter_order = (filter_order == 2 || filter_order == 4) ? filter_order : 2; // Default to order 2
    pFilterBank->initialized = 0;
    pFilterBank->channelType = LEFT_CHANNEL; // Default to left channel

    // Initialize filter coefficients for the specified order
    initializeFilterBank(pFilterBank->bands, pFilterBank->filter_order);
    resetFilterState();

    pFilterBank->mic_constant = 120.0f; // Default microphone constant in dB
    pFilterBank->integrationTime = 125; // Default integration time in milliseconds
    pFilterBank->samplesCount = 0;
    // Initialize temporal sum and volume leven in dBSPL to zero for all bands
    for (int i = 0; i < NUM_BANDS; ++i) {
        pFilterBank->temporal_sum[i] = 0.0f;
        pFilterBank->volume_level[i] = 0.0f; //Principal array for storing volume levels in dBSPL
        pFilterBank->smoothed_level[i] = 0;
    }
    pFilterBank->maxNumberOfSamples = (uint64_t)((pFilterBank->sample_rate / 1000) * pFilterBank->integrationTime);

    calculateSamplesIntegration();
}

void third_octave_filter_free() {
    if (pFilterBank != NULL) {
        free(pFilterBank);
        pFilterBank = NULL;
    }
}

// Updated processing function using cascaded biquads
void third_octave_filter_process(float *data, int buffer_size) {
    const int byPass = pFilterBank->bypass == 1 ? 1 : 0;
    if (byPass) return;

    double buffer_sum = 0;
    int samples_processed = 0;

    // Process samples
    switch (pFilterBank->channelType) {
        case LEFT_CHANNEL:
            for (int i = 0; i < buffer_size; i += 2) {
                auto input = (double)data[i];
                double output = 0.0;

                // Process through all bands and accumulate the output for temporal integration
                for (int band = 0; band < NUM_BANDS; band++) {
                    output = processFilterBand(&pFilterBank->bands[band], input);
                    // Accumulate the output for temporal integration
                    pFilterBank->temporal_sum[band] += (float)pow(output, 2);

                }

                pFilterBank->samplesCount += 1;

                // If the number of samples processed reaches the maximum, calculate the level
                if (pFilterBank->samplesCount >= pFilterBank->maxNumberOfSamples) {
                    calculateLevel();
                }
            }
            break;

        case RIGHT_CHANNEL:
            for (int i = 1; i < buffer_size; i += 2) {
                auto input = (double)data[i];
                double output = 0.0;

                // Process through all bands and accumulate the output for temporal integration
                for (int band = 0; band < NUM_BANDS; band++) {
                    output = processFilterBand(&pFilterBank->bands[band], input);
                    // Accumulate the output for temporal integration
                    pFilterBank->temporal_sum[band] += (float)pow(output, 2);

                }

                pFilterBank->samplesCount += 1;

                // If the number of samples processed reaches the maximum, calculate the level
                if (pFilterBank->samplesCount >= pFilterBank->maxNumberOfSamples) {
                    calculateLevel();
                }
            }
            break;

        case STEREO_CHANNEL:
            for (int i = 0; i < buffer_size; i += 2) {
                auto input_left = (double)data[i];
                auto input_right = (double)data[i+1];
                double output_left = 0.0, output_right = 0.0;

                for (int band = 0; band < NUM_BANDS; band++) {
                    output_left = processFilterBand(&pFilterBank->bands[band], input_left);
                    output_right = processFilterBand(&pFilterBank->bands[band], input_right);
                    // Accumulate the output for temporal integration
                    pFilterBank->temporal_sum[band] += (float)pow(output_left, 2);
                    pFilterBank->temporal_sum[band] += (float)pow(output_right, 2);
                }

                pFilterBank->samplesCount += 2; // Two samples processed (left and right)

                // If the number of samples processed reaches the maximum, calculate the level
                if (pFilterBank->samplesCount >= pFilterBank->maxNumberOfSamples) {
                    calculateLevel();
                }
            }
            break;

        default:
            break;
    }
}