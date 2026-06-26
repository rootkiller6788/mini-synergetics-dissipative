/**
 * @file test_oscillation.c
 * @brief Tests for oscillation detection, Fourier analysis, phase plane.
 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../include/bz_oscillation.h"
#include "../include/bz_reaction.h"

static int test_detect_oscillation_sine(void) {
    int i, n = 200;
    double *data = malloc(n * sizeof(double));
    for (i = 0; i < n; i++) data[i] = sin(2.0*3.14159*i/20.0);
    int result = bz_detect_oscillation(data, n, 0.1);
    free(data);
    assert(result == 1);
    return 1;
}

static int test_detect_oscillation_constant(void) {
    int i, n = 100;
    double *data = malloc(n * sizeof(double));
    for (i = 0; i < n; i++) data[i] = 1.0;
    int result = bz_detect_oscillation(data, n, 0.1);
    free(data);
    assert(result == 0);
    return 1;
}

static int test_count_peaks(void) {
    int i, n = 100;
    double *data = malloc(n * sizeof(double));
    int peaks[50], n_peaks;
    for (i = 0; i < n; i++) data[i] = sin(2.0*3.14159*i/10.0);
    n_peaks = bz_count_peaks(data, n, peaks, 50);
    free(data);
    assert(n_peaks >= 5);
    assert(n_peaks <= 15);
    return 1;
}

static int test_estimate_period(void) {
    int n_peaks = 5;
    int peaks[] = {10, 30, 50, 70, 90};
    double period, std_dev;
    int valid = bz_estimate_period(peaks, n_peaks, 0.1, &period, &std_dev);
    assert(valid == 4);
    assert(fabs(period - 2.0) < 0.01);
    return 1;
}

static int test_dft_magnitude(void) {
    int i, n = 64, n_freq;
    double *data = malloc(n * sizeof(double));
    double *mag;
    for (i = 0; i < n; i++) data[i] = sin(2.0*3.14159*i*2.0/n);
    mag = malloc((n/2+1) * sizeof(double));
    bz_dft_magnitude(mag, data, n, &n_freq);
    assert(n_freq == n/2 + 1);
    /* Peak should be at bin 2 */
    assert(mag[2] > mag[1]);
    assert(mag[2] > mag[3]);
    free(data);
    free(mag);
    return 1;
}

static int test_dominant_frequency(void) {
    int i, n = 64, n_freq;
    double *data = malloc(n * sizeof(double));
    double *mag;
    double freq;
    int idx;
    for (i = 0; i < n; i++) data[i] = sin(2.0*3.14159*i*3.0/n);
    mag = malloc((n/2+1) * sizeof(double));
    bz_dft_magnitude(mag, data, n, &n_freq);
    idx = bz_dominant_frequency(&freq, mag, n_freq);
    assert(idx > 0);
    assert(freq > 0.0);
    free(data);
    free(mag);
    return 1;
}

static int test_hanning_window(void) {
    int i, n = 64;
    double *data = malloc(n * sizeof(double));
    for (i = 0; i < n; i++) data[i] = 1.0;
    bz_hanning_window(data, n);
    assert(fabs(data[0]) < 1e-10);
    assert(fabs(data[n-1]) < 1e-10);
    assert(data[n/2] > 0.9);
    free(data);
    return 1;
}

static int test_oscillation_amplitude(void) {
    double data[] = {0.0, 1.0, 0.5, 2.0, -1.0, 0.5};
    double peak, trough, amp;
    amp = bz_oscillation_amplitude(data, 6, &peak, &trough);
    assert(fabs(peak - 2.0) < 1e-10);
    assert(fabs(trough + 1.0) < 1e-10);
    assert(fabs(amp - 3.0) < 1e-10);
    return 1;
}

static int test_duty_cycle(void) {
    int i, n = 100;
    double *data = malloc(n * sizeof(double));
    for (i = 0; i < n; i++) data[i] = (i < 30) ? 1.0 : 0.0;
    double duty = bz_duty_cycle(data, n, 0.5);
    free(data);
    assert(fabs(duty - 0.3) < 0.05);
    return 1;
}

static int test_classify_oscillation(void) {
    int i, n = 200;
    double *traj_x = malloc(n * sizeof(double));
    double *traj_z = malloc(n * sizeof(double));
    for (i = 0; i < n; i++) {
        traj_x[i] = sin(2.0*3.14159*i/20.0);
        traj_z[i] = cos(2.0*3.14159*i/20.0);
    }
    int osc_type = bz_classify_oscillation(traj_x, traj_z, n, 0.01);
    assert(osc_type >= 0 && osc_type <= 3);
    free(traj_x);
    free(traj_z);
    return 1;
}

int main(void) {
    int passed = 0, total = 0;
#define RUN(t) do { total++; printf("  %s... ", #t); \
    if (t()) { printf("PASS\n"); passed++; } else { printf("FAIL\n"); } \
} while(0)

    printf("=== Oscillation Tests ===\n");
    RUN(test_detect_oscillation_sine);
    RUN(test_detect_oscillation_constant);
    RUN(test_count_peaks);
    RUN(test_estimate_period);
    RUN(test_dft_magnitude);
    RUN(test_dominant_frequency);
    RUN(test_hanning_window);
    RUN(test_oscillation_amplitude);
    RUN(test_duty_cycle);
    RUN(test_classify_oscillation);

    printf("\n%d/%d tests passed\n", passed, total);
    return (passed == total) ? 0 : 1;
}
