
/**
 * @file bz_oscillation.c
 * @brief Oscillation detection and analysis: peak counting, period estimation,
 *        Fourier analysis, phase plane, Poincare section.
 *
 * Implements:
 *   L2: Oscillation detection, peak/trough counting
 *   L3: Phase plane trajectory, Poincare section
 *   L5: DFT-based period analysis, Hanning window
 *   L6: Oscillation classification (sinusoidal/relaxation/bursting)
 */

#include "bz_oscillation.h"
#include "bz_reaction.h"
#include <math.h>
#include <string.h>

/* =========================================================================
 * L2: Oscillation detection
 * ========================================================================= */

int bz_detect_oscillation(const double *data, int n, double threshold)
{
    /* Detect oscillation by counting sign changes of the derivative.
     * Requires at least 4 sign changes (2 full cycles) for sustained oscillation.
     * threshold is the minimum amplitude of the overall signal to consider. */
    int i, crossings;
    double deriv, prev_deriv;
    double data_min, data_max;

    if (!data || n < 5) return 0;

    /* First, check overall amplitude */
    data_min = data_max = data[0];
    for (i = 1; i < n; i++) {
        if (data[i] < data_min) data_min = data[i];
        if (data[i] > data_max) data_max = data[i];
    }
    if (data_max - data_min < threshold) return 0;

    /* Count sign changes of the derivative */
    crossings = 0;
    prev_deriv = data[1] - data[0];

    for (i = 2; i < n; i++) {
        deriv = data[i] - data[i-1];
        if (prev_deriv * deriv < 0.0) {
            crossings++;
        }
        if (fabs(deriv) > 1e-15) {
            prev_deriv = deriv;
        }
    }

    return (crossings >= 4) ? 1 : 0;
}

int bz_count_peaks(const double *data, int n, int *peaks, int max_peaks)
{
    /* Find local maxima (peaks):
     *   data[i-1] < data[i] > data[i+1] */
    int i, count = 0;

    if (!data || n < 3 || !peaks || max_peaks <= 0) return 0;

    for (i = 1; i < n - 1 && count < max_peaks; i++) {
        if (data[i] > data[i-1] && data[i] > data[i+1]) {
            peaks[count] = i;
            count++;
        }
    }
    return count;
}

int bz_estimate_period(const int *peaks, int n_peaks, double dt,
                        double *period, double *std_dev)
{
    /* Estimate average period and standard deviation from peak intervals. */
    double sum, sum_sq, mean;
    int i, valid = 0;

    if (!peaks || n_peaks < 2) {
        if (period) *period = 0.0;
        if (std_dev) *std_dev = 0.0;
        return 0;
    }

    sum = 0.0;
    for (i = 1; i < n_peaks; i++) {
        double interval = (peaks[i] - peaks[i-1]) * dt;
        sum += interval;
        valid++;
    }

    if (valid == 0) {
        if (period) *period = 0.0;
        if (std_dev) *std_dev = 0.0;
        return 0;
    }

    mean = sum / (double)valid;
    if (period) *period = mean;

    if (std_dev && valid > 1) {
        sum_sq = 0.0;
        for (i = 1; i < n_peaks; i++) {
            double interval = (peaks[i] - peaks[i-1]) * dt;
            double diff = interval - mean;
            sum_sq += diff * diff;
        }
        *std_dev = sqrt(sum_sq / (double)(valid - 1));
    } else if (std_dev) {
        *std_dev = 0.0;
    }

    return valid;
}

/* =========================================================================
 * L3: Phase plane and Poincare section
 * ========================================================================= */

int bz_phase_plane_trajectory(double *traj_x, double *traj_z,
                               const bz_oregonator_state_t *s0,
                               const bz_oregonator_params_t *p,
                               int n_steps, double dt)
{
    /* Integrate Oregonator and record (x,z) projection. */
    bz_oregonator_state_t s, dsdt;
    int i;

    if (!traj_x || !traj_z || !s0 || !p || n_steps <= 0) return 0;

    s = *s0;
    traj_x[0] = s.x;
    traj_z[0] = s.z;

    for (i = 1; i < n_steps; i++) {
        bz_oregonator_rhs(&dsdt, &s, p);
        s.x += dt * dsdt.x;
        s.y += dt * dsdt.y;
        s.z += dt * dsdt.z;
        traj_x[i] = s.x;
        traj_z[i] = s.z;
    }

    return n_steps;
}

int bz_poincare_section(double *crossings_z, int max_cross,
                         const double *traj_x, const double *traj_z __attribute__((unused)),
                         int n, double x_section)
{
    /* Find crossings of x = x_section with dx/dt > 0.
     * A crossing is detected when traj_x[i-1] < x_section < traj_x[i]. */
    int cross = 0;
    int i;

    if (!crossings_z || !traj_x || !traj_z || n < 2) return 0;

    for (i = 1; i < n && cross < max_cross; i++) {
        double x_prev = traj_x[i-1];
        double x_curr = traj_x[i];

        if (x_prev < x_section && x_curr >= x_section) {
            /* Linear interpolation for z at crossing */
            double frac = (x_section - x_prev) / (x_curr - x_prev);
            if (frac >= 0.0 && frac <= 1.0) {
                crossings_z[cross] = traj_z[i-1] + frac * (traj_z[i] - traj_z[i-1]);
                cross++;
            }
        }
    }
    return cross;
}

/* =========================================================================
 * L5: Fourier analysis
 * ========================================================================= */

void bz_dft_magnitude(double *magnitude, const double *data, int n,
                       int *n_freq)
{
    /* Compute magnitude spectrum via DFT.
     * X[k] = sum_{n=0}^{N-1} x[n] * exp(-2*pi*i*k*n/N)
     * magnitude[k] = |X[k]| */
    int k, j;
    int n_out;
    double pi = 3.141592653589793;
    double omega;

    if (!magnitude || !data || n <= 0) {
        if (n_freq) *n_freq = 0;
        return;
    }

    n_out = n / 2 + 1;
    if (n_freq) *n_freq = n_out;

    for (k = 0; k < n_out; k++) {
        double re = 0.0, im = 0.0;
        omega = -2.0 * pi * k / (double)n;
        for (j = 0; j < n; j++) {
            double angle = omega * j;
            re += data[j] * cos(angle);
            im += data[j] * sin(angle);
        }
        magnitude[k] = sqrt(re * re + im * im);
    }
}

int bz_dominant_frequency(double *freq, const double *magnitude, int n_freq)
{
    /* Find k with maximum magnitude (skip DC, k=0).
     * Dominant frequency = k / N (cycles per sample). */
    int k, max_k = 0;
    double max_val = 0.0;

    if (!freq || !magnitude || n_freq <= 1) return 0;

    for (k = 1; k < n_freq; k++) {
        if (magnitude[k] > max_val) {
            max_val = magnitude[k];
            max_k = k;
        }
    }

    if (max_k == 0) {
        *freq = 0.0;
        return 0;
    }

    *freq = (double)max_k / (double)((n_freq - 1) * 2);
    return max_k;
}

void bz_hanning_window(double *data, int n)
{
    /* Hanning window: w[j] = 0.5*(1 - cos(2*pi*j/(N-1))) */
    int j;
    double pi = 3.141592653589793;
    double denom;

    if (!data || n <= 1) return;

    denom = (double)(n - 1);
    for (j = 0; j < n; j++) {
        double w = 0.5 * (1.0 - cos(2.0 * pi * j / denom));
        data[j] *= w;
    }
}

/* =========================================================================
 * L6: Oscillation characterization
 * ========================================================================= */

int bz_classify_oscillation(const double *traj_x, const double *traj_z __attribute__((unused)),
                             int n, double dt __attribute__((unused)))
{
    /* Classify oscillation type:
     *   0: steady state
     *   1: sinusoidal (near Hopf onset)
     *   2: relaxation (sharp spikes)
     *   3: bursting (mixed timescales) */

    double x_min, x_max, x_mean, x_var, duty;
    int peaks[100], n_peaks;
    int i;

    if (!traj_x || n < 10) return 0;

    /* Compute basic statistics */
    x_min = x_max = traj_x[0];
    x_mean = 0.0;
    for (i = 0; i < n; i++) {
        x_mean += traj_x[i];
        if (traj_x[i] < x_min) x_min = traj_x[i];
        if (traj_x[i] > x_max) x_max = traj_x[i];
    }
    x_mean /= (double)n;

    /* Check if oscillatory */
    if (x_max - x_min < 0.01) return 0;  /* steady state */

    /* Count peaks */
    n_peaks = bz_count_peaks(traj_x, n, peaks, 100);

    if (n_peaks < 2) return 0;  /* transients */

    /* Compute duty cycle (fraction above mean) */
    {
        int above = 0;
        for (i = 0; i < n; i++) {
            if (traj_x[i] > x_mean) above++;
        }
        duty = (double)above / (double)n;
    }

    /* Compute variance */
    x_var = 0.0;
    for (i = 0; i < n; i++) {
        double d = traj_x[i] - x_mean;
        x_var += d * d;
    }
    x_var /= (double)n;

    /* Classification logic */
    if (duty < 0.2 || duty > 0.8) {
        return 2;  /* relaxation: asymmetric duty cycle */
    } else if (x_var / (x_max - x_min + 1e-10) < 0.1) {
        return 1;  /* sinusoidal: low variance relative to amplitude */
    } else {
        return 3;  /* bursting */
    }
}

double bz_oscillation_amplitude(const double *data, int n,
                                 double *peak, double *trough)
{
    double min_val, max_val;
    int i;

    if (!data || n <= 0) {
        if (peak) *peak = 0.0;
        if (trough) *trough = 0.0;
        return 0.0;
    }

    min_val = max_val = data[0];
    for (i = 1; i < n; i++) {
        if (data[i] < min_val) min_val = data[i];
        if (data[i] > max_val) max_val = data[i];
    }

    if (peak) *peak = max_val;
    if (trough) *trough = min_val;
    return max_val - min_val;
}

double bz_duty_cycle(const double *data, int n, double mid)
{
    int above = 0;
    int i;

    if (!data || n <= 0) return 0.0;

    for (i = 0; i < n; i++) {
        if (data[i] > mid) above++;
    }
    return (double)above / (double)n;
}
