/**
 * @file bz_oscillation.h
 * @brief Oscillation detection, period analysis, phase plane,
 *        Fourier spectrum, and Poincare section for BZ reaction.
 *
 * Knowledge coverage:
 *   L2: Oscillation onset, limit cycles, period estimation
 *   L3: Phase plane, Poincare map, Fourier decomposition
 *   L5: Peak/trough detection, FFT-based period analysis
 *   L6: Oregonator oscillation characterization
 */

#ifndef BZ_OSCILLATION_H
#define BZ_OSCILLATION_H

#include "bz_reaction.h"

/* ---------------------------------------------------------------------------
 * L2: Oscillation detection and characterization
 * --------------------------------------------------------------------------- */

/** Detect if a time series is oscillatory by counting zero-crossings
 *  of the derivative. Returns 1 if sustained oscillation detected.
 *  @param data     Time series array
 *  @param n        Length of time series
 *  @param threshold Minimum amplitude for a crossing to count
 *  @return         1 if oscillatory, 0 otherwise */
int bz_detect_oscillation(const double *data, int n, double threshold);

/** Count oscillation peaks in a time series.
 *  A peak is a local maximum with amplitude above the threshold.
 *  @param data      Time series
 *  @param n         Number of points
 *  @param peaks     Output array of peak indices (caller allocates >= n/2)
 *  @param max_peaks Maximum peaks to find
 *  @return          Number of peaks found */
int bz_count_peaks(const double *data, int n, int *peaks, int max_peaks);

/** Estimate oscillation period from peak spacing.
 *  @param peaks     Array of peak indices
 *  @param n_peaks   Number of peaks
 *  @param dt        Time step between samples
 *  @param period    Output: average period
 *  @param std_dev   Output: standard deviation of period
 *  @return          Number of valid inter-peak intervals used */
int bz_estimate_period(const int *peaks, int n_peaks, double dt,
                        double *period, double *std_dev);

/* ---------------------------------------------------------------------------
 * L3: Phase plane analysis
 * --------------------------------------------------------------------------- */

/** Compute phase plane trajectory from Oregonator integration.
 *  Stores (x, z) projection of the 3D trajectory.
 *  @param traj_x   Output x trajectory
 *  @param traj_z   Output z trajectory
 *  @param s0       Initial state
 *  @param p        Parameters
 *  @param n_steps  Number of integration steps
 *  @param dt       Time step
 *  @return         Number of steps successfully integrated */
int bz_phase_plane_trajectory(double *traj_x, double *traj_z,
                               const bz_oregonator_state_t *s0,
                               const bz_oregonator_params_t *p,
                               int n_steps, double dt);

/** Detect limit cycle crossing of Poincare section (x = x_section).
 *  Records crossings with dx/dt > 0 only.
 *  @param crossings_z   Output: z values at crossing
 *  @param max_cross     Max crossings to record
 *  @param traj_x, traj_z Trajectory data
 *  @param n              Number of points
 *  @param x_section      Poincare section x value
 *  @return               Number of crossings found */
int bz_poincare_section(double *crossings_z, int max_cross,
                         const double *traj_x, const double *traj_z,
                         int n, double x_section);

/* ---------------------------------------------------------------------------
 * L5: Fourier analysis of oscillations
 * --------------------------------------------------------------------------- */

/** Compute discrete Fourier transform (real-valued, magnitude spectrum).
 *  Simple DFT: X[k] = sum_n x[n] * exp(-2*pi*i*k*n/N)
 *  @param magnitude  Output: |X[k]| for k=0..N/2
 *  @param data       Input time series
 *  @param n          Number of points
 *  @param n_freq     Output: number of frequency bins (n/2+1) */
void bz_dft_magnitude(double *magnitude, const double *data, int n,
                       int *n_freq);

/** Find dominant frequency from DFT magnitude spectrum.
 *  @param freq      Output: dominant frequency (cycles per sample)
 *  @param magnitude DFT magnitude spectrum
 *  @param n_freq    Number of frequency bins
 *  @return          Index of maximum magnitude (excluding DC) */
int bz_dominant_frequency(double *freq, const double *magnitude, int n_freq);

/** Apply Hanning window to a time series.
 *  w[n] = 0.5 * (1 - cos(2*pi*n/(N-1)))
 *  @param data   Time series (modified in-place)
 *  @param n      Number of points */
void bz_hanning_window(double *data, int n);

/* ---------------------------------------------------------------------------
 * L6: Oscillation characterization
 * --------------------------------------------------------------------------- */

/** Characterize Oregonator oscillation type.
 *  Determines if the oscillation is:
 *    0: steady state (no oscillation)
 *    1: sinusoidal (near Hopf onset)
 *    2: relaxation (large spikes with slow recovery)
 *    3: bursting (mixed time scales)
 *  @param traj_x    x trajectory
 *  @param traj_z    z trajectory
 *  @param n         Number of points
 *  @param dt        Time step
 *  @return          Oscillation type (0-3) */
int bz_classify_oscillation(const double *traj_x, const double *traj_z,
                             int n, double dt);

/** Compute oscillation amplitude (peak-to-trough).
 *  @param data   Time series
 *  @param n      Number of points
 *  @param peak   Output: maximum value
 *  @param trough Output: minimum value
 *  @return       Amplitude (peak - trough) */
double bz_oscillation_amplitude(const double *data, int n,
                                 double *peak, double *trough);

/** Compute duty cycle: fraction of period spent above mid-point.
 *  @param data   Time series (one period)
 *  @param n      Points in one period
 *  @param mid    Mid-point value
 *  @return       Duty cycle (0 to 1) */
double bz_duty_cycle(const double *data, int n, double mid);

#endif /* BZ_OSCILLATION_H */
