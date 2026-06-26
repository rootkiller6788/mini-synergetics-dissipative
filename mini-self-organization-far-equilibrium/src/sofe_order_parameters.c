#include "sofe_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

double sofe_landau_free_energy(double psi, double a, double b, double Tc,
                                double T) {
    double r = a * (T - Tc);
    double F0 = 0.0;
    return F0 + 0.5 * r * psi * psi + 0.25 * b * psi * psi * psi * psi;
}

double sofe_landau_equilibrium_amplitude(double a, double b, double Tc,
                                          double T) {
    double r = a * (T - Tc);
    if (r >= 0.0 || fabs(b) < 1e-15) return 0.0;
    return sqrt(-r / b);
}

void sofe_slaving_eliminate(double y_u, double* y_s,
                             double (*g)(double), double tau,
                             double* y_u_dot, double* y_s_dot) {
    double ys_eq = g(y_u);
    *y_s = ys_eq;
    *y_s_dot = -(ys_eq - ys_eq) / tau;
    (void)y_u_dot;
}

double sofe_critical_exponent_beta(double* temperatures, double* amplitudes,
                                     int n, double Tc) {
    if (!temperatures || !amplitudes || n < 3) return 0.5;
    double sum_x = 0.0, sum_y = 0.0, sum_xx = 0.0, sum_xy = 0.0;
    int count = 0;
    for (int i = 0; i < n; i++) {
        double dt = Tc - temperatures[i];
        if (dt > 1e-15 && amplitudes[i] > 1e-15) {
            double x = log(dt);
            double y = log(amplitudes[i]);
            sum_x += x; sum_y += y;
            sum_xx += x * x; sum_xy += x * y;
            count++;
        }
    }
    if (count < 2) return 0.5;
    return (count * sum_xy - sum_x * sum_y) /
           (count * sum_xx - sum_x * sum_x + 1e-15);
}

double sofe_susceptibility(double psi, double a, double b, double Tc,
                            double T) {
    (void)psi; (void)b;
    double r = a * (T - Tc);
    if (r <= 0.0) return 1.0 / (2.0 * fabs(r) + 1e-15);
    return 1.0 / (r + 1e-15);
}

double sofe_correlation_length_mean_field(double xi0, double Tc, double T) {
    double t = (T - Tc) / Tc;
    if (fabs(t) < 1e-15) return xi0 * 1e6;
    return xi0 / sqrt(fabs(t));
}

double sofe_lmc_complexity(double* probabilities, int n) {
    if (!probabilities || n <= 0) return 0.0;
    double H = 0.0;
    for (int i = 0; i < n; i++) {
        if (probabilities[i] > 1e-15) H -= probabilities[i] * log(probabilities[i]);
    }
    double H_max = log((double)n);
    double H_norm = (H_max > 1e-15) ? H / H_max : 0.0;
    double D = 0.0;
    double uniform = 1.0 / n;
    for (int i = 0; i < n; i++) {
        double diff = probabilities[i] - uniform;
        D += diff * diff;
    }
    return H_norm * D;
}

double sofe_effective_measure_complexity(double* time_series, int n,
                                           int block_size) {
    if (!time_series || n < block_size || block_size <= 0) return 0.0;
    int n_blocks = n / block_size;
    if (n_blocks <= 1) return 0.0;
    double* block_means = (double*)malloc(n_blocks * sizeof(double));
    for (int i = 0; i < n_blocks; i++) {
        double sum = 0.0;
        for (int j = 0; j < block_size; j++)
            sum += time_series[i * block_size + j];
        block_means[i] = sum / block_size;
    }
    double mean = 0.0;
    for (int i = 0; i < n_blocks; i++) mean += block_means[i];
    mean /= n_blocks;
    double variance = 0.0;
    for (int i = 0; i < n_blocks; i++)
        variance += (block_means[i] - mean) * (block_means[i] - mean);
    variance /= n_blocks;
    free(block_means);
    return variance;
}

double sofe_order_parameter_correlation_time(double* timeseries, int n,
                                               double dt) {
    if (!timeseries || n < 2 || dt <= 0.0) return 0.0;
    double mean = 0.0;
    for (int i = 0; i < n; i++) mean += timeseries[i];
    mean /= n;
    double* fluct = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) fluct[i] = timeseries[i] - mean;
    double C0 = 0.0;
    for (int i = 0; i < n; i++) C0 += fluct[i] * fluct[i];
    C0 /= n;
    if (fabs(C0) < 1e-15) { free(fluct); return 0.0; }
    double integral = C0;
    int max_lag = n / 4;
    for (int lag = 1; lag < max_lag; lag++) {
        double C = 0.0;
        for (int i = 0; i < n - lag; i++)
            C += fluct[i + lag] * fluct[i];
        C /= (n - lag);
        integral += C;
        if (C < 0.01 * C0) break;
    }
    free(fluct);
    return integral * dt / C0;
}
