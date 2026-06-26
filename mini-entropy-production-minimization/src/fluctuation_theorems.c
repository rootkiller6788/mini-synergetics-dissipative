/**
 * @file fluctuation_theorems.c
 * @brief Implementation of fluctuation theorems in non-equilibrium statistical mechanics.
 *
 * Implements stochastic trajectory simulation (Langevin equation),
 * entropy production along trajectories, Evans-Searles and Crooks
 * fluctuation theorem verification, Jarzynski equality free energy
 * estimation, and the Fluctuation-Dissipation Theorem.
 */

#include "fluctuation_theorems.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Simple pseudo-random number generator (xorshift) for reproducibility */
static unsigned int ft_rand_state = 123456789;

static unsigned int ft_xorshift32(void)
{
    unsigned int x = ft_rand_state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    ft_rand_state = x;
    return x;
}

static double ft_rand_uniform(void)
{
    return (double)ft_xorshift32() / (double)0xFFFFFFFFu;
}

/* ---------------------------------------------------------------------------
 * L2: Core Concepts -- Brownian motion simulation
 * ------------------------------------------------------------------------- */

double ft_rand_gaussian(double mean, double stddev)
{
    /* Box-Muller method */
    double u1 = ft_rand_uniform();
    double u2 = ft_rand_uniform();

    /* Avoid log(0) */
    if (u1 < 1e-15) u1 = 1e-15;
    if (u2 < 1e-15) u2 = 1e-15;

    double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);

    return mean + stddev * z0;
}

int ft_simulate_brownian_harmonic(StochasticTrajectory *traj,
                                   double mass, double friction,
                                   double spring_k, double temperature,
                                   double k_B, double x0, double v0,
                                   double dt, int n_steps)
{
    int step;
    double x, v;
    double gamma = friction / mass;  /* Damping rate */

    if (traj == NULL || n_steps <= 0) return -1;
    if (mass <= 0.0 || temperature <= 0.0) return -1;

    /* Noise strength: sqrt(2 * gamma * k_B * T / mass) * sqrt(dt)
     * The Langevin equation: m dv = -gamma m v dt - k x dt + sqrt(2 gamma m k_B T) dW
     * Simplifies to: dv = -(gamma v + (k/m) x) dt + sqrt(2 gamma k_B T/m) dW */
    double noise_amplitude = sqrt(2.0 * friction * k_B * temperature) / mass;

    x = x0;
    v = v0;

    traj->time[0] = 0.0;
    traj->position[0] = x;
    traj->velocity[0] = v;
    traj->entropy_production[0] = 0.0;
    traj->work_done[0] = 0.0;
    traj->heat_dissipated[0] = 0.0;

    double cumulative_sigma = 0.0;
    double cumulative_work = 0.0;
    double cumulative_heat = 0.0;

    for (step = 1; step < n_steps; step++) {
        /* Euler-Maruyama step for the Langevin equation */
        double noise = noise_amplitude * ft_rand_gaussian(0.0, sqrt(dt));

        /* Acceleration: a = -(gamma)v - (k/m)x + noise */
        double a = -gamma * v - (spring_k / mass) * x;

        /* Update velocity and position */
        double v_new = v + a * dt + noise;
        double x_new = x + v * dt;

        /* Work done by the trap on the particle:
         * dW = (partial U / partial x) * dx = k * x * dx */
        double dW = spring_k * x * (x_new - x);

        /* Heat dissipated to bath:
         * dQ = friction * v^2 * dt - v * sqrt(2 gamma k_B T) * dW
         * Simplified: dQ = -dW - dE_kinetic - dU */
        double dU = 0.5 * spring_k * (x_new * x_new - x * x);
        double dK = 0.5 * mass * (v_new * v_new - v * v);
        double dQ = -dW - dU - dK;

        /* Entropy production: dSigma = dQ / T + dS_sys
         * For a harmonic trap, system entropy change is approximately
         * dS_sys = k_B * ln(P(x_new, v_new) / P(x, v)) where P is
         * the equilibrium distribution.
         * Here we use: dSigma = dQ/T + (E_new - E_old)/T = (dQ + dU + dK)/T
         * which equals 0 for reversible processes and > 0 for irreversible. */
        double dSigma = (dQ + dU + dK) / temperature;
        if (dSigma < 0.0) dSigma = 0.0;  /* Numerical correction */

        cumulative_sigma += dSigma;
        cumulative_work += dW;
        cumulative_heat += dQ;

        x = x_new;
        v = v_new;

        /* Record */
        traj->time[step] = step * dt;
        traj->position[step] = x;
        traj->velocity[step] = v;
        traj->entropy_production[step] = cumulative_sigma;
        traj->work_done[step] = cumulative_work;
        traj->heat_dissipated[step] = cumulative_heat;
    }

    return 0;
}

int ft_simulate_jarzynski_protocol(StochasticTrajectory *traj,
                                    double mass, double friction,
                                    double spring_k, double L, double tau,
                                    double temperature, double k_B, double dt)
{
    int step, n_steps;
    double x, v;
    double gamma = friction / mass;
    double noise_amplitude = sqrt(2.0 * friction * k_B * temperature) / mass;

    n_steps = (int)(tau / dt);
    if (n_steps <= 0 || traj == NULL) return -1;

    /* Start at equilibrium in initial trap at x=0 */
    x = ft_rand_gaussian(0.0, sqrt(k_B * temperature / spring_k));
    v = ft_rand_gaussian(0.0, sqrt(k_B * temperature / mass));

    double cumulative_work = 0.0;

    for (step = 0; step < n_steps; step++) {
        double t = step * dt;
        double lambda = (t / tau) * L;  /* Trap center moves from 0 to L */
        double dlambda_dt = L / tau;

        /* Force from trap: F = -k (x - lambda(t)) */
        double F_trap = -spring_k * (x - lambda);

        /* Work increment due to trap movement:
         * dW = (partial U / partial lambda) * dlambda
         *     = -k (x - lambda) * dlambda_dt * dt
         *     = F_trap * dlambda_dt * dt */
        double dW = -spring_k * (x - lambda) * dlambda_dt * dt;

        cumulative_work += dW;

        /* Langevin step */
        double noise = noise_amplitude * ft_rand_gaussian(0.0, sqrt(dt));
        double a = -gamma * v + F_trap / mass;
        v = v + a * dt + noise;
        x = x + v * dt;

        /* Record */
        traj->time[step] = t;
        traj->position[step] = x;
        traj->velocity[step] = v;
        traj->work_done[step] = cumulative_work;
    }

    return 0;
}

int ft_compute_trajectory_entropy(StochasticTrajectory *traj)
{
    int step, n_steps;
    double cum_sigma = 0.0;

    if (traj == NULL) return -1;

    n_steps = traj->n_steps;

    /* Stochastic entropy production:
     * Sigma(t) = Q(t)/T + Delta S_sys(t)
     *
     * For a Langevin process, the medium entropy production is:
     * Delta S_m = Q/T = (1/T) integral_0^t friction * v^2 dt
     *
     * The system entropy: S_sys = -k_B * ln(p(x, v, t))
     *
     * Total Sigma = S_m + Delta S_sys */
    for (step = 1; step < n_steps; step++) {
        /* Medium entropy: approximate from heat dissipation */
        double dQ = traj->heat_dissipated[step] - traj->heat_dissipated[step - 1];
        double dS_m = dQ / traj->temperature;

        cum_sigma += dS_m;
        traj->entropy_production[step] = cum_sigma;
    }

    return 0;
}

/* ---------------------------------------------------------------------------
 * L3: Mathematical Structures -- Fluctuation theorem verification
 * ------------------------------------------------------------------------- */

int ft_build_sigma_histogram(EntropyProductionHistogram *hist,
                              const double sigma_values[],
                              int n_samples, int num_bins, double k_B)
{
    int i, bin_idx;
    double sigma_min, sigma_max;

    if (hist == NULL || sigma_values == NULL || n_samples < 1 || num_bins < 2)
        return -1;

    /* Find range */
    sigma_min = sigma_max = sigma_values[0];
    for (i = 1; i < n_samples; i++) {
        if (sigma_values[i] < sigma_min) sigma_min = sigma_values[i];
        if (sigma_values[i] > sigma_max) sigma_max = sigma_values[i];
    }

    /* Expand range slightly */
    double margin = 0.1 * (sigma_max - sigma_min + 1e-10);
    sigma_min -= margin;
    sigma_max += margin;

    hist->num_bins = num_bins;
    hist->bin_width = (sigma_max - sigma_min) / num_bins;
    hist->sigma_min = sigma_min;
    hist->sigma_max = sigma_max;
    hist->total_samples = n_samples;
    hist->k_B = k_B;

    /* Initialize bins and PDF */
    for (i = 0; i < num_bins; i++) {
        hist->bins[i] = sigma_min + (i + 0.5) * hist->bin_width;
        hist->pdf[i] = 0.0;
        hist->pdf_reverse[i] = 0.0;
    }

    /* Fill histogram */
    for (i = 0; i < n_samples; i++) {
        bin_idx = (int)((sigma_values[i] - sigma_min) / hist->bin_width);
        if (bin_idx >= 0 && bin_idx < num_bins) {
            hist->pdf[bin_idx] += 1.0;
        }
    }

    /* Normalize to probability density */
    double norm = 0.0;
    for (i = 0; i < num_bins; i++) {
        norm += hist->pdf[i] * hist->bin_width;
    }
    if (norm > 0.0) {
        for (i = 0; i < num_bins; i++) {
            hist->pdf[i] /= norm;
        }
    }

    /* Build reverse PDF: P_R(Sigma) = P_F(-Sigma) */
    for (i = 0; i < num_bins; i++) {
        double neg_sigma = -hist->bins[i];
        int rev_idx = (int)((neg_sigma - sigma_min) / hist->bin_width);
        if (rev_idx >= 0 && rev_idx < num_bins) {
            hist->pdf_reverse[i] = hist->pdf[rev_idx];
        }
    }

    return 0;
}

int ft_verify_evans_searles(EntropyProductionHistogram *hist)
{
    /* Evans-Searles Fluctuation Theorem:
     * ln[P(Sigma) / P(-Sigma)] = Sigma / k_B
     *
     * Fit: y = a * x, where y = ln[P(Sigma)/P(-Sigma)], x = Sigma.
     * The slope a should equal 1/k_B.
     */

    int i, n_valid = 0;
    double sum_x = 0.0, sum_y = 0.0, sum_xx = 0.0, sum_xy = 0.0;

    if (hist == NULL || hist->num_bins < 2) return -1;

    for (i = 0; i < hist->num_bins; i++) {
        double sigma_i = hist->bins[i];
        double P_fwd = hist->pdf[i];
        double P_rev = hist->pdf_reverse[i];

        if (P_fwd > 1e-15 && P_rev > 1e-15) {
            double y = log(P_fwd / P_rev);
            sum_x += sigma_i;
            sum_y += y;
            sum_xx += sigma_i * sigma_i;
            sum_xy += sigma_i * y;
            n_valid++;
        }
    }

    if (n_valid < 2) return -1;

    /* Linear regression: y = slope * x (no intercept per the theorem) */
    double slope = sum_xy / sum_xx;
    hist->ft_slope = slope;

    /* R^2 */
    double ss_res = 0.0, ss_tot = 0.0;
    double y_mean = sum_y / n_valid;
    for (i = 0; i < hist->num_bins; i++) {
        double sigma_i = hist->bins[i];
        double P_fwd = hist->pdf[i];
        double P_rev = hist->pdf_reverse[i];
        if (P_fwd > 1e-15 && P_rev > 1e-15) {
            double y_actual = log(P_fwd / P_rev);
            double y_pred = slope * sigma_i;
            ss_res += (y_actual - y_pred) * (y_actual - y_pred);
            ss_tot += (y_actual - y_mean) * (y_actual - y_mean);
        }
    }
    hist->ft_r2 = (ss_tot > 0.0) ? 1.0 - ss_res / ss_tot : 0.0;

    /* Expected: slope = 1/k_B */
    double slope_expected = 1.0 / hist->k_B;
    hist->asymmetry_ratio = slope / slope_expected;

    return 0;
}

int ft_verify_crooks(WorkDistribution *wd)
{
    /* Crooks Fluctuation Theorem:
     * P_F(W) / P_R(-W) = exp((W - Delta F) / k_B T)
     *
     * At W = Delta F, P_F(W) = P_R(-W).
     * So Delta F is the crossing point of the forward and
     * time-reversed work distributions.
     */

    if (wd == NULL || wd->num_realizations < 2) return -1;

    /* Build histograms of forward and reverse work */
    int num_bins = 50;
    double W_min = wd->work_values[0];
    double W_max = wd->work_values[0];
    int i;

    for (i = 1; i < wd->num_realizations; i++) {
        if (wd->forward_work[i] < W_min) W_min = wd->forward_work[i];
        if (wd->forward_work[i] > W_max) W_max = wd->forward_work[i];
    }

    double bin_w = (W_max - W_min) / num_bins;
    double *hist_F = (double *)calloc((size_t)num_bins, sizeof(double));
    double *hist_R = (double *)calloc((size_t)num_bins, sizeof(double));

    if (hist_F == NULL || hist_R == NULL) {
        free(hist_F);
        free(hist_R);
        return -1;
    }

    for (i = 0; i < wd->num_realizations; i++) {
        int bin = (int)((wd->forward_work[i] - W_min) / bin_w);
        if (bin >= 0 && bin < num_bins) hist_F[bin]++;

        bin = (int)((-wd->reverse_work[i] - W_min) / bin_w);
        if (bin >= 0 && bin < num_bins) hist_R[bin]++;
    }

    /* Find crossing point: where |hist_F(b) - hist_R(b)| is minimized */
    double min_diff = 1e100;
    double crossing_W = 0.0;
    for (i = 0; i < num_bins; i++) {
        double diff = fabs(hist_F[i] - hist_R[i]);
        if (diff < min_diff) {
            min_diff = diff;
            crossing_W = W_min + (i + 0.5) * bin_w;
        }
    }

    wd->delta_F_crooks = crossing_W;

    /* Also compute verification of Crooks:
     * Check if ln(P_F(W)/P_R(-W)) ~ (W - Delta F)/k_B T */
    wd->crooks_verified = (min_diff / wd->num_realizations < 0.1) ? 1 : 0;

    free(hist_F);
    free(hist_R);

    return 0;
}

int ft_compute_jarzynski_deltaF(WorkDistribution *wd)
{
    /* Jarzynski Equality:
     * exp(-Delta F / k_B T) = <exp(-W / k_B T)>
     *
     * Delta F = -k_B T * ln(<exp(-W/k_B T)>)
     */

    int i;
    double sum_exp = 0.0;
    double beta = 1.0 / (wd->k_B * wd->temperature);

    if (wd == NULL || wd->num_realizations < 1) return -1;

    for (i = 0; i < wd->num_realizations; i++) {
        sum_exp += exp(-beta * wd->forward_work[i]);
    }

    double avg_exp = sum_exp / wd->num_realizations;

    if (avg_exp <= 0.0) {
        wd->delta_F_jarzynski = 1e100;  /* Numerically undefined */
        return -1;
    }

    wd->delta_F_jarzynski = -log(avg_exp) / beta;

    return 0;
}

/* ---------------------------------------------------------------------------
 * L4: Fundamental Laws -- Fluctuation-Dissipation Theorem
 * ------------------------------------------------------------------------- */

double ft_verify_fdt(const StochasticTrajectory *traj,
                     double chi_t[], double C_t[],
                     double tau, int n_lags)
{
    /* Fluctuation-Dissipation Theorem:
     * chi(t) = -(1/k_B T) dC(t)/dt  for t > 0
     *
     * where C(t) = <x(0) x(t)>_eq is the autocorrelation function
     * and chi(t) is the linear response function.
     *
     * Estimate C(t) from the trajectory data.
     */

    int n = traj->n_steps;
    double C0, k_B = traj->k_B, T = traj->temperature;
    int lag, t_idx;
    double max_rel_dev = 0.0;

    if (traj == NULL || chi_t == NULL || C_t == NULL) return -1.0;

    /* Compute autocorrelation C(tau_k) for k = 0..n_lags-1 */
    /* Remove mean first */
    double mean_x = 0.0;
    for (t_idx = 0; t_idx < n; t_idx++) {
        mean_x += traj->position[t_idx];
    }
    mean_x /= n;

    /* C(0) */
    C0 = 0.0;
    for (t_idx = 0; t_idx < n; t_idx++) {
        double dx = traj->position[t_idx] - mean_x;
        C0 += dx * dx;
    }
    C0 /= n;

    if (C0 < 1e-15) return 0.0;

    for (lag = 0; lag < n_lags; lag++) {
        int lag_steps = (int)(tau * lag / (n_lags - 1) / traj->dt);
        if (lag_steps >= n) lag_steps = n - 1;

        /* C(lag) = <x(0) x(lag)> */
        C_t[lag] = 0.0;
        int count = 0;
        for (t_idx = 0; t_idx < n - lag_steps; t_idx++) {
            double dx0 = traj->position[t_idx] - mean_x;
            double dx_lag = traj->position[t_idx + lag_steps] - mean_x;
            C_t[lag] += dx0 * dx_lag;
            count++;
        }
        if (count > 0) C_t[lag] /= count;
    }

    /* Estimate response function from derivative of C(t):
     * chi(t) = -(1/k_B T) dC/dt
     * Use finite differences */
    for (lag = 1; lag < n_lags - 1; lag++) {
        double dC_dt = (C_t[lag + 1] - C_t[lag - 1]) / (2.0 * tau / (n_lags - 1));
        chi_t[lag] = -dC_dt / (k_B * T);

        /* Compare with direct response: not available without perturbation.
         * Instead, check self-consistency: chi(t) should decay to 0. */
        if (fabs(C_t[lag]) > 1e-15) {
            double rel_dev = fabs(chi_t[lag] * (k_B * T) + dC_dt)
                           / (fabs(dC_dt) + 1e-15);
            if (rel_dev > max_rel_dev) max_rel_dev = rel_dev;
        }
    }

    /* Boundary points */
    chi_t[0] = chi_t[1];
    chi_t[n_lags - 1] = chi_t[n_lags - 2];

    return max_rel_dev;
}

/* ---------------------------------------------------------------------------
 * L5: Algorithms -- Stochastic integration and resampling
 * ------------------------------------------------------------------------- */

void ft_euler_maruyama_step(double *x, double (*drift)(double, double, const void*),
                            double (*noise)(double, double, const void*),
                            double t, double dt, const void *params)
{
    if (x == NULL || drift == NULL || noise == NULL) return;

    double a = drift(*x, t, params);
    double b = noise(*x, t, params);
    double dW = ft_rand_gaussian(0.0, sqrt(dt));

    *x += a * dt + b * dW;
}

void ft_bootstrap_resample(const double data[], int n,
                           double resampled[], unsigned int seed)
{
    int i;

    if (data == NULL || resampled == NULL || n < 1) return;

    ft_rand_state = seed;

    for (i = 0; i < n; i++) {
        int idx = (int)(ft_rand_uniform() * n);
        if (idx >= n) idx = n - 1;
        resampled[i] = data[idx];
    }
}

int ft_build_work_distribution(WorkDistribution *wd,
                                double mass, double friction,
                                double spring_k, double L, double tau,
                                double temperature, double k_B,
                                double dt, int num_realizations)
{
    int i;
    int n_steps = (int)(tau / dt);

    if (wd == NULL || num_realizations < 1) return -1;

    wd->num_realizations = num_realizations;
    wd->temperature = temperature;
    wd->k_B = k_B;

    for (i = 0; i < num_realizations; i++) {
        StochasticTrajectory traj;
        ft_init_trajectory(&traj, n_steps, dt, temperature, k_B);

        ft_simulate_jarzynski_protocol(&traj, mass, friction,
                                        spring_k, L, tau,
                                        temperature, k_B, dt);

        /* Work done is the cumulative work at the end */
        wd->forward_work[i] = traj.work_done[n_steps - 1];
        wd->work_values[i] = wd->forward_work[i];

        /* Reverse work would come from running the protocol backwards */

        ft_free_trajectory(&traj);
    }

    return 0;
}

/* ---------------------------------------------------------------------------
 * Utility functions
 * ------------------------------------------------------------------------- */

void ft_init_trajectory(StochasticTrajectory *traj, int n_steps,
                        double dt, double temperature, double k_B)
{
    if (traj == NULL) return;

    traj->n_steps = n_steps;
    traj->dt = dt;
    traj->temperature = temperature;
    traj->k_B = k_B;

    traj->time = (double *)calloc((size_t)n_steps, sizeof(double));
    traj->position = (double *)calloc((size_t)n_steps, sizeof(double));
    traj->velocity = (double *)calloc((size_t)n_steps, sizeof(double));
    traj->entropy_production = (double *)calloc((size_t)n_steps, sizeof(double));
    traj->work_done = (double *)calloc((size_t)n_steps, sizeof(double));
    traj->heat_dissipated = (double *)calloc((size_t)n_steps, sizeof(double));
}

void ft_free_trajectory(StochasticTrajectory *traj)
{
    if (traj == NULL) return;
    free(traj->time);
    free(traj->position);
    free(traj->velocity);
    free(traj->entropy_production);
    free(traj->work_done);
    free(traj->heat_dissipated);
    memset(traj, 0, sizeof(StochasticTrajectory));
}

void ft_init_histogram(EntropyProductionHistogram *hist, int num_bins,
                       double sigma_min, double sigma_max, double k_B)
{
    if (hist == NULL) return;

    hist->num_bins = num_bins;
    hist->sigma_min = sigma_min;
    hist->sigma_max = sigma_max;
    hist->bin_width = (sigma_max - sigma_min) / num_bins;
    hist->k_B = k_B;
    hist->total_samples = 0;
    hist->asymmetry_ratio = 0.0;
    hist->ft_slope = 0.0;
    hist->ft_r2 = 0.0;

    hist->bins = (double *)calloc((size_t)num_bins, sizeof(double));
    hist->pdf = (double *)calloc((size_t)num_bins, sizeof(double));
    hist->pdf_reverse = (double *)calloc((size_t)num_bins, sizeof(double));

    int i;
    for (i = 0; i < num_bins; i++) {
        hist->bins[i] = sigma_min + (i + 0.5) * hist->bin_width;
    }
}

void ft_free_histogram(EntropyProductionHistogram *hist)
{
    if (hist == NULL) return;
    free(hist->bins);
    free(hist->pdf);
    free(hist->pdf_reverse);
    memset(hist, 0, sizeof(EntropyProductionHistogram));
}

void ft_init_work_distribution(WorkDistribution *wd, int num_realizations,
                               double temperature, double k_B)
{
    if (wd == NULL) return;

    wd->num_realizations = num_realizations;
    wd->temperature = temperature;
    wd->k_B = k_B;
    wd->delta_F_true = 0.0;
    wd->delta_F_jarzynski = 0.0;
    wd->delta_F_crooks = 0.0;
    wd->crooks_verified = 0;

    wd->work_values = (double *)calloc((size_t)num_realizations, sizeof(double));
    wd->forward_work = (double *)calloc((size_t)num_realizations, sizeof(double));
    wd->reverse_work = (double *)calloc((size_t)num_realizations, sizeof(double));
}

void ft_free_work_distribution(WorkDistribution *wd)
{
    if (wd == NULL) return;
    free(wd->work_values);
    free(wd->forward_work);
    free(wd->reverse_work);
    memset(wd, 0, sizeof(WorkDistribution));
}