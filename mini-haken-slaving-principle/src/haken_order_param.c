/* ==============================================================
 * haken_order_param.c — Order Parameter Identification & GL Analysis
 *
 * The order parameter concept is the cornerstone of Haken's
 * synergetics. This module identifies order parameters from
 * spectral decomposition, principal component analysis, mean-field
 * averaging, and Fourier mode amplitudes. It also fits and analyzes
 * the Ginzburg-Landau equation that governs order parameter dynamics.
 *
 * Ginzburg-Landau theory (1950):
 *   F(ξ) = F₀ + (a/2)ξ² + (b/4)ξ⁴   (b > 0 for stability)
 *   dξ/dt = -∂F/∂ξ = -a·ξ - b·ξ³
 *
 * Critical exponents (mean-field):
 *   β = 1/2    ξ_eq ~ |α - α_c|^β
 *   γ = 1      χ ~ |α - α_c|^{-γ}
 *   ν = 1/2    correlation length
 *
 * References:
 *   Haken (1977) Synergetics, §2.2-2.4, §6
 *   Ma (1976) Modern Theory of Critical Phenomena
 *   Cross & Greenside (2009) Pattern Formation, §5
 *   Landau & Lifshitz (1980) Statistical Physics, Part 1
 * ============================================================== */

#include "haken_order_param.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==============================================================
 * L5: Spectral order parameter identification
 * ============================================================== */

void haken_identify_spectral(const Haken_System* sys, int n_ops,
                              double eps, Haken_OrderParam* result) {
    (void)eps;
    if (!sys || !result) return;

    int n = sys->n_dim;
    /* Find modes with smallest |Re(λ)| */
    typedef struct { int idx; double abs_re; } Candidate;
    Candidate* cand = (Candidate*)malloc((size_t)n * sizeof(Candidate));
    for (int k = 0; k < n; k++) {
        cand[k].idx = k;
        cand[k].abs_re = fabs(sys->modes[k].lambda_re);
    }
    /* Sort by absolute real part ascending */
    for (int i = 0; i < n - 1; i++)
        for (int j = i + 1; j < n; j++)
            if (cand[j].abs_re < cand[i].abs_re) {
                Candidate tmp = cand[i]; cand[i] = cand[j]; cand[j] = tmp;
            }

    int n_keep = n_ops < n ? n_ops : n;
    result->n_op = n_keep;
    result->amplitude = (double*)calloc((size_t)n_keep, sizeof(double));
    result->effective_lambda = (double*)calloc((size_t)n_keep, sizeof(double));
    result->projection_matrix = (double*)calloc((size_t)(n_keep * n), sizeof(double));

    for (int a = 0; a < n_keep; a++) {
        int k = cand[a].idx;
        result->effective_lambda[a] = sys->modes[k].lambda_re;
        result->amplitude[a] = 0.0;
        /* Projection: row a of projection matrix is right eigenvector of mode k */
        for (int i = 0; i < n; i++)
            result->projection_matrix[a * n + i] = sys->modes[k].eigenvector_re[i];
    }

    /* Set GL parameters */
    result->gl_potential_coeff = 1.0;
    result->gl_cubic_coeff     = 1.0;
    result->gl_correlation_len  = 1.0 / sqrt(fabs(sys->control_param - sys->critical_param) + 1e-6);
    result->gl_relaxation_time  = 1.0 / (fabs(sys->control_param - sys->critical_param) + 1e-6);

    free(cand);
}

/* ==============================================================
 * PCA-based order parameter identification
 * ============================================================== */

void haken_identify_pca(const double* trajectory, int n_steps,
                         int n_dim, int n_ops, Haken_OrderParam* result) {
    if (!trajectory || n_steps < 2 || n_dim <= 0 || !result) return;

    /* Compute mean */
    double* mean = (double*)calloc((size_t)n_dim, sizeof(double));
    for (int t = 0; t < n_steps; t++)
        for (int i = 0; i < n_dim; i++)
            mean[i] += trajectory[t * n_dim + i];
    for (int i = 0; i < n_dim; i++) mean[i] /= n_steps;

    /* Compute covariance matrix C_{ij} = <(q_i - μ_i)(q_j - μ_j)> */
    double* cov = (double*)calloc((size_t)(n_dim * n_dim), sizeof(double));
    for (int t = 0; t < n_steps; t++) {
        for (int i = 0; i < n_dim; i++) {
            double di = trajectory[t * n_dim + i] - mean[i];
            for (int j = 0; j < n_dim; j++) {
                double dj = trajectory[t * n_dim + j] - mean[j];
                cov[i * n_dim + j] += di * dj;
            }
        }
    }
    for (int i = 0; i < n_dim * n_dim; i++) cov[i] /= (n_steps - 1);

    /* Power iteration for top n_ops eigenvectors */
    int n_keep = n_ops < n_dim ? n_ops : n_dim;
    result->n_op = n_keep;
    result->projection_matrix = (double*)calloc((size_t)(n_keep * n_dim), sizeof(double));
    result->effective_lambda = (double*)calloc((size_t)n_keep, sizeof(double));

    for (int a = 0; a < n_keep; a++) {
        double* v = (double*)calloc((size_t)n_dim, sizeof(double));
        for (int i = 0; i < n_dim; i++) v[i] = (double)(rand()) / RAND_MAX;

        /* Power iteration with deflation */
        for (int iter = 0; iter < 50; iter++) {
            double* v_new = (double*)calloc((size_t)n_dim, sizeof(double));
            for (int i = 0; i < n_dim; i++)
                for (int j = 0; j < n_dim; j++)
                    v_new[i] += cov[i * n_dim + j] * v[j];
            /* Deflate against previous eigenvectors */
            for (int b = 0; b < a; b++)
                for (int i = 0; i < n_dim; i++)
                    v_new[i] -= result->projection_matrix[b * n_dim + i] * v[i];
            /* Normalize */
            double norm = 0.0;
            for (int i = 0; i < n_dim; i++) norm += v_new[i] * v_new[i];
            norm = sqrt(norm);
            if (norm < 1e-15) { free(v_new); break; }
            for (int i = 0; i < n_dim; i++) v[i] = v_new[i] / norm;
            free(v_new);
        }
        /* Rayleigh quotient for eigenvalue */
        double rayleigh = 0.0;
        for (int i = 0; i < n_dim; i++)
            for (int j = 0; j < n_dim; j++)
                rayleigh += v[i] * cov[i * n_dim + j] * v[j];
        result->effective_lambda[a] = rayleigh;
        for (int i = 0; i < n_dim; i++) result->projection_matrix[a * n_dim + i] = v[i];
        free(v);
    }

    result->amplitude = (double*)calloc((size_t)n_keep, sizeof(double));
    free(mean); free(cov);
}

/* ==============================================================
 * Mean-field order parameter
 * ============================================================== */

double haken_mean_field_op(const double* q, int n) {
    if (!q || n <= 0) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < n; i++) sum += q[i];
    return sum / n;
}

/* ==============================================================
 * Fourier mode amplitude
 * ============================================================== */

double haken_fourier_mode_amplitude(const double* q, int n, int k) {
    if (!q || n <= 0 || k < 0 || k >= n) return 0.0;
    double re = 0.0, im = 0.0;
    for (int j = 0; j < n; j++) {
        double theta = -2.0 * M_PI * j * k / n;
        re += q[j] * cos(theta);
        im += q[j] * sin(theta);
    }
    return (re * re + im * im) / (n * n);
}

void haken_spatial_correlation(const double* field, int n, int r,
                                double* corr) {
    if (!field || n <= 0 || !corr) return;
    double mean = 0.0;
    for (int i = 0; i < n; i++) mean += field[i];
    mean /= n;

    double var = 0.0;
    for (int i = 0; i < n; i++) { double d = field[i] - mean; var += d * d; }
    if (var < 1e-15) { for (int lag = 0; lag <= r; lag++) corr[lag] = 0.0; return; }

    for (int lag = 0; lag <= r && lag < n; lag++) {
        double sum = 0.0;
        int count = 0;
        for (int i = 0; i < n - lag; i++) {
            sum += (field[i] - mean) * (field[i + lag] - mean);
            count++;
        }
        corr[lag] = (count > 0) ? (sum / count) / (var / n) : 0.0;
    }
}

double haken_correlation_length(const double* correlation, int r_max) {
    if (!correlation || r_max < 1) return 0.0;
    /* Fit log(C(r)) ≈ -r/ξ → linear regression of log(C) vs r */
    double sum_x = 0, sum_y = 0, sum_xx = 0, sum_xy = 0;
    int n = 0;
    for (int r = 1; r <= r_max; r++) {
        if (correlation[r] <= 0) continue;
        double x = (double)r;
        double y = log(correlation[r]);
        sum_x += x; sum_y += y; sum_xx += x * x; sum_xy += x * y;
        n++;
    }
    if (n < 2) return 0.0;
    double denom = n * sum_xx - sum_x * sum_x;
    if (fabs(denom) < 1e-15) return 0.0;
    double slope = (n * sum_xy - sum_x * sum_y) / denom;
    return (slope < 0) ? (-1.0 / slope) : INFINITY;
}

/* ==============================================================
 * Ginzburg-Landau fitting
 * ============================================================== */

void haken_fit_gl_parameters(const Haken_OPTimeSeries* ts, double dt,
                              Haken_GLParameters* params) {
    if (!ts || !params || ts->n_samples < 3) return;

    /* Fit dξ/dt = -a·ξ - b·ξ³ via least squares */
    int n = ts->n_samples - 1;
    double sum_y2 = 0, sum_yy3 = 0, sum_y3y3 = 0, sum_rhs_y = 0, sum_rhs_y3 = 0;

    for (int t = 0; t < n; t++) {
        double xi = ts->data[t];
        double xi3 = xi * xi * xi;
        double dxi_dt = (ts->data[t + 1] - ts->data[t]) / dt;
        sum_y2 += xi * xi;
        sum_yy3 += xi * xi3;
        sum_y3y3 += xi3 * xi3;
        sum_rhs_y += dxi_dt * xi;
        sum_rhs_y3 += dxi_dt * xi3;
    }

    /* Solve [sum_y2, sum_yy3; sum_yy3, sum_y3y3] * [a; b] = [-sum_rhs_y; -sum_rhs_y3] */
    double det = sum_y2 * sum_y3y3 - sum_yy3 * sum_yy3;
    if (fabs(det) < 1e-15) {
        params->a = 0.0; params->b = 1.0;
    } else {
        params->a = (-sum_rhs_y * sum_y3y3 + sum_rhs_y3 * sum_yy3) / det;
        params->b = (sum_y2 * (-sum_rhs_y3) + sum_yy3 * sum_rhs_y) / det;
    }
    params->critical_alpha = 0.0;
    params->xi0 = 1.0;
    params->tau0 = 1.0;
}

double haken_gl_free_energy(double xi, const Haken_GLParameters* params) {
    if (!params) return 0.0;
    return -(params->a / 2.0) * xi * xi - (params->b / 4.0) * xi * xi * xi * xi;
}

void haken_gl_potential_landscape(const Haken_GLParameters* params,
                                   double xi_min, double xi_max,
                                   int n_points,
                                   double* xi_grid, double* F_grid) {
    if (!params || !xi_grid || !F_grid || n_points <= 0) return;
    double step = (xi_max - xi_min) / (n_points - 1);
    for (int i = 0; i < n_points; i++) {
        double xi = xi_min + step * i;
        xi_grid[i] = xi;
        F_grid[i] = haken_gl_free_energy(xi, params);
    }
}

int haken_gl_fixed_points(const Haken_GLParameters* params,
                           double* fixed_xi) {
    if (!params || !fixed_xi) return 0;
    double a = params->a;
    double b = params->b;
    if (fabs(b) < 1e-15) { fixed_xi[0] = 0.0; return 1; }
    fixed_xi[0] = 0.0;
    if (a * b < 0) {
        double val = sqrt(-a / b);
        fixed_xi[1] = val;
        fixed_xi[2] = -val;
        return 3;
    }
    return 1;
}

void haken_gl_stability(const Haken_GLParameters* params,
                         const double* fixed_xi, int n_fixed,
                         double* stability) {
    if (!params || !fixed_xi || !stability) return;
    for (int k = 0; k < n_fixed; k++) {
        double xi = fixed_xi[k];
        /* λ = -a - 3b·ξ² */
        stability[k] = -params->a - 3.0 * params->b * xi * xi;
    }
}

/* ==============================================================
 * Critical exponents estimation
 * ============================================================== */

void haken_critical_exponents(const Haken_ParameterSweep* sweep,
                               Haken_GLParameters* params) {
    if (!sweep || !params || sweep->n_points < 4) return;

    /* Estimate β from order parameter amplitude vs |α - α_c| */
    double sum_log_eps = 0, sum_log_amp = 0, sum_log_eps2 = 0, sum_log_prod = 0;
    int count = 0;
    for (int i = 0; i < sweep->n_points; i++) {
        double alpha = sweep->alpha_start + i * sweep->alpha_step;
        double eps = fabs(alpha - params->critical_alpha);
        if (eps < 1e-10) continue;
        double amp = 0.0;
        if (sweep->order_param_amp) amp = sweep->order_param_amp[i];
        if (amp < 1e-10) continue;
        double log_eps = log(eps), log_amp = log(amp);
        sum_log_eps += log_eps; sum_log_amp += log_amp;
        sum_log_eps2 += log_eps * log_eps; sum_log_prod += log_eps * log_amp;
        count++;
    }
    if (count >= 3) {
        double denom = count * sum_log_eps2 - sum_log_eps * sum_log_eps;
        if (fabs(denom) > 1e-15) {
            params->beta = (count * sum_log_prod - sum_log_eps * sum_log_amp) / denom;
        }
    }
}

/* ==============================================================
 * Order parameter validation
 * ============================================================== */

bool haken_validate_order_parameter(const Haken_System* sys,
                                     int candidate_mode,
                                     double tolerance) {
    if (!sys || candidate_mode < 0 || candidate_mode >= sys->n_dim) return false;

    double lambda_re = sys->modes[candidate_mode].lambda_re;
    double alpha = sys->control_param;
    double alpha_c = sys->critical_param;

    /* Condition 1: Re(λ) ≈ 0 near critical point */
    if (fabs(alpha - alpha_c) < tolerance) {
        if (fabs(lambda_re) > tolerance) return false;
    }

    /* Condition 2: Re(λ) < 0 below threshold (stable disordered state) */
    if (alpha < alpha_c - tolerance) {
        if (lambda_re > -tolerance) return false;
    }

    /* Condition 3: Spectral gap exists */
    double gap = haken_spectral_gap(sys);
    if (gap < tolerance) return false;

    return true;
}

double haken_variance_explained(const Haken_System* sys,
                                 const Haken_OrderParam* op) {
    if (!sys || !op || op->n_op <= 0) return 0.0;

    /* Fraction of variance captured = Σ_k (eigenvalue of OP mode)² / Σ_j λ_j² */
    double total_var = 0.0, op_var = 0.0;
    for (int k = 0; k < sys->n_dim; k++) {
        double lam = sys->modes[k].lambda_re;
        total_var += lam * lam;
    }
    for (int a = 0; a < op->n_op; a++) {
        op_var += op->effective_lambda[a] * op->effective_lambda[a];
    }
    return (total_var > 1e-15) ? (op_var / total_var) : 0.0;
}

/* ==============================================================
 * Time series analysis
 * ============================================================== */

/* Compute power spectrum via simple periodogram */
void haken_power_spectrum(double* data, int n, double* ps) {
    for (int k = 0; k < n / 2 + 1; k++) {
        double re = 0.0, im = 0.0;
        for (int t = 0; t < n; t++) {
            double theta = -2.0 * M_PI * k * t / n;
            re += data[t] * cos(theta);
            im += data[t] * sin(theta);
        }
        ps[k] = (re * re + im * im) / n;
    }
}
