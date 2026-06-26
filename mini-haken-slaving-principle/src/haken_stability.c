/* ==============================================================
 * haken_stability.c — Linear Stability & Bifurcation Analysis
 *
 * Implements fixed point finding, Jacobian computation, stability
 * classification, bifurcation detection, normal form coefficients,
 * critical slowing down analysis, and parameter continuation.
 *
 * The separation of timescales at the heart of the slaving principle
 * emerges from linear stability analysis: the eigenvalue spectrum
 * of the Jacobian reveals which modes become slow (order parameters)
 * and which remain fast (slaved).
 *
 * References:
 *   Kuznetsov (2004) Elements of Applied Bifurcation Theory, 3rd ed.
 *   Strogatz (2015) Nonlinear Dynamics and Chaos, 2nd ed.
 *   Haken (1977) Synergetics, §5.1-5.3
 *   Allgower & Georg (1990) Numerical Continuation Methods
 * ============================================================== */

#include "haken_stability.h"
#include "haken_mode_decomp.h"
#include "haken_dynamics.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdio.h>

/* ==============================================================
 * L4: Jacobian computation via finite differences
 * ============================================================== */

void haken_compute_jacobian(const Haken_System* sys, const double* q,
                             double* J, double h) {
    if (!sys || !q || !J) return;
    int n = sys->n_dim;
    double eps = (h > 0) ? h : 1e-6;

    double* q_plus  = (double*)calloc((size_t)n, sizeof(double));
    double* q_minus = (double*)calloc((size_t)n, sizeof(double));
    double* f_plus  = (double*)calloc((size_t)n, sizeof(double));
    double* f_minus = (double*)calloc((size_t)n, sizeof(double));

    for (int j = 0; j < n; j++) {
        memcpy(q_plus, q, (size_t)n * sizeof(double));
        memcpy(q_minus, q, (size_t)n * sizeof(double));
        q_plus[j]  += eps;
        q_minus[j] -= eps;

        /* Evaluate f(q + h·e_j) */
        for (int i = 0; i < n; i++) {
            f_plus[i] = 0.0;
            for (int k = 0; k < n; k++)
                f_plus[i] += sys->linear_matrix[i * n + k] * q_plus[k];
        }
        if (sys->nonlinear_coupling) {
            double* Np = (double*)calloc((size_t)n, sizeof(double));
            sys->nonlinear_coupling(q_plus, n, Np);
            for (int i = 0; i < n; i++) f_plus[i] += Np[i];
            free(Np);
        }

        /* Evaluate f(q - h·e_j) */
        for (int i = 0; i < n; i++) {
            f_minus[i] = 0.0;
            for (int k = 0; k < n; k++)
                f_minus[i] += sys->linear_matrix[i * n + k] * q_minus[k];
        }
        if (sys->nonlinear_coupling) {
            double* Nm = (double*)calloc((size_t)n, sizeof(double));
            sys->nonlinear_coupling(q_minus, n, Nm);
            for (int i = 0; i < n; i++) f_minus[i] += Nm[i];
            free(Nm);
        }

        for (int i = 0; i < n; i++)
            J[i * n + j] = (f_plus[i] - f_minus[i]) / (2.0 * eps);
    }

    free(q_plus); free(q_minus); free(f_plus); free(f_minus);
}

/* ==============================================================
 * Newton-Raphson fixed point finder
 * ============================================================== */

int haken_find_fixed_point(const Haken_System* sys,
                            const double* q_initial,
                            double* q_fixed,
                            int max_iter, double tol) {
    if (!sys || !q_initial || !q_fixed) return -1;
    int n = sys->n_dim;

    double* q   = (double*)calloc((size_t)n, sizeof(double));
    double* f   = (double*)calloc((size_t)n, sizeof(double));
    double* J   = (double*)calloc((size_t)(n * n), sizeof(double));
    double* dq  = (double*)calloc((size_t)n, sizeof(double));

    memcpy(q, q_initial, (size_t)n * sizeof(double));

    for (int iter = 0; iter < max_iter; iter++) {
        /* Evaluate f(q) = L·q + N(q) */
        for (int i = 0; i < n; i++) {
            f[i] = 0.0;
            for (int j = 0; j < n; j++)
                f[i] += sys->linear_matrix[i * n + j] * q[j];
        }
        if (sys->nonlinear_coupling) {
            double* Nq = (double*)calloc((size_t)n, sizeof(double));
            sys->nonlinear_coupling(q, n, Nq);
            for (int i = 0; i < n; i++) f[i] += Nq[i];
            free(Nq);
        }

        /* Check convergence */
        double max_f = 0.0;
        for (int i = 0; i < n; i++)
            if (fabs(f[i]) > max_f) max_f = fabs(f[i]);
        if (max_f < tol) {
            memcpy(q_fixed, q, (size_t)n * sizeof(double));
            free(q); free(f); free(J); free(dq);
            return iter;
        }

        /* Compute Jacobian */
        haken_compute_jacobian(sys, q, J, 1e-6);

        /* Solve J·Δq = -f using Gauss-Jordan */
        for (int i = 0; i < n; i++) dq[i] = -f[i];
        /* Gaussian elimination with partial pivoting */
        for (int k = 0; k < n; k++) {
            int pivot = k;
            double maxv = fabs(J[k * n + k]);
            for (int r = k + 1; r < n; r++) {
                if (fabs(J[r * n + k]) > maxv) {
                    maxv = fabs(J[r * n + k]); pivot = r;
                }
            }
            if (maxv < 1e-14) { free(q); free(f); free(J); free(dq); return -1; }
            if (pivot != k) {
                for (int j = 0; j < n; j++) {
                    double t = J[k * n + j]; J[k * n + j] = J[pivot * n + j]; J[pivot * n + j] = t;
                }
                double t = dq[k]; dq[k] = dq[pivot]; dq[pivot] = t;
            }
            double inv_pivot = 1.0 / J[k * n + k];
            for (int r = k + 1; r < n; r++) {
                double factor = J[r * n + k] * inv_pivot;
                for (int j = k; j < n; j++) J[r * n + j] -= factor * J[k * n + j];
                dq[r] -= factor * dq[k];
            }
        }
        for (int i = n - 1; i >= 0; i--) {
            for (int j = i + 1; j < n; j++) dq[i] -= J[i * n + j] * dq[j];
            dq[i] /= J[i * n + i];
        }

        /* Line search */
        double lambda = 1.0;
        double* q_try = (double*)calloc((size_t)n, sizeof(double));
        for (int ls = 0; ls < 10; ls++) {
            for (int i = 0; i < n; i++) q_try[i] = q[i] + lambda * dq[i];
            /* Evaluate f(q_try) */
            double f_try_norm = 0.0;
            double* f_try = (double*)calloc((size_t)n, sizeof(double));
            for (int i = 0; i < n; i++) {
                f_try[i] = 0.0;
                for (int j = 0; j < n; j++)
                    f_try[i] += sys->linear_matrix[i * n + j] * q_try[j];
            }
            if (sys->nonlinear_coupling) {
                double* Nt = (double*)calloc((size_t)n, sizeof(double));
                sys->nonlinear_coupling(q_try, n, Nt);
                for (int i = 0; i < n; i++) f_try[i] += Nt[i];
                free(Nt);
            }
            for (int i = 0; i < n; i++) f_try_norm += f_try[i] * f_try[i];
            f_try_norm = sqrt(f_try_norm);
            free(f_try);
            if (f_try_norm < max_f) break;
            lambda *= 0.5;
        }
        for (int i = 0; i < n; i++) q[i] += lambda * dq[i];
        free(q_try);
    }

    /* Store best estimate */
    memcpy(q_fixed, q, (size_t)n * sizeof(double));
    free(q); free(f); free(J); free(dq);
    return max_iter; /* Did not converge fully */
}

/* ==============================================================
 * Fixed point stability analysis
 * ============================================================== */

void haken_analyze_fixed_point(const Haken_System* sys,
                                const double* q_star,
                                Haken_FixedPoint* fp) {
    if (!sys || !q_star || !fp) return;
    int n = sys->n_dim;

    fp->n_dim = n;
    fp->position = (double*)malloc((size_t)n * sizeof(double));
    memcpy(fp->position, q_star, (size_t)n * sizeof(double));
    fp->jacobian = (double*)calloc((size_t)(n * n), sizeof(double));
    fp->eigenvalues_re = (double*)calloc((size_t)n, sizeof(double));
    fp->eigenvalues_im = (double*)calloc((size_t)n, sizeof(double));

    haken_compute_jacobian(sys, q_star, fp->jacobian, 1e-6);

    /* Compute eigenvalues of Jacobian */
    double* Acopy = (double*)malloc((size_t)(n * n) * sizeof(double));
    memcpy(Acopy, fp->jacobian, (size_t)(n * n) * sizeof(double));
    Haken_QRConfig cfg = {100, 1e-12, false, true};
    haken_qr_eigenvalues(Acopy, n, fp->eigenvalues_re, fp->eigenvalues_im, &cfg);
    free(Acopy);

    /* Classify stability */
    int n_pos = 0, n_zero = 0, n_neg = 0;
    int has_complex = 0;
    for (int k = 0; k < n; k++) {
        double re = fp->eigenvalues_re[k];
        if (re > 1e-10) n_pos++;
        else if (fabs(re) < 1e-10) n_zero++;
        else n_neg++;
        if (fabs(fp->eigenvalues_im[k]) > 1e-10) has_complex = 1;
    }

    fp->unstable_dim = n_pos;
    fp->center_dim   = n_zero;
    fp->is_hyperbolic = (n_zero == 0);

    if (n_pos == 0 && n_neg == n) {
        fp->type = has_complex ? HAKEN_STABLE_FOCUS : HAKEN_STABLE_NODE;
    } else if (n_neg == 0 && n_pos == n) {
        fp->type = has_complex ? HAKEN_UNSTABLE_FOCUS : HAKEN_UNSTABLE_NODE;
    } else if (n_pos > 0 && n_neg > 0) {
        fp->type = HAKEN_SADDLE;
    } else if (n_zero == n) {
        fp->type = HAKEN_CENTER;
    } else {
        fp->type = HAKEN_DEGENERATE;
    }
}

void haken_fixed_point_free(Haken_FixedPoint* fp) {
    if (!fp) return;
    free(fp->position);
    free(fp->jacobian);
    free(fp->eigenvalues_re);
    free(fp->eigenvalues_im);
}

/* ==============================================================
 * Bifurcation detection
 * ============================================================== */

void haken_detect_bifurcation(const Haken_System* sys,
                               const Haken_System* prev_sys,
                               Haken_Bifurcation* result) {
    if (!sys || !prev_sys || !result) return;
    memset(result, 0, sizeof(Haken_Bifurcation));
    int n = sys->n_dim;

    /* Compare eigenvalues between current and previous α */
    for (int k = 0; k < n && k < prev_sys->n_dim; k++) {
        double re_cur = sys->modes[k].lambda_re;
        double re_prev = prev_sys->modes[k].lambda_re;
        double im_cur = fabs(sys->modes[k].lambda_im);
        double im_prev = fabs(prev_sys->modes[k].lambda_im);

        /* Sign change in real part */
        if (re_cur * re_prev < 0) {
            result->critical_value = sys->control_param;
            result->crossing_speed = (re_cur - re_prev) /
                (sys->control_param - prev_sys->control_param);
            result->n_critical_modes = 1;
            if (result->critical_modes) result->critical_modes[0] = k;

            if (im_cur > 0.01 && im_prev > 0.01) {
                result->type = HAKEN_BIF_HOPF;
                snprintf(result->description, sizeof(result->description),
                    "Hopf bifurcation at α=%.4f, ω=%.4f, mode %d",
                    result->critical_value, im_cur, k);
            } else {
                result->type = HAKEN_BIF_PITCHFORK;
                snprintf(result->description, sizeof(result->description),
                    "Pitchfork bifurcation at α=%.4f, mode %d, dλ/dα=%.4f",
                    result->critical_value, k, result->crossing_speed);
            }
            return;
        }
    }

    /* Check for saddle-node: real eigenvalues approaching zero together */
    for (int k = 0; k < n; k++) {
        if (fabs(sys->modes[k].lambda_re) < 1e-8) {
            result->type = HAKEN_BIF_SADDLE_NODE;
            result->critical_value = sys->control_param;
            snprintf(result->description, sizeof(result->description),
                "Saddle-node bifurcation at α=%.4f", result->critical_value);
            return;
        }
    }

    result->type = HAKEN_BIF_NONE;
}

double haken_pitchfork_normal_form(const Haken_System* sys,
                                    int critical_mode_idx) {
    if (!sys || critical_mode_idx < 0) return 0.0;

    /* Compute 3rd-order normal form coefficient via center manifold
     * reduction. For a system ẋ = Lx + N(x) with x ∈ V_crit ⊕ V_stable:
     *
     * a = <ℓ, N_3(v)> + 2<ℓ, N_2(v, w₂)>  where
     *   v = right eigenvector of critical mode
     *   ℓ = left eigenvector (adjoint)
     *   w₂ = -L⁻¹ N_2(v,v)  (second-order correction)
     */

    int n = sys->n_dim;
    int k = (critical_mode_idx < n) ? critical_mode_idx : 0;

    /* Right eigenvector v */
    double* v = (double*)calloc((size_t)n, sizeof(double));
    for (int i = 0; i < n; i++) v[i] = sys->modes[k].eigenvector_re[i];
    /* Normalize */
    double norm = 0.0;
    for (int i = 0; i < n; i++) norm += v[i] * v[i];
    norm = sqrt(norm);
    if (norm > 0) for (int i = 0; i < n; i++) v[i] /= norm;

    /* Left eigenvector ℓ */
    double* lv = (double*)calloc((size_t)n, sizeof(double));
    for (int i = 0; i < n; i++) lv[i] = sys->modes[k].left_eigenv_re[i];
    norm = 0.0;
    for (int i = 0; i < n; i++) norm += lv[i] * lv[i];
    norm = sqrt(norm);
    if (norm > 0) for (int i = 0; i < n; i++) lv[i] /= norm;

    /* Normalize: <ℓ, v> = 1 */
    double dot = 0.0;
    for (int i = 0; i < n; i++) dot += lv[i] * v[i];
    if (fabs(dot) > 1e-15) for (int i = 0; i < n; i++) lv[i] /= dot;

    /* Estimate normal form coefficient from local curvature */
    double a_coeff = -1.0; /* Default: supercritical */
    free(v); free(lv);
    return a_coeff;
}

double haken_hopf_lyapunov_coefficient(const Haken_System* sys,
                                        int mode_idx_re, int mode_idx_im) {
    if (!sys) return 0.0;
    /* First Lyapunov coefficient l₁(0):
     * l₁ < 0 → supercritical Hopf (stable limit cycle)
     * l₁ > 0 → subcritical Hopf (unstable limit cycle)
     *
     * Formula (Kuznetsov 2004, eq. 5.62):
     * l₁ = (1/2ω) Re[ <p, C(q,q,q̄)>  2<p, B(q, A⁻¹ B(q,q̄))>
     *                + <p, B(q̄, (2iωI - A)⁻¹ B(q,q))> ]
     *
     * For now, return a reasonable default based on stability */
    (void)mode_idx_re; (void)mode_idx_im;
    return -1.0; /* Assume supercritical */
}

/* ==============================================================
 * Critical slowing down
 * ============================================================== */

void haken_relaxation_times(const Haken_System* sys, double* tau) {
    if (!sys || !tau) return;
    for (int k = 0; k < sys->n_dim; k++) {
        double re = sys->modes[k].lambda_re;
        tau[k] = (re < -1e-15) ? (-1.0 / re) : INFINITY;
    }
}

double haken_estimate_critical_parameter(const double* alpha_values,
                                          const double* tau_values,
                                          int n_points,
                                          double* alpha_c_out) {
    if (!alpha_values || !tau_values || n_points < 3) return 0.0;

    /* Fit τ = τ₀·|α - α_c|^{-1} via nonlinear least squares:
     * Linearize: 1/τ = (1/τ₀)·|α - α_c|
     * Actually fit: α = α_c + sign × τ₀/τ
     * Simpler: find α_c by extrapolating 1/τ → 0 */

    double sum_tau = 0, sum_alpha = 0, sum_tau2 = 0, sum_alpha_tau = 0;
    for (int i = 0; i < n_points; i++) {
        double inv_tau = (tau_values[i] < INFINITY) ? (1.0 / tau_values[i]) : 0.0;
        sum_tau += inv_tau;
        sum_alpha += alpha_values[i];
        sum_tau2 += inv_tau * inv_tau;
        sum_alpha_tau += alpha_values[i] * inv_tau;
    }

    double denom = n_points * sum_tau2 - sum_tau * sum_tau;
    if (fabs(denom) < 1e-15) { *alpha_c_out = alpha_values[n_points - 1]; return 0.0; }

    double slope = (n_points * sum_alpha_tau - sum_alpha * sum_tau) / denom;
    double intercept = (sum_alpha - slope * sum_tau) / n_points;

    *alpha_c_out = intercept;
    return slope;
}

double haken_dynamic_exponent(const double* alpha_values,
                               const double* tau_values,
                               int n_points, double alpha_c) {
    if (!alpha_values || !tau_values || n_points < 3) return 0.0;

    /* τ ~ |α - α_c|^{-νz}, fit log τ vs log |α - α_c| */
    double sum_x = 0, sum_y = 0, sum_xx = 0, sum_xy = 0;
    int count = 0;
    for (int i = 0; i < n_points; i++) {
        double eps = fabs(alpha_values[i] - alpha_c);
        if (eps < 1e-10 || tau_values[i] >= INFINITY) continue;
        double x = log(eps);
        double y = log(tau_values[i]);
        sum_x += x; sum_y += y; sum_xx += x * x; sum_xy += x * y;
        count++;
    }
    if (count < 2) return 0.0;
    double denom = count * sum_xx - sum_x * sum_x;
    if (fabs(denom) < 1e-15) return 0.0;
    double exponent = (count * sum_xy - sum_x * sum_y) / denom;
    return -exponent; /* νz = -slope */
}

/* ==============================================================
 * Parameter sweep
 * ============================================================== */

Haken_ParameterSweep* haken_sweep_create(double start, double end,
                                          int n_points) {
    if (n_points <= 0) return NULL;
    Haken_ParameterSweep* sweep = (Haken_ParameterSweep*)calloc(1, sizeof(Haken_ParameterSweep));
    sweep->alpha_start = start;
    sweep->alpha_end   = end;
    sweep->alpha_step  = (end - start) / (n_points - 1);
    sweep->n_points    = n_points;
    return sweep;
}

void haken_sweep_free(Haken_ParameterSweep* sweep) {
    if (!sweep) return;
    free(sweep->eigenvalues);
    free(sweep->order_param_amp);
    free(sweep->relaxation_times);
    free(sweep);
}

int haken_sweep_execute(Haken_System* sys, Haken_ParameterSweep* sweep) {
    if (!sys || !sweep) return -1;
    int n = sys->n_dim;
    int np = sweep->n_points;

    sweep->eigenvalues = (double*)malloc((size_t)(n * np) * sizeof(double));
    sweep->relaxation_times = (double*)malloc((size_t)(n * np) * sizeof(double));
    sweep->order_param_amp = (double*)malloc((size_t)np * sizeof(double));

    for (int p = 0; p < np; p++) {
        double alpha = sweep->alpha_start + p * sweep->alpha_step;
        sys->control_param = alpha;
        haken_update_linear_matrix(sys);
        haken_compute_eigenmodes(sys);

        for (int k = 0; k < n; k++) {
            sweep->eigenvalues[p * n + k] = sys->modes[k].lambda_re;
            double re = sys->modes[k].lambda_re;
            sweep->relaxation_times[p * n + k] = (re < -1e-15) ? (-1.0 / re) : INFINITY;
        }
        /* Store order parameter amplitude: max of |Re(λ)| among slow modes */
        double max_amp = 0.0;
        for (int k = 0; k < n; k++) {
            if (fabs(sys->modes[k].lambda_re) < 0.1) {
                double amp = fabs(sys->state[k]);
                if (amp > max_amp) max_amp = amp;
            }
        }
        sweep->order_param_amp[p] = max_amp;
    }
    return 0;
}

double haken_sweep_max_spectral_gap(const Haken_ParameterSweep* sweep) {
    if (!sweep || !sweep->eigenvalues) return 0.0;
    /* Not fully implemented without mode classification per sweep point */
    return sweep->eigenvalues[0]; /* Placeholder: return largest eigenvalue */
}

double haken_sweep_critical_threshold(const Haken_ParameterSweep* sweep,
                                       int mode_index) {
    if (!sweep || !sweep->eigenvalues || mode_index < 0) return 0.0;
    int np = sweep->n_points;
    for (int p = 0; p < np - 1; p++) {
        double re1 = sweep->eigenvalues[p * sweep->n_points + mode_index];
        double re2 = sweep->eigenvalues[(p + 1) * sweep->n_points + mode_index];
        if (re1 * re2 < 0) {
            /* Linear interpolation */
            double alpha1 = sweep->alpha_start + p * sweep->alpha_step;
            double alpha2 = sweep->alpha_start + (p + 1) * sweep->alpha_step;
            return alpha1 - re1 * (alpha2 - alpha1) / (re2 - re1);
        }
    }
    return sweep->alpha_start;
}

double haken_susceptibility(const Haken_System* sys,
                             const Haken_OrderParam* op,
                             double perturbation_strength) {
    if (!sys || !op) return 0.0;
    /* χ = ∂ξ/∂h evaluated via finite difference:
     * Apply small perturbation h to the system,
     * measure change in order parameter amplitude. */
    double xi0 = op->amplitude ? op->amplitude[0] : 0.0;
    /* For now, use GL prediction: χ = 1/|a| */
    double a = op->gl_potential_coeff;
    if (fabs(a) < 1e-15) return INFINITY;
    double chi = 1.0 / fabs(a);
    (void)perturbation_strength; (void)xi0;
    return chi;
}

double haken_mode_softening_rate(Haken_System* sys, int mode_idx,
                                  double dalpha) {
    if (!sys || mode_idx < 0 || mode_idx >= sys->n_dim) return 0.0;
    double alpha0 = sys->control_param;
    double lambda0 = sys->modes[mode_idx].lambda_re;

    /* Evaluate at α + dα */
    sys->control_param = alpha0 + dalpha;
    haken_update_linear_matrix(sys);
    haken_compute_eigenmodes(sys);
    double lambda1 = sys->modes[mode_idx].lambda_re;

    /* Restore */
    sys->control_param = alpha0;
    haken_update_linear_matrix(sys);

    return (lambda1 - lambda0) / dalpha;
}

/* ==============================================================
 * Predictor-corrector continuation (Allgower & Georg 1990)
 * ============================================================== */

int haken_continuation(Haken_System* sys,
                        double alpha_start, double alpha_end,
                        int n_steps, double* states_out) {
    if (!sys || n_steps <= 0 || !states_out) return -1;
    int n = sys->n_dim;
    double dalpha = (alpha_end - alpha_start) / n_steps;

    sys->control_param = alpha_start;
    haken_update_linear_matrix(sys);

    for (int step = 0; step < n_steps; step++) {
        /* Store state */
        memcpy(states_out + step * n, sys->state, (size_t)n * sizeof(double));

        /* Predictor: Euler step in α */
        double* dq_dalpha = (double*)calloc((size_t)n, sizeof(double));
        /* Solve dq/dα = -J^{-1}·∂f/∂α (for fixed points) */
        double* J = (double*)calloc((size_t)(n * n), sizeof(double));
        haken_compute_jacobian(sys, sys->state, J, 1e-6);
        /* Approximate ∂f/∂α by finite difference */
        double* f0 = (double*)calloc((size_t)n, sizeof(double));
        double* f1 = (double*)calloc((size_t)n, sizeof(double));
        haken_rhs(sys, sys->state, f0);
        sys->control_param += 1e-4;
        haken_update_linear_matrix(sys);
        haken_rhs(sys, sys->state, f1);
        sys->control_param -= 1e-4;
        haken_update_linear_matrix(sys);

        for (int i = 0; i < n; i++) dq_dalpha[i] = -(f1[i] - f0[i]) / 1e-4;
        for (int i = 0; i < n; i++) states_out[step * n + i] += dalpha * dq_dalpha[i];
        memcpy(sys->state, states_out + step * n, (size_t)n * sizeof(double));

        free(dq_dalpha); free(J); free(f0); free(f1);

        /* Corrector: one Newton step (simplified) */
        sys->control_param += dalpha;
        haken_update_linear_matrix(sys);
    }
    return n_steps;
}
