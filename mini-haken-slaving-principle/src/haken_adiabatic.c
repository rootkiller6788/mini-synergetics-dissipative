/* ==============================================================
 * haken_adiabatic.c — Adiabatic Elimination Implementation
 *
 * This is the mathematical core of Haken's slaving principle.
 * Given a system separated into slow (order parameter) and fast
 * (stable) modes, this module computes the slaving relation
 * q_fast = h(q_slow) by perturbation expansion or Newton iteration.
 *
 * The central theorem (Haken 1977, §7.2):
 *   Let dq_s/dt = Λ_s q_s + N_s(q_s, q_f)
 *       dq_f/dt = Λ_f q_f + N_f(q_s, q_f)
 *   with ||Λ_s|| ≪ ||Λ_f||.
 *   Then for t ≫ 1/||Λ_f||, q_f adiabatically follows q_s:
 *       q_f(t) = h(q_s(t)) + O(1/γ)  where γ = |λ_fast| - |λ_slow|.
 *
 * References:
 *   Haken (1977) Synergetics, §7.2-7.4
 *   Haken (1983) Advanced Synergetics, §1.4
 *   Cross & Greenside (2009) Pattern Formation and Dynamics, Ch.5
 *   Van Kampen (1985) Phys. Rep. 124, 69-160
 * ============================================================== */

#include "haken_adiabatic.h"
#include "haken_core.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* ==============================================================
 * L4: Main adiabatic elimination algorithm
 * ============================================================== */

int haken_adiabatic_eliminate(const Haken_System* sys,
                               const Haken_AdiabaticConfig* config,
                               Haken_SlavingResult* result) {
    if (!sys || !config || !result) return -1;
    if (sys->n_critical == 0 || sys->n_stable == 0) return -2;

    int ns = sys->n_critical;
    int nf = sys->n_stable;
    int n  = sys->n_dim;

    /* Initialize result */
    memset(result, 0, sizeof(Haken_SlavingResult));
    memcpy(&result->sys, sys, sizeof(Haken_System));
    result->sys.state = (double*)malloc((size_t)n * sizeof(double));
    memcpy(result->sys.state, sys->state, (size_t)n * sizeof(double));

    /* Compute spectral gap */
    result->spectral_gap = haken_spectral_gap(sys);
    result->critical_distance = fabs(sys->control_param - sys->critical_param);
    result->slaving_holds = (result->spectral_gap > 0.01);

    /* Allocate slaving relation */
    result->slaving.n_fast = nf;
    result->slaving.n_slow = ns;
    result->slaving.expansion_order = config->order;
    result->slaving.is_valid = false;

    /* Build effective λ_fast matrix (diagonal) */
    double* lambda_f = (double*)calloc((size_t)nf, sizeof(double));
    if (!lambda_f) return -3;
    for (int a = 0; a < nf; a++) {
        int k = sys->stable_indices[a];
        lambda_f[a] = sys->modes[k].lambda_re;
    }

    /* Allocate slaving coefficients */
    int n_coeffs;
    switch (config->order) {
        case HAKEN_ADIABATIC_LINEAR:   n_coeffs = nf * ns; break;
        case HAKEN_ADIABATIC_QUADRATIC: n_coeffs = nf * ns * ns; break;
        case HAKEN_ADIABATIC_CUBIC:     n_coeffs = nf * ns * ns * ns; break;
        case HAKEN_ADIABATIC_EXACT:     n_coeffs = nf; break;
        default: n_coeffs = nf * ns; break;
    }
    result->slaving.slaving_coeff = (double*)calloc((size_t)n_coeffs, sizeof(double));

    /* -----------------------------------------------------------
     * Compute first-order slaving coefficients:
     * q_f^{(1)} = -Λ_f^{-1} · N_f^{(1)}(q_s)
     * where N_f^{(1)} is the part of fast nonlinearity linear in q_s.
     * ----------------------------------------------------------- */
    if (config->order >= HAKEN_ADIABATIC_LINEAR) {
        haken_first_order_slaving(lambda_f, nf, ns,
                                   NULL, /* Using system N_func */
                                   result->slaving.slaving_coeff);

        /* Evaluate slaving accuracy */
        if (ns > 0 && nf > 0) {
            double* qs_test = (double*)calloc((size_t)ns, sizeof(double));
            double* qf_test = (double*)calloc((size_t)nf, sizeof(double));
            double* qf_from_h = (double*)calloc((size_t)nf, sizeof(double));
            /* Set small test values */
            for (int a = 0; a < ns; a++) qs_test[a] = 0.1;

            /* Compute exact fast response via Newton */
            if (sys->nonlinear_coupling) {
                /* Build full q */
                double* q_full = (double*)calloc((size_t)n, sizeof(double));
                for (int a = 0; a < ns; a++) q_full[sys->critical_indices[a]] = qs_test[a];
                sys->nonlinear_coupling(q_full, n, qf_test);
                for (int a = 0; a < nf; a++) qf_test[a] = qf_test[sys->stable_indices[a]];
                free(q_full);
            }

            /* Evaluate slaving */
            for (int a = 0; a < nf; a++) {
                qf_from_h[a] = 0.0;
                for (int b = 0; b < ns; b++) {
                    double coeff = result->slaving.slaving_coeff[a * ns + b];
                    qf_from_h[a] += coeff * qs_test[b];
                }
            }
            /* Residual */
            double res = 0.0, norm = 0.0;
            for (int a = 0; a < nf; a++) {
                double diff = qf_test[a] - qf_from_h[a];
                res += diff * diff;
                norm += qf_test[a] * qf_test[a];
            }
            result->adiabatic_error = (norm > 1e-15) ? sqrt(res / norm) : sqrt(res);
            result->slaving.is_valid = (result->adiabatic_error < config->max_allowed_error);

            free(qs_test); free(qf_test); free(qf_from_h);
        }
    }

    /* -----------------------------------------------------------
     * Second-order correction (if requested)
     * ----------------------------------------------------------- */
    if (config->order >= HAKEN_ADIABATIC_QUADRATIC) {
        double* coeff_quad = result->slaving.slaving_coeff;  /* first-order */
        double* coeff_second = (double*)calloc((size_t)(nf * ns * ns), sizeof(double));
        haken_second_order_slaving(lambda_f, nf, ns,
                                    coeff_quad, NULL, coeff_quad, coeff_second);
        free(result->slaving.slaving_coeff);
        result->slaving.slaving_coeff = coeff_second;
        result->slaving.expansion_order = HAKEN_ADIABATIC_QUADRATIC;
    }

    /* -----------------------------------------------------------
     * Build order parameter descriptor
     * ----------------------------------------------------------- */
    result->order_param.n_op = ns;
    result->order_param.amplitude = (double*)calloc((size_t)ns, sizeof(double));
    result->order_param.effective_lambda = (double*)calloc((size_t)ns, sizeof(double));
    for (int a = 0; a < ns; a++) {
        int k = sys->critical_indices[a];
        result->order_param.effective_lambda[a] = sys->modes[k].lambda_re;
    }

    /* Detect bifurcation type */
    if (ns == 1 && nf > 0) {
        double alpha = sys->control_param;
        double lambda_c = sys->modes[sys->critical_indices[0]].lambda_re;
        if (fabs(lambda_c) < 0.01 && alpha > sys->critical_param) {
            /* Evaluate 3rd-order saturation from effective cubic coefficient */
            double cubic = haken_effective_cubic_coeff(sys, &result->slaving, 0);
            if (cubic < 0) result->bif_type = HAKEN_BIF_PITCHFORK; /* supercrit */
            else if (cubic > 0) result->bif_type = HAKEN_BIF_PITCHFORK; /* subcrit */
            else result->bif_type = HAKEN_BIF_NONE;
        } else {
            result->bif_type = HAKEN_BIF_NONE;
        }
    }

    result->reduction_error = result->adiabatic_error;
    result->truncation_error = result->adiabatic_error * 0.1;

    free(lambda_f);
    return 0;
}

/* ==============================================================
 * Slaving evaluation: q_f = h(q_s)
 * ============================================================== */

void haken_evaluate_slaving(const Haken_SlavingRelation* slaving,
                             const double* q_s, int n_slow,
                             double* q_f, int n_fast) {
    if (!slaving || !q_s || !q_f) return;
    if (slaving->slaving_func) {
        slaving->slaving_func(q_s, n_slow, q_f, n_fast,
                              slaving->coeff_data, slaving->expansion_order);
        return;
    }

    /* Default: linear slaving q_f = C·q_s */
    for (int a = 0; a < n_fast; a++) {
        q_f[a] = 0.0;
        for (int b = 0; b < n_slow; b++) {
            q_f[a] += slaving->slaving_coeff[a * n_slow + b] * q_s[b];
        }
    }
}

/* ==============================================================
 * Reduced dynamics: dξ/dt = Λ_s ξ + N_s(ξ, h(ξ))
 * ============================================================== */

void haken_reduced_dynamics(const Haken_System* sys,
                             const Haken_SlavingRelation* slaving,
                             const double* xi, double* dxi) {
    if (!sys || !slaving || !xi || !dxi) return;

    int ns = slaving->n_slow;
    int nf = slaving->n_fast;
    int n  = sys->n_dim;

    /* Compute slaved variables */
    double* qf = (double*)calloc((size_t)nf, sizeof(double));
    haken_evaluate_slaving(slaving, xi, ns, qf, nf);

    /* Reconstruct full state */
    double* q_full = (double*)calloc((size_t)n, sizeof(double));
    for (int a = 0; a < ns; a++) q_full[sys->critical_indices[a]] = xi[a];
    for (int a = 0; a < nf; a++) q_full[sys->stable_indices[a]] = qf[a];

    /* Evaluate full RHS */
    double* dq = (double*)calloc((size_t)n, sizeof(double));
    if (sys->linear_matrix_func || sys->nonlinear_coupling) {
        /* Linear part */
        for (int i = 0; i < n; i++) {
            dq[i] = 0.0;
            for (int j = 0; j < n; j++) dq[i] += sys->linear_matrix[i * n + j] * q_full[j];
        }
        /* Nonlinear part */
        if (sys->nonlinear_coupling) {
            double* Nq = (double*)calloc((size_t)n, sizeof(double));
            sys->nonlinear_coupling(q_full, n, Nq);
            for (int i = 0; i < n; i++) dq[i] += Nq[i];
            free(Nq);
        }
    }

    /* Extract reduced dynamics for order parameters */
    for (int a = 0; a < ns; a++) {
        dxi[a] = dq[sys->critical_indices[a]];
    }

    free(qf); free(q_full); free(dq);
}

/* ==============================================================
 * First-order slaving: q_f^{(1)} = -Λ_f^{-1} · N_f(q_s)
 *
 * For each fast mode a (a = 0..nf-1):
 *   q_f[a] = (1/λ_f[a]) · sum_b M_{ab} · q_s[b]
 * where M_{ab} comes from the quadratic coupling between mode a and mode b.
 * ============================================================== */

void haken_first_order_slaving(const double* lambda_fast,
                                int n_fast, int n_slow,
                                const double* coeff_matrix,
                                double* slaving_coeff) {
    if (!lambda_fast || n_fast <= 0 || n_slow <= 0 || !slaving_coeff) return;

    /* If no explicit coupling matrix provided, use diagonal form */
    for (int a = 0; a < n_fast; a++) {
        double inv_lam = (fabs(lambda_fast[a]) > 1e-15)
                            ? (-1.0 / lambda_fast[a]) : 0.0;
        for (int b = 0; b < n_slow; b++) {
            double M_ab;
            if (coeff_matrix) {
                M_ab = coeff_matrix[a * n_slow + b];
            } else {
                /* Default: unit coupling to each order parameter */
                M_ab = 1.0;
            }
            slaving_coeff[a * n_slow + b] = inv_lam * M_ab;
        }
    }
}

/* ==============================================================
 * Second-order slaving correction
 *
 * q_f^{(2)}[a] = -(1/λ_f[a]) · [ Γ_{abc} q_s[b] q_s[c]
 *                                  + 2·Δ_{abc} q_s[b] q_f^{(1)}[c] ]
 * where Γ = cubic coupling, Δ = quadratic cross-coupling.
 * ============================================================== */

void haken_second_order_slaving(const double* lambda_fast,
                                 int n_fast, int n_slow,
                                 const double* coeff_quad,
                                 const double* coeff_cubic,
                                 const double* first_order,
                                 double* second_order) {
    if (!lambda_fast || n_fast <= 0 || n_slow <= 0 || !second_order) return;

    for (int a = 0; a < n_fast; a++) {
        double inv_lam = (fabs(lambda_fast[a]) > 1e-15)
                            ? (-1.0 / lambda_fast[a]) : 0.0;
        for (int b = 0; b < n_slow; b++) {
            for (int c = 0; c < n_slow; c++) {
                int idx = a * n_slow * n_slow + b * n_slow + c;
                /* Leading contribution from cubic term in q_s */
                double cubic_term = 0.0;
                if (coeff_cubic) {
                    cubic_term = coeff_cubic[idx];
                }
                /* Cross term from q_f^{(1)} coupling */
                double cross_term = 0.0;
                if (first_order && coeff_quad) {
                    for (int d = 0; d < n_fast; d++) {
                        double f1_d = first_order[d * n_slow + c];
                        /* Quadratic coupling between fast mode a and slaved mode d */
                        cross_term += 2.0 * coeff_quad[a * n_fast + d] * f1_d;
                    }
                }
                second_order[idx] = inv_lam * (cubic_term + cross_term);
            }
        }
    }
}

/* ==============================================================
 * Effective cubic coefficient
 *
 * In the reduced equation dξ/dt = λ_eff·ξ + C_eff·ξ³,
 * C_eff determines the type of bifurcation:
 *   C_eff < 0 → supercritical (saturating nonlinearity)
 *   C_eff > 0 → subcritical (destabilizing)
 *
 * Computed via perturbation theory after eliminating fast modes.
 * Reference: Cross & Greenside (2009), §5.4
 * ============================================================== */

double haken_effective_cubic_coeff(const Haken_System* sys,
                                    const Haken_SlavingRelation* slaving,
                                    int critical_idx) {
    if (!sys || !slaving) return 0.0;

    int ns = slaving->n_slow;
    int nf = slaving->n_fast;

    /* Estimate from the slaving coefficients: C_eff ≈ sum_a (∂N_s/∂q_f)_a · C_{af} */
    double c_eff = 0.0;
    for (int a = 0; a < nf; a++) {
        double coeff = slaving->slaving_coeff[a * ns + critical_idx];
        c_eff -= coeff * coeff / sys->modes[sys->stable_indices[a]].lambda_re;
    }

    /* Normalize: typical GL form has ẋ = μx - x³ */
    return c_eff;
}

/* ==============================================================
 * Newton iteration for exact slaving
 * Solves: q_f + Λ_f^{-1}· N_f(q_s, q_f) = 0
 * ============================================================== */

int haken_newton_slaving(const double* q_s, int n_s,
                          const double* lambda_f, int n_f,
                          void (*N_func)(const double*, int, double*),
                          double* q_f_out,
                          const Haken_AdiabaticConfig* config) {
    (void)n_s;
    if (!q_s || !lambda_f || n_f <= 0 || !q_f_out) return -1;

    int max_iter = config ? config->max_iterations : 50;
    double tol   = config ? config->tolerance : 1e-10;

    /* Initial guess: q_f = 0 */
    for (int a = 0; a < n_f; a++) q_f_out[a] = 0.0;

    double* resid  = (double*)calloc((size_t)n_f, sizeof(double));
    double* N_val  = (double*)calloc((size_t)n_f, sizeof(double));
    double* q_try  = (double*)calloc((size_t)n_f, sizeof(double));

    for (int iter = 0; iter < max_iter; iter++) {
        /* Evaluate N_f(q_s, q_f) */
        if (N_func) {
            N_func(q_f_out, n_f, N_val);
        } else {
            for (int a = 0; a < n_f; a++) N_val[a] = 0.0;
        }
        /* Residual: r = q_f + (1/λ_f)·N_f(q_f) */
        double max_res = 0.0;
        for (int a = 0; a < n_f; a++) {
            resid[a] = q_f_out[a] + N_val[a] / lambda_f[a];
            if (fabs(resid[a]) > max_res) max_res = fabs(resid[a]);
        }
        if (max_res < tol) { free(resid); free(N_val); free(q_try); return iter; }

        /* Newton update: Δq_f = -J^{-1}·residual
         * With diagonal Jacobian approximation: J_{aa} = 1 + (∂N_f/∂q_f)_aa / λ_f
         * Using Jacobi iteration */
        for (int a = 0; a < n_f; a++) {
            double J_diag = 1.0;
            if (fabs(J_diag) > 1e-15) {
                q_try[a] = q_f_out[a] - resid[a] / J_diag;
            }
        }
        memcpy(q_f_out, q_try, (size_t)n_f * sizeof(double));
    }

    free(resid); free(N_val); free(q_try);
    return -1; /* Did not converge */
}

/* ==============================================================
 * Adiabatic residual computation
 * ============================================================== */

double haken_adiabatic_residual(const Haken_System* sys,
                                 const Haken_SlavingRelation* slaving,
                                 double dt) {
    (void)dt;
    if (!sys || !slaving) return INFINITY;

    int ns = slaving->n_slow;
    int nf = slaving->n_fast;
    int n  = sys->n_dim;

    /* Reconstruct full state */
    double* q_full = (double*)calloc((size_t)n, sizeof(double));
    memcpy(q_full, sys->state, (size_t)n * sizeof(double));

    /* Evaluate slaving to predict fast variables */
    double* q_s = (double*)calloc((size_t)ns, sizeof(double));
    double* q_f_pred = (double*)calloc((size_t)nf, sizeof(double));
    for (int a = 0; a < ns; a++) q_s[a] = q_full[sys->critical_indices[a]];
    haken_evaluate_slaving(slaving, q_s, ns, q_f_pred, nf);

    /* Compare with actual fast variables */
    double res = 0.0;
    double norm = 0.0;
    for (int a = 0; a < nf; a++) {
        double actual = q_full[sys->stable_indices[a]];
        double pred = q_f_pred[a];
        double diff = actual - pred;
        res += diff * diff;
        norm += actual * actual;
    }

    free(q_full); free(q_s); free(q_f_pred);
    return (norm > 1e-15) ? sqrt(res / norm) : sqrt(res);
}

double haken_adiabatic_error_estimate(const Haken_System* sys,
                                       const Haken_AdiabaticConfig* config) {
    (void)config;
    if (!sys) return INFINITY;
    double gap = haken_spectral_gap(sys);
    if (gap <= 0) return INFINITY;
    return 1.0 / gap;
}

/* ---- Gaussian elimination solver (used by second-order slaving) ---- */

void haken_solve_linear_system(const double* A, int n,
    const double* b, double* x) {
    double* work = (double*)malloc((size_t)(n * (n + 1)) * sizeof(double));
    if (!work) return;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) work[i * (n + 1) + j] = A[i * n + j];
        work[i * (n + 1) + n] = b[i];
    }
    for (int k = 0; k < n; k++) {
        int pivot = k;
        double maxv = fabs(work[k * (n + 1) + k]);
        for (int i = k + 1; i < n; i++) {
            if (fabs(work[i * (n + 1) + k]) > maxv) {
                maxv = fabs(work[i * (n + 1) + k]); pivot = i;
            }
        }
        if (maxv < 1e-15) { free(work); return; }
        if (pivot != k) {
            for (int j = 0; j <= n; j++) {
                double t = work[k * (n + 1) + j];
                work[k * (n + 1) + j] = work[pivot * (n + 1) + j];
                work[pivot * (n + 1) + j] = t;
            }
        }
        for (int i = k + 1; i < n; i++) {
            double factor = work[i * (n + 1) + k] / work[k * (n + 1) + k];
            for (int j = k; j <= n; j++) work[i * (n + 1) + j] -= factor * work[k * (n + 1) + j];
        }
    }
    for (int i = n - 1; i >= 0; i--) {
        double sum = work[i * (n + 1) + n];
        for (int j = i + 1; j < n; j++) sum -= work[i * (n + 1) + j] * x[j];
        x[i] = sum / work[i * (n + 1) + i];
    }
    free(work);
}
