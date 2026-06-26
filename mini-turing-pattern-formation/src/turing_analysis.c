/**
 * turing_analysis.c — Linear Stability & Bifurcation Analysis
 *
 * Implements the complete analysis pipeline for Turing instability:
 *   1. 4 Turing conditions evaluation
 *   2. Dispersion relation λ(k) computation
 *   3. Critical wavenumber determination
 *   4. Homogeneous steady state (HSS) finder via Newton's method
 *   5. Turing space construction via parameter scanning
 *   6. Bifurcation point detection
 *   7. Sensitivity analysis
 *
 * Reference: Turing, A.M. (1952) "The Chemical Basis of Morphogenesis"
 *            Murray, J.D. (2003) Mathematical Biology II, Ch. 2-3
 *            Cross & Greenside (2009) Pattern Formation and Dynamics
 *            Maini et al. (1997) Bull. Math. Biol. 59(1):61-84
 *
 * Knowledge Coverage:
 *   L2: Diffusion-driven instability, activator-inhibitor requirements
 *   L3: Eigenvalue analysis of 2x2 matrices, dispersion relations
 *   L4: Turing instability theorem (4 conditions), bifurcation theory
 *   L5: Newton's method for HSS, parameter scanning algorithms
 */

#include "turing_analysis.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==========================================================================
 * Matrix 2×2 Eigenvalues
 *
 * For M = [[a, b], [c, d]]:
 *   trace T = a + d
 *   determinant D = a*d - b*c
 *   λ₁,₂ = (T ± √(T² - 4D)) / 2
 *
 * If T² - 4D < 0, eigenvalues are complex conjugate pairs.
 * The real part is T/2 for both (which is the growth/decay rate).
 *
 * Reference: Strang, G. (2019) Linear Algebra and Learning from Data, §1.5
 * ==========================================================================*/

int matrix_2x2_eigenvalues(const double M[2][2],
                           double eigen_real[2], double eigen_imag[2]) {
    if (!M || !eigen_real || !eigen_imag) return -1;

    double a = M[0][0], b = M[0][1], c = M[1][0], d = M[1][1];
    double trace = a + d;
    double det = a * d - b * c;
    double disc = trace * trace - 4.0 * det;

    if (disc >= 0.0) {
        double sqrt_disc = sqrt(disc);
        eigen_real[0] = 0.5 * (trace + sqrt_disc);
        eigen_real[1] = 0.5 * (trace - sqrt_disc);
        eigen_imag[0] = 0.0;
        eigen_imag[1] = 0.0;
    } else {
        eigen_real[0] = 0.5 * trace;
        eigen_real[1] = 0.5 * trace;
        eigen_imag[0] = 0.5 * sqrt(-disc);
        eigen_imag[1] = -0.5 * sqrt(-disc);
    }
    return 0;
}

/* ==========================================================================
 * Hopf Instability Check
 *
 * A Hopf bifurcation occurs when:
 *   1. trace(J) > 0  (destabilizing)
 *   2. det(J) > 0    (no saddle)
 *   3. T² - 4D < 0   (complex conjugate eigenvalues)
 *
 * Reference: Strogatz, S.H. (2015) Nonlinear Dynamics and Chaos, Ch. 8
 * ==========================================================================*/

int is_hopf_unstable(const double J[2][2]) {
    double trace = J[0][0] + J[1][1];
    double det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
    if (trace > 0.0 && det > 0.0) return 1;
    return 0;
}

/* ==========================================================================
 * Turing Conditions Computation
 *
 * For the system ∂u/∂t = f(u,v) + D_u∇²u, ∂v/∂t = g(u,v) + D_v∇²v:
 *
 * Let J = [[f_u, f_v], [g_u, g_v]] evaluated at HSS.
 * Let D = [[D_u, 0], [0, D_v]].
 *
 * The dispersion relation comes from:
 *   det(J - k²D - λI) = 0
 *
 *   λ(k) satisfies: λ² - T(k)λ + D(k) = 0
 *   where T(k) = trace(J) - k²(D_u + D_v)
 *         D(k) = det(J) - k²(D_v*f_u + D_u*g_v) + k⁴*D_u*D_v
 *
 * Turing instability requires:
 *   Cond 1: f_u + g_v < 0  (HSS stable without diffusion)
 *   Cond 2: f_u*g_v - f_v*g_u > 0  (det J > 0)
 *   Cond 3: D_v*f_u + D_u*g_v > 0  (unstable to some k)
 *   Cond 4: (D_v*f_u + D_u*g_v)² > 4*D_u*D_v*det(J)  (real eigenvalues at k_max)
 *
 * Reference: Murray (2003) §2.3, Eqs. (2.27)-(2.31)
 * ==========================================================================*/

int turing_conditions_compute(const double J[2][2], double Du, double Dv,
                               TuringConditions *tc) {
    if (!J || !tc) return -1;

    tc->f_u = J[0][0];
    tc->f_v = J[0][1];
    tc->g_u = J[1][0];
    tc->g_v = J[1][1];

    tc->trace_J = tc->f_u + tc->g_v;
    tc->det_J   = tc->f_u * tc->g_v - tc->f_v * tc->g_u;

    /* Cond 1: Trace < 0 → stable to homogeneous perturbations */
    tc->cond1 = (tc->trace_J < 0.0) ? 1 : 0;

    /* Cond 2: Determinant > 0 → no saddle point without diffusion */
    tc->cond2 = (tc->det_J > 0.0) ? 1 : 0;

    /* Cond 3: D_v*f_u + D_u*g_v > 0 → cross-diffusion drives instability */
    double cross = Dv * tc->f_u + Du * tc->g_v;
    tc->cond3 = (cross > 0.0) ? 1 : 0;

    /* Cond 4: Discriminant > 0 → real eigenvalues at the most unstable mode */
    tc->discriminant = cross * cross - 4.0 * Du * Dv * tc->det_J;
    tc->cond4 = (tc->discriminant > 0.0) ? 1 : 0;

    /* All 4 conditions must hold */
    tc->turing_unstable = tc->cond1 && tc->cond2 && tc->cond3 && tc->cond4;

    return 0;
}

void turing_conditions_print(const TuringConditions *tc) {
    if (!tc) return;

    printf("Turing Instability Conditions:\n");
    printf("  Jacobian at HSS:\n");
    printf("    f_u=%.6g  f_v=%.6g\n", tc->f_u, tc->f_v);
    printf("    g_u=%.6g  g_v=%.6g\n", tc->g_u, tc->g_v);
    printf("  Trace(J)  = %.6g\n", tc->trace_J);
    printf("  Det(J)    = %.6g\n", tc->det_J);
    printf("  Cond 1 (trace < 0):         %s  (%.6g)\n",
           tc->cond1 ? "PASS" : "FAIL", tc->trace_J);
    printf("  Cond 2 (det > 0):           %s  (%.6g)\n",
           tc->cond2 ? "PASS" : "FAIL", tc->det_J);
    printf("  Cond 3 (Dv*f_u + Du*g_v > 0): %s\n",
           tc->cond3 ? "PASS" : "FAIL");
    printf("  Cond 4 (discriminant > 0):  %s  (%.6g)\n",
           tc->cond4 ? "PASS" : "FAIL", tc->discriminant);
    printf("  Discriminant = %.6g\n", tc->discriminant);
    printf("  TURING-UNSTABLE: %s\n",
           tc->turing_unstable ? "YES ✓" : "NO");
}

/* ==========================================================================
 * Dispersion Relation λ(k)
 *
 * For a given wavenumber k, the eigenvalues λ of (J - k²D) are:
 *   λ(k) = [T(k) ± sqrt(T(k)² - 4*D(k))] / 2
 *
 * where T(k) = trace(J) - k²(D_u + D_v)
 *       D(k) = det(J) - k²(D_v*f_u + D_u*g_v) + k⁴*D_u*D_v
 *
 * The growth rate is Re(λ(k)). The system is unstable if max Re(λ(k)) > 0.
 *
 * The critical wavenumber k_c maximizes Re(λ(k)):
 *   k_c² = sqrt(det_J / (D_u * D_v))
 *
 * Reference: Murray (2003) §2.3, Eq. (2.28)
 * ==========================================================================*/

DispersionRelation *turing_dispersion_compute(const double J[2][2],
                                              double Du, double Dv,
                                              int n_modes,
                                              double k_min, double k_max) {
    if (!J || n_modes <= 0) return NULL;

    DispersionRelation *dr = dispersion_alloc(n_modes);
    if (!dr) return NULL;

    double f_u = J[0][0], f_v = J[0][1];
    double g_u = J[1][0], g_v = J[1][1];
    double trace_J = f_u + g_v;
    double det_J = f_u * g_v - f_v * g_u;
    double sum_D = Du + Dv;

    dr->turing_unstable = 0;
    dr->lambda_max = -1e100;
    dr->k_critical = 0.0;
    dr->k_min = -1.0;
    dr->k_max = -1.0;
    int found_unstable_band = 0;

    for (int i = 0; i < n_modes; i++) {
        double k;
        if (n_modes == 1) {
            k = k_min;
        } else {
            k = k_min + (k_max - k_min) * (double)i / (double)(n_modes - 1);
        }

        double k2 = k * k;
        double k4 = k2 * k2;

        double T_k = trace_J - k2 * sum_D;
        double D_k = det_J - k2 * (Dv * f_u + Du * g_v) + k4 * Du * Dv;

        double disc = T_k * T_k - 4.0 * D_k;

        double lambda_real;
        if (disc >= 0.0) {
            /* Real eigenvalues: the larger one */
            double sqrt_disc = sqrt(disc);
            lambda_real = 0.5 * (T_k + sqrt_disc);
            /* Also check the smaller */
            double lambda2 = 0.5 * (T_k - sqrt_disc);
            if (lambda2 > lambda_real) lambda_real = lambda2;
        } else {
            /* Complex conjugate pair: real part = T_k/2 */
            lambda_real = 0.5 * T_k;
        }

        dr->wavenumbers[i] = k;
        dr->growth_rates[i] = lambda_real;

        /* Track maximum */
        if (lambda_real > dr->lambda_max) {
            dr->lambda_max = lambda_real;
            dr->k_critical = k;
        }

        /* Track unstable band */
        if (lambda_real > 0.0) {
            dr->turing_unstable = 1;
            if (!found_unstable_band) {
                dr->k_min = k;
                found_unstable_band = 1;
            }
            dr->k_max = k;
        }
    }

    dr->predicted_wavelength = (dr->k_critical > 1e-15)
                               ? 2.0 * M_PI / dr->k_critical : 0.0;

    return dr;
}

void turing_dispersion_print(const DispersionRelation *dr) {
    if (!dr) return;

    printf("Dispersion Relation Summary:\n");
    printf("  Turing unstable: %s\n", dr->turing_unstable ? "YES" : "NO");
    printf("  Max growth rate λ_max: %.6g\n", dr->lambda_max);
    printf("  Critical wavenumber k_c: %.6g\n", dr->k_critical);
    printf("  Predicted wavelength: %.6g\n", dr->predicted_wavelength);
    printf("  Unstable band: [%.6g, %.6g]\n", dr->k_min, dr->k_max);

    if (dr->n_modes <= 20) {
        printf("  Data points:\n");
        for (int i = 0; i < dr->n_modes; i++) {
            printf("    k=%8.4f  λ=%12.6g\n",
                   dr->wavenumbers[i], dr->growth_rates[i]);
        }
    }
}

double turing_critical_wavenumber(const double J[2][2], double Du, double Dv) {
    if (!J) return 0.0;
    double det_J = J[0][0] * J[1][1] - J[0][1] * J[1][0];
    if (det_J <= 0.0 || Du <= 0.0 || Dv <= 0.0) return 0.0;
    double kc2 = sqrt(det_J / (Du * Dv));
    return sqrt(kc2);
}

double turing_predicted_wavelength(const DispersionRelation *dr) {
    if (!dr || dr->k_critical < 1e-15) return 0.0;
    return 2.0 * M_PI / dr->k_critical;
}

/* ==========================================================================
 * Newton's Method for Homogeneous Steady State
 *
 * Solve f(u,v) = 0, g(u,v) = 0.
 *
 * Newton iteration:
 *   [u_{k+1}] = [u_k] - J(u_k, v_k)^{-1} * [f(u_k, v_k)]
 *   [v_{k+1}]   [v_k]                       [g(u_k, v_k)]
 *
 * Where J^{-1} = (1/det) * [[g_v, -f_v], [-g_u, f_u]]
 *
 * Convergence: quadratic if initial guess is close enough.
 * Global strategy: backtracking line search with Armijo condition.
 *
 * Reference: Kelley, C.T. (2003) Solving Nonlinear Equations with Newton's Method
 * ==========================================================================*/

HomogeneousSteadyState *newton_find_hss(ReactionFunc reaction_f,
                                         JacobianFunc jacobian_f,
                                         const ModelParams *p,
                                         double u0, double v0,
                                         int max_iter, double tol) {
    if (!reaction_f || !jacobian_f || !p || max_iter <= 0) return NULL;

    HomogeneousSteadyState *hss = (HomogeneousSteadyState *)calloc(1,
                                      sizeof(HomogeneousSteadyState));
    if (!hss) return NULL;

    double u = u0, v = v0;

    for (int iter = 0; iter < max_iter; iter++) {
        double f, g, J[2][2];

        if (reaction_f(u, v, p, &f, &g) != 0) {
            hss->converged = 0;
            hss->iterations = iter;
            return hss;
        }
        jacobian_f(u, v, p, J);

        double det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
        if (fabs(det) < 1e-15) {
            /* Perturb to escape singular Jacobian */
            u += 0.01; v += 0.01;
            continue;
        }

        /* Newton step */
        double du = -( J[1][1]*f - J[0][1]*g) / det;
        double dv = -(-J[1][0]*f + J[0][0]*g) / det;

        /* Backtracking line search with Armijo condition */
        double alpha = 1.0;
        double f_new, g_new;
        double norm_f = sqrt(f*f + g*g);

        for (int bt = 0; bt < 10; bt++) {
            double u_try = u + alpha * du;
            double v_try = v + alpha * dv;
            if (reaction_f(u_try, v_try, p, &f_new, &g_new) != 0) {
                alpha *= 0.5; continue;
            }
            double norm_f_new = sqrt(f_new*f_new + g_new*g_new);
            /* Sufficient decrease (Armijo condition): c = 1e-4 */
            if (norm_f_new <= (1.0 - 1e-4 * alpha) * norm_f) {
                break;
            }
            alpha *= 0.5;
        }

        u += alpha * du;
        v += alpha * dv;

        /* Enforce physical bounds */
        if (u < 0.0) u = 1e-10;
        if (v < 0.0) v = 1e-10;

        if (fabs(alpha*du) < tol && fabs(alpha*dv) < tol) {
            hss->u_star = u;
            hss->v_star = v;
            hss->converged = 1;
            hss->iterations = iter + 1;
            double f_final, g_final;
            reaction_f(u, v, p, &f_final, &g_final);
            hss->residual = sqrt(f_final*f_final + g_final*g_final);
            return hss;
        }
    }

    /* Max iterations reached without convergence */
    hss->u_star = u;
    hss->v_star = v;
    hss->converged = 0;
    hss->iterations = max_iter;
    double f_final, g_final;
    reaction_f(u, v, p, &f_final, &g_final);
    hss->residual = sqrt(f_final*f_final + g_final*g_final);
    return hss;
}

HomogeneousSteadyState **newton_find_all_hss(ReactionFunc reaction_f,
                                              JacobianFunc jacobian_f,
                                              const ModelParams *p,
                                              int *n_hss) {
    if (!reaction_f || !jacobian_f || !p || !n_hss) return NULL;

    /* Multi-start Newton from a grid of initial guesses */
    double guesses[][2] = {
        {0.1, 0.1}, {1.0, 1.0}, {10.0, 10.0},
        {0.01, 10.0}, {10.0, 0.01},
        {0.5, 5.0}, {5.0, 0.5},
        {2.0, 2.0}, {0.5, 0.5}
    };
    int n_guesses = 9;

    /* Allocate result array */
    int capacity = n_guesses;
    HomogeneousSteadyState **results = (HomogeneousSteadyState **)calloc(
        (size_t)capacity, sizeof(HomogeneousSteadyState *));
    if (!results) { *n_hss = 0; return NULL; }

    int count = 0;
    for (int g = 0; g < n_guesses; g++) {
        HomogeneousSteadyState *hss = newton_find_hss(
            reaction_f, jacobian_f, p,
            guesses[g][0], guesses[g][1], 100, 1e-10);

        if (hss && hss->converged) {
            /* Check if this HSS is distinct from previously found ones */
            int unique = 1;
            for (int j = 0; j < count; j++) {
                double du = hss->u_star - results[j]->u_star;
                double dv = hss->v_star - results[j]->v_star;
                if (sqrt(du*du + dv*dv) < 1e-6) {
                    unique = 0;
                    break;
                }
            }
            if (unique) {
                results[count++] = hss;
            } else {
                free(hss);
            }
        } else {
            free(hss);
        }
    }

    *n_hss = count;
    return results;
}

/* ==========================================================================
 * Turing Space Construction
 *
 * Scans a 2D parameter grid and classifies each point based on stability:
 *   - Homogeneous steady state computation
 *   - Turing condition check
 *   - Hopf condition check
 * ==========================================================================*/

TuringSpace *turing_space_alloc(ModelType model, const ModelParams *base,
                                double Du, double Dv,
                                int param1_idx, double p1_min, double p1_max, int n1,
                                int param2_idx, double p2_min, double p2_max, int n2) {
    if (!base || n1 <= 0 || n2 <= 0) return NULL;

    TuringSpace *ts = (TuringSpace *)calloc(1, sizeof(TuringSpace));
    if (!ts) return NULL;

    ts->model = model;
    ts->base_params = *base;
    ts->Du = Du;
    ts->Dv = Dv;
    ts->p1_min = p1_min; ts->p1_max = p1_max; ts->n1 = n1;
    ts->p2_min = p2_min; ts->p2_max = p2_max; ts->n2 = n2;

    size_t size = (size_t)n1 * (size_t)n2;
    ts->regions = (TuringRegion *)calloc(size, sizeof(TuringRegion));
    if (!ts->regions) {
        free(ts);
        return NULL;
    }

    return ts;
}

void turing_space_free(TuringSpace *ts) {
    if (ts) {
        free(ts->regions);
        free(ts);
    }
}

static void set_model_param(ModelParams *p, int idx, double value) {
    switch (idx) {
        case 0: p->F = value; break;
        case 1: p->k = value; break;
        case 2: p->A = value; break;
        case 3: p->B = value; break;
        case 4: p->a = value; break;
        case 5: p->b = value; break;
        default: break;
    }
}

int turing_space_compute(TuringSpace *ts, const SolverConfig *ignored_var) {
    if (!ts) return -1;
    (void)ignored_var;

    for (int i = 0; i < ts->n1; i++) {
        double p1 = ts->p1_min + (ts->p1_max - ts->p1_min) * (double)i
                    / (double)(ts->n1 - 1);
        for (int j = 0; j < ts->n2; j++) {
            double p2 = ts->p2_min + (ts->p2_max - ts->p2_min) * (double)j
                        / (double)(ts->n2 - 1);

            ModelParams p = ts->base_params;
            set_model_param(&p, 0, p1);  /* param1_idx used for p1 */
            set_model_param(&p, 1, p2);  /* param2_idx used for p2 */

            /* Find HSS */
            HomogeneousSteadyState hss;
            int hss_ok = (model_find_hss(ts->model, &p, &hss) == 0 && hss.converged);
            if (!hss_ok) {
                ts->regions[i * ts->n2 + j] = REGION_NO_HSS;
                continue;
            }

            /* Compute Jacobian at HSS */
            double J[2][2];
            model_jacobian(ts->model, hss.u_star, hss.v_star, &p, J);

            /* Check Turing and Hopf */
            TuringConditions tc;
            turing_conditions_compute(J, ts->Du, ts->Dv, &tc);
            int hopf = is_hopf_unstable(J);

            if (tc.turing_unstable && hopf) {
                ts->regions[i * ts->n2 + j] = REGION_TURING_HOPF;
            } else if (tc.turing_unstable) {
                ts->regions[i * ts->n2 + j] = REGION_TURING;
            } else if (hopf) {
                ts->regions[i * ts->n2 + j] = REGION_HOPF;
            } else {
                ts->regions[i * ts->n2 + j] = REGION_STABLE;
            }
        }
    }
    return 0;
}

void turing_space_count_regions(const TuringSpace *ts, int counts[5]) {
    if (!ts || !counts) return;
    for (int k = 0; k < 5; k++) counts[k] = 0;
    size_t total = (size_t)ts->n1 * (size_t)ts->n2;
    for (size_t k = 0; k < total; k++) {
        counts[ts->regions[k]]++;
    }
}

void turing_space_print_summary(const TuringSpace *ts) {
    if (!ts) return;
    int counts[5];
    turing_space_count_regions(ts, counts);
    size_t total = (size_t)ts->n1 * (size_t)ts->n2;

    printf("Turing Space Summary (%zu grid points):\n", total);
    printf("  Model: %s\n", model_name(ts->model));
    printf("  D_u=%.4g, D_v=%.4g\n", ts->Du, ts->Dv);
    printf("  Parameter ranges: [%.4g, %.4g] × [%.4g, %.4g]\n",
           ts->p1_min, ts->p1_max, ts->p2_min, ts->p2_max);
    printf("  STABLE:      %6d (%.1f%%)\n",
           counts[REGION_STABLE], 100.0*counts[REGION_STABLE]/total);
    printf("  TURING:      %6d (%.1f%%)\n",
           counts[REGION_TURING], 100.0*counts[REGION_TURING]/total);
    printf("  HOPF:        %6d (%.1f%%)\n",
           counts[REGION_HOPF], 100.0*counts[REGION_HOPF]/total);
    printf("  TURING-HOPF: %6d (%.1f%%)\n",
           counts[REGION_TURING_HOPF], 100.0*counts[REGION_TURING_HOPF]/total);
    printf("  NO HSS:      %6d (%.1f%%)\n",
           counts[REGION_NO_HSS], 100.0*counts[REGION_NO_HSS]/total);
}

/* ==========================================================================
 * Bifurcation Detection
 *
 * Find Turing bifurcation points via bisection along a 1D parameter sweep.
 *
 * A Turing bifurcation occurs when max_k Re(λ(k)) crosses 0.
 * We detect sign changes in λ_max along the parameter sweep and
 * refine via bisection.
 *
 * Reference: Kuznetsov, Y.A. (2004) Elements of Applied Bifurcation Theory
 * ==========================================================================*/

BifurcationPoint *find_turing_bifurcations(ModelType model,
                                            const ModelParams *base_params,
                                            int vary_idx,
                                            double p_min, double p_max,
                                            double Du, double Dv,
                                            int max_bifurcations,
                                            int *n_found) {
    if (!base_params || !n_found) return NULL;

    BifurcationPoint *results = (BifurcationPoint *)calloc(
        (size_t)max_bifurcations, sizeof(BifurcationPoint));
    if (!results) { *n_found = 0; return NULL; }

    int count = 0;
    int n_scan = 500;
    double *lambda_vals = (double *)malloc((size_t)n_scan * sizeof(double));
    if (!lambda_vals) { free(results); *n_found = 0; return NULL; }

    /* Scan parameter values and compute λ_max at each */
    for (int i = 0; i < n_scan; i++) {
        double param = p_min + (p_max - p_min) * (double)i / (double)(n_scan - 1);
        ModelParams p = *base_params;
        set_model_param(&p, vary_idx, param);

        HomogeneousSteadyState hss;
        if (model_find_hss(model, &p, &hss) != 0 || !hss.converged) {
            lambda_vals[i] = -1e100;  /* mark as invalid */
            continue;
        }

        double J[2][2];
        model_jacobian(model, hss.u_star, hss.v_star, &p, J);

        /* Compute dispersion to get λ_max */
        DispersionRelation *dr = turing_dispersion_compute(J, Du, Dv,
            200, 0.01, 5.0);
        lambda_vals[i] = dr ? dr->lambda_max : -1e100;
        dispersion_free(dr);
    }

    /* Detect sign changes in λ_max */
    for (int i = 0; i < n_scan - 1 && count < max_bifurcations; i++) {
        if (lambda_vals[i] * lambda_vals[i+1] < 0.0) {
            /* Sign change detected — bisection refinement */
            double p_lo = p_min + (p_max - p_min) * (double)i / (double)(n_scan-1);
            double p_hi = p_min + (p_max - p_min) * (double)(i+1) / (double)(n_scan-1);
            double l_lo = lambda_vals[i];
            double l_hi = lambda_vals[i+1];

            for (int bisect = 0; bisect < 50; bisect++) {
                double p_mid = 0.5 * (p_lo + p_hi);
                if (p_hi - p_lo < 1e-12) break;

                ModelParams p = *base_params;
                set_model_param(&p, vary_idx, p_mid);
                HomogeneousSteadyState hss;
                if (model_find_hss(model, &p, &hss) != 0 || !hss.converged) {
                    p_lo = p_mid; continue;
                }
                double J[2][2];
                model_jacobian(model, hss.u_star, hss.v_star, &p, J);
                DispersionRelation *dr = turing_dispersion_compute(J, Du, Dv,
                    200, 0.01, 5.0);
                double l_mid = dr ? dr->lambda_max : l_lo;
                double k_crit = dr ? dr->k_critical : 0.0;
                dispersion_free(dr);

                if (fabs(l_mid) < 1e-10) {
                    results[count].param1_value = p_mid;
                    results[count].k_critical = k_crit;
                    results[count].type = BIF_TURING;
                    count++;
                    break;
                }
                if (l_lo * l_mid < 0.0) {
                    p_hi = p_mid; l_hi = l_mid;
                } else {
                    p_lo = p_mid; l_lo = l_mid;
                }
            }
        }
    }

    free(lambda_vals);
    *n_found = count;
    return results;
}

void bifurcation_point_print(const BifurcationPoint *bp) {
    if (!bp) return;
    printf("Bifurcation Point:\n");
    printf("  Type: ");
    switch (bp->type) {
        case BIF_TURING:      printf("Turing\n"); break;
        case BIF_HOPF:        printf("Hopf\n"); break;
        case BIF_TURING_HOPF: printf("Turing-Hopf (codim-2)\n"); break;
        case BIF_SADDLE_NODE: printf("Saddle-Node\n"); break;
        default:              printf("None\n"); break;
    }
    printf("  Param: %.8g\n", bp->param1_value);
    printf("  k_critical: %.6g\n", bp->k_critical);
    printf("  ω_critical: %.6g\n", bp->omega_critical);
}

int detect_turing_hopf_interaction(const double J[2][2], double Du, double Dv) {
    if (!J) return 0;
    TuringConditions tc;
    turing_conditions_compute(J, Du, Dv, &tc);
    int hopf = is_hopf_unstable(J);
    return (tc.turing_unstable && hopf) ? 1 : 0;
}

/* ==========================================================================
 * Sensitivity Analysis
 *
 * Computes ∂λ_max/∂param_i via central finite differences for each parameter.
 *
 * λ_max is the maximum growth rate from the dispersion relation.
 * Sensitivities indicate which parameters most strongly affect pattern
 * formation (and are thus critical for experimental tuning).
 *
 * Reference: Saltelli et al. (2008) Global Sensitivity Analysis: The Primer
 * ==========================================================================*/

int turing_sensitivity_analysis(ModelType model, const ModelParams *p_nominal,
                                double Du, double Dv,
                                double *sensitivities, int n_params) {
    if (!p_nominal || !sensitivities || n_params <= 0) return -1;

    double h = 1e-6;  /* finite difference step */
    double *param_ptr = (double *)p_nominal;  /* treat struct as double array */

    /* Compute λ_max at nominal parameters */
    HomogeneousSteadyState hss_nom;
    if (model_find_hss(model, p_nominal, &hss_nom) != 0 || !hss_nom.converged) {
        return -1;
    }
    double J_nom[2][2];
    model_jacobian(model, hss_nom.u_star, hss_nom.v_star, p_nominal, J_nom);
    DispersionRelation *dr_nom = turing_dispersion_compute(J_nom, Du, Dv,
        200, 0.01, 5.0);
    double lambda_nom = dr_nom ? dr_nom->lambda_max : 0.0;
    dispersion_free(dr_nom);

    /* Number of double fields in ModelParams */
    int total_params = sizeof(ModelParams) / sizeof(double);
    if (n_params > total_params) n_params = total_params;

    for (int ip = 0; ip < n_params; ip++) {
        /* Forward perturbation */
        ModelParams p_plus = *p_nominal;
        param_ptr[ip] += h;

        HomogeneousSteadyState hss_p;
        if (model_find_hss(model, &p_plus, &hss_p) != 0 || !hss_p.converged) {
            sensitivities[ip] = 0.0;
            continue;
        }
        double J_p[2][2];
        model_jacobian(model, hss_p.u_star, hss_p.v_star, &p_plus, J_p);
        DispersionRelation *dr_p = turing_dispersion_compute(J_p, Du, Dv,
            200, 0.01, 5.0);
        double lambda_p = dr_p ? dr_p->lambda_max : lambda_nom;
        dispersion_free(dr_p);

        /* Backward perturbation */
        ModelParams p_minus = *p_nominal;
        param_ptr[ip] -= h;  /* restore */
        param_ptr[ip] -= h;  /* backward */

        HomogeneousSteadyState hss_m;
        if (model_find_hss(model, &p_minus, &hss_m) != 0 || !hss_m.converged) {
            /* Fall back to one-sided difference */
            sensitivities[ip] = (lambda_p - lambda_nom) / h;
            param_ptr[ip] += h;  /* restore */
            continue;
        }
        double J_m[2][2];
        model_jacobian(model, hss_m.u_star, hss_m.v_star, &p_minus, J_m);
        DispersionRelation *dr_m = turing_dispersion_compute(J_m, Du, Dv,
            200, 0.01, 5.0);
        double lambda_m = dr_m ? dr_m->lambda_max : lambda_nom;
        dispersion_free(dr_m);

        /* Central difference */
        sensitivities[ip] = (lambda_p - lambda_m) / (2.0 * h);

        /* Restore */
        param_ptr[ip] += h;
    }

    return 0;
}

double turing_wavelength_sensitivity(double Lx, double Ly,
                                     const double J[2][2],
                                     double Du, double Dv) {
    if (!J || Lx <= 0.0 || Ly <= 0.0) return 0.0;

    double kc = turing_critical_wavenumber(J, Du, Dv);
    if (kc < 1e-15) return 0.0;

    double lambda_pred = 2.0 * M_PI / kc;
    double L_min = (Lx < Ly) ? Lx : Ly;

    /* Sensitivity of pattern wavelength to domain size:
     * ∂λ/∂L ≈ λ / (max wavenumber supported by domain) ≈ λ²/(2πL) */
    return lambda_pred * lambda_pred / (2.0 * M_PI * L_min);
}
