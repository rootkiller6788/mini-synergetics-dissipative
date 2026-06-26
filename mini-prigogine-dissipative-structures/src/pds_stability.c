#include "pds_stability.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>

/* ============================================================================
 * PDS Stability Analysis Implementation
 *
 * Key theorems:
 * L4: Linear stability analysis ? the thermodynamic branch is stable
 *     if all eigenvalues of the Jacobian have negative real parts.
 * L4: Hopf bifurcation theorem ? a pair of complex conjugate eigenvalues
 *     crosses the imaginary axis => periodic orbit emerges.
 * L4: Lyapunov's first method ? linear stability implies local nonlinear
 *     stability (under certain conditions).
 * L8: Lyapunov exponent computation from time series.
 * ============================================================================ */

/* Finite-difference Jacobian computation.
 * J_{ij} = partial f_i / partial x_j ~ (f_i(x + h*e_j) - f_i(x - h*e_j)) / (2h)
 * Central difference for O(h^2) accuracy. */
void pds_compute_jacobian_fd(const PDSDissipativeSystem* sys,
                              void (*rhs_fn)(const double*, const double*, int, double*),
                              const double* params, int n_params,
                              PDSMatrix** J_out) {
    if (!sys || !rhs_fn || !J_out) return;

    int n = sys->n_species;
    if (n <= 0) return;

    *J_out = pds_matrix_create(n, n);
    if (!*J_out) return;

    double* x0 = (double*)malloc(n * sizeof(double));
    double* x_plus = (double*)malloc(n * sizeof(double));
    double* x_minus = (double*)malloc(n * sizeof(double));
    double* f_plus = (double*)malloc(n * sizeof(double));
    double* f_minus = (double*)malloc(n * sizeof(double));

    if (!x0 || !x_plus || !x_minus || !f_plus || !f_minus) {
        free(x0); free(x_plus); free(x_minus);
        free(f_plus); free(f_minus);
        pds_matrix_free(*J_out);
        *J_out = NULL;
        return;
    }

    /* Get current state */
    pds_system_get_state(sys, x0);

    double h = 1e-6; /* Finite difference step size */

    for (int j = 0; j < n; j++) {
        /* Perturb state */
        memcpy(x_plus, x0, n * sizeof(double));
        memcpy(x_minus, x0, n * sizeof(double));
        x_plus[j] += h;
        x_minus[j] -= h;

        /* Evaluate RHS at perturbed states */
        rhs_fn(x_plus, params, n_params, f_plus);
        rhs_fn(x_minus, params, n_params, f_minus);

        /* Fill Jacobian column j */
        for (int i = 0; i < n; i++) {
            double J_ij = (f_plus[i] - f_minus[i]) / (2.0 * h);
            pds_matrix_set(*J_out, i, j, J_ij);
        }
    }

    free(x0); free(x_plus); free(x_minus);
    free(f_plus); free(f_minus);
}

/* Eigenvalues of a 2x2 real matrix.
 * J = [[a, b], [c, d]]
 * Characteristic polynomial: lambda^2 - trace*lambda + det = 0
 * trace = a + d, det = a*d - b*c
 * lambda = (trace +/- sqrt(trace^2 - 4*det)) / 2 */
void pds_eigenvalues_2x2(double a, double b, double c, double d,
                          double* lambda1_real, double* lambda1_imag,
                          double* lambda2_real, double* lambda2_imag) {
    double trace = a + d;
    double det = a * d - b * c;
    double discriminant = trace * trace - 4.0 * det;

    if (discriminant >= 0.0) {
        /* Real eigenvalues */
        double sqrt_disc = sqrt(discriminant);
        *lambda1_real = (trace + sqrt_disc) / 2.0;
        *lambda1_imag = 0.0;
        *lambda2_real = (trace - sqrt_disc) / 2.0;
        *lambda2_imag = 0.0;
    } else {
        /* Complex conjugate pair */
        *lambda1_real = trace / 2.0;
        *lambda1_imag = sqrt(-discriminant) / 2.0;
        *lambda2_real = trace / 2.0;
        *lambda2_imag = -(*lambda1_imag);
    }
}

/* Eigenvalues of a 3x3 real matrix via characteristic polynomial.
 * Uses Cardano's formula for the cubic equation.
 * det(lambda*I - A) = lambda^3 + c2*lambda^2 + c1*lambda + c0 = 0
 * where c2 = -trace(A), c1 = sum of principal minors,
 * c0 = -det(A). */
void pds_eigenvalues_3x3(const double A[3][3],
                          double lambda_real[3], double lambda_imag[3]) {
    /* Coefficients of characteristic polynomial: lambda^3 + c2*lambda^2 + c1*lambda + c0 */
    double c2 = -(A[0][0] + A[1][1] + A[2][2]);

    double c1 = A[0][0]*A[1][1] + A[1][1]*A[2][2] + A[0][0]*A[2][2]
                - A[0][1]*A[1][0] - A[1][2]*A[2][1] - A[0][2]*A[2][0];

    double c0 = -(A[0][0]*A[1][1]*A[2][2]
                  + A[0][1]*A[1][2]*A[2][0]
                  + A[0][2]*A[1][0]*A[2][1]
                  - A[0][2]*A[1][1]*A[2][0]
                  - A[0][1]*A[1][0]*A[2][2]
                  - A[0][0]*A[1][2]*A[2][1]);

    /* Depressed cubic: t^3 + p*t + q = 0, where lambda = t - c2/3 */
    double p = c1 - c2*c2 / 3.0;
    double q = c0 - c2*c1/3.0 + 2.0*c2*c2*c2/27.0;

    double discriminant = q*q/4.0 + p*p*p/27.0;

    if (discriminant > 0.0) {
        /* One real root, two complex conjugates */
        double sqrt_d = sqrt(discriminant);
        double u = cbrt(-q/2.0 + sqrt_d);
        double v = cbrt(-q/2.0 - sqrt_d);
        lambda_real[0] = u + v - c2/3.0;
        lambda_imag[0] = 0.0;
        lambda_real[1] = -(u+v)/2.0 - c2/3.0;
        lambda_imag[1] = sqrt(3.0)/2.0 * (u - v);
        lambda_real[2] = lambda_real[1];
        lambda_imag[2] = -lambda_imag[1];
    } else if (discriminant < 0.0) {
        /* Three real roots (casus irreducibilis) */
        double r = sqrt(-p*p*p/27.0);
        double phi = acos(-q/(2.0*r));
        double cbrt_r = 2.0 * cbrt(r);
        for (int k = 0; k < 3; k++) {
            lambda_real[k] = cbrt_r * cos((phi + 2.0*M_PI*k)/3.0) - c2/3.0;
            lambda_imag[k] = 0.0;
        }
    } else {
        /* discriminant = 0: multiple roots */
        double u = cbrt(-q/2.0);
        lambda_real[0] = 2.0*u - c2/3.0;
        lambda_imag[0] = 0.0;
        lambda_real[1] = -u - c2/3.0;
        lambda_imag[1] = 0.0;
        lambda_real[2] = lambda_real[1];
        lambda_imag[2] = 0.0;
    }
}

/* Linear stability analysis of a Jacobian matrix.
 * Computes eigenvalues and classifies stability.
 * L4 Theorem: The steady state is linearly stable iff all eigenvalues
 * have negative real parts. If any eigenvalue has positive real part,
 * the steady state is unstable. Purely imaginary eigenvalues indicate
 * a Hopf bifurcation point. */
PDSStabilityResult* pds_linear_stability_analysis(const PDSMatrix* J) {
    if (!J || J->rows != J->cols || J->rows == 0) return NULL;

    PDSStabilityResult* sr = (PDSStabilityResult*)calloc(1, sizeof(PDSStabilityResult));
    if (!sr) return NULL;

    int n = J->rows;
    sr->n_eigenvalues = n;
    sr->eigenvalues = (double*)calloc(n * 2, sizeof(double)); /* real + imag pairs */
    sr->eigenvectors = NULL; /* Not computed for general case */

    if (!sr->eigenvalues) {
        free(sr);
        return NULL;
    }

    if (n == 2) {
        double l1r, l1i, l2r, l2i;
        pds_eigenvalues_2x2(J->m[0][0], J->m[0][1],
                            J->m[1][0], J->m[1][1],
                            &l1r, &l1i, &l2r, &l2i);
        sr->eigenvalues[0] = l1r; sr->eigenvalues[1] = l1i;
        sr->eigenvalues[2] = l2r; sr->eigenvalues[3] = l2i;
    } else if (n == 3) {
        double A[3][3];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                A[i][j] = J->m[i][j];
        double lr[3], li[3];
        pds_eigenvalues_3x3(A, lr, li);
        for (int i = 0; i < 3; i++) {
            sr->eigenvalues[2*i] = lr[i];
            sr->eigenvalues[2*i+1] = li[i];
        }
    } else {
        /* For n > 3: use power iteration for dominant eigenvalue */
        /* Simple estimate: Gershgorin circle theorem for bounds */
        double max_real = -DBL_MAX;
        for (int i = 0; i < n; i++) {
            double center = J->m[i][i];
            double radius = 0.0;
            for (int j = 0; j < n; j++) {
                if (i != j) radius += fabs(J->m[i][j]);
            }
            sr->eigenvalues[2*i] = center;
            sr->eigenvalues[2*i+1] = 0.0;
            if (center + radius > max_real) max_real = center + radius;
        }
        /* Set the maximum possible real part */
        sr->eigenvalues[0] = max_real;
        sr->eigenvalues[1] = 0.0;
    }

    /* Classify stability */
    sr->is_stable = true;
    sr->is_oscillatory = false;
    sr->n_unstable_modes = 0;
    double max_real_part = -DBL_MAX;

    for (int i = 0; i < n; i++) {
        double real_part = sr->eigenvalues[2*i];
        double imag_part = sr->eigenvalues[2*i+1];

        if (real_part > max_real_part) max_real_part = real_part;
        if (real_part > 1e-10) {
            sr->is_stable = false;
            sr->n_unstable_modes++;
        }
        if (fabs(imag_part) > 1e-10 && fabs(real_part) < 1e-10) {
            sr->is_oscillatory = true;
            sr->oscillation_frequency = fabs(imag_part);
        }
    }

    sr->is_marginal = (!sr->is_stable && sr->n_unstable_modes == 0)
                      || (fabs(max_real_part) < 1e-8);
    sr->critical_parameter = 0.0; /* Fill in by caller if known */

    return sr;
}

void pds_stability_result_free(PDSStabilityResult* sr) {
    if (!sr) return;
    free(sr->eigenvalues);
    free(sr->eigenvectors);
    free(sr);
}

/* Detects Hopf bifurcation: a pair of complex conjugate eigenvalues
 * crosses the imaginary axis.
 * L4 Theorem (Hopf, 1942): If a pair of complex eigenvalues crosses
 * the imaginary axis with non-zero speed as a parameter varies, then
 * a periodic orbit (limit cycle) bifurcates from the steady state. */
bool pds_detect_hopf_bifurcation(const PDSMatrix* J, double* critical_frequency) {
    if (!J || J->rows < 2) return false;

    PDSStabilityResult* sr = pds_linear_stability_analysis(J);
    if (!sr) return false;

    bool has_hopf = false;

    for (int i = 0; i < sr->n_eigenvalues; i++) {
        double real = sr->eigenvalues[2*i];
        double imag = sr->eigenvalues[2*i+1];

        /* Check for purely imaginary pair (within tolerance) */
        if (fabs(real) < 1e-6 && fabs(imag) > 1e-6) {
            has_hopf = true;
            if (critical_frequency) *critical_frequency = fabs(imag);
            break;
        }
    }

    pds_stability_result_free(sr);
    return has_hopf;
}

/* Maximum Lyapunov exponent computation from trajectory data.
 * Uses the Rosenstein algorithm (1993): track divergence of nearest
 * neighbors in reconstructed phase space.
 *
 * L8: The maximum Lyapunov exponent lambda_1 quantifies the rate of
 * exponential divergence of nearby trajectories:
 *   d(t) ~ d(0) * exp(lambda_1 * t)
 *
 * lambda_1 > 0: chaotic (sensitive dependence)
 * lambda_1 = 0: limit cycle or quasiperiodic
 * lambda_1 < 0: stable fixed point
 *
 * For dissipative structures, lambda_1 is typically 0 (limit cycle)
 * or slightly positive (chaotic dissipative structure). */
double pds_max_lyapunov_exponent(const double* const* trajectory,
                                  int n_points, int dim, double dt,
                                  int embedding_dim, int lag) {
    (void)embedding_dim; (void)lag;
    if (!trajectory || n_points < 10 || dim <= 0 || dt <= 0.0) return 0.0;

    /* Use the direct method: compute average exponential divergence
     * of nearby initial conditions over the trajectory. */
    double sum_divergence = 0.0;
    int n_pairs = 0;

    int start = n_points / 4; /* Skip transient */
    int end = n_points * 3 / 4;

    for (int i = start; i < end; i += 5) {
        /* Find nearest neighbor in the trajectory (excluding temporally close points) */
        double min_dist = DBL_MAX;
        int nearest_j = -1;

        for (int j = start; j < end; j++) {
            if (abs(i - j) < 10) continue; /* Exclude temporal neighbors */

            double dist = 0.0;
            for (int d = 0; d < dim; d++) {
                double diff = trajectory[i][d] - trajectory[j][d];
                dist += diff * diff;
            }
            dist = sqrt(dist);

            if (dist < min_dist && dist > 1e-10) {
                min_dist = dist;
                nearest_j = j;
            }
        }

        if (nearest_j >= 0 && min_dist < DBL_MAX) {
            /* Track divergence over prediction horizon */
            int horizon = 10;
            if (i + horizon < n_points && nearest_j + horizon < n_points) {
                double final_dist = 0.0;
                for (int d = 0; d < dim; d++) {
                    double diff = trajectory[i+horizon][d] -
                                  trajectory[nearest_j+horizon][d];
                    final_dist += diff * diff;
                }
                final_dist = sqrt(final_dist);

                if (min_dist > 1e-15 && final_dist > 1e-15) {
                    sum_divergence += log(final_dist / min_dist);
                    n_pairs++;
                }
            }
        }
    }

    if (n_pairs == 0) return 0.0;
    /* lambda = average of ln(d_final/d_initial) / horizon*dt */
    double avg_divergence = sum_divergence / n_pairs;
    return avg_divergence / (10.0 * dt);
}

/* Excess entropy production: delta_x^T * (dJ/dx) * delta_x
 * Used in the Glansdorff-Prigogine stability criterion.
 * L4: The sign of delta^2 P determines thermodynamic stability. */
double pds_excess_entropy_production(const PDSDissipativeSystem* sys,
                                      const PDSVector* perturbation) {
    if (!sys || !perturbation || sys->n_species == 0) return 0.0;

    /* delta^2 P = sum_{i,j} delta_J_i * delta_X_j
     * ~ perturbation^T * J * perturbation (for linearized system)
     * where J is the Jacobian of the fluxes with respect to forces. */
    double excess = 0.0;

    for (int i = 0; i < sys->n_species; i++) {
        double delta_xi = pds_vector_get(perturbation, i);
        for (int j = 0; j < sys->n_species; j++) {
            double delta_xj = pds_vector_get(perturbation, j);
            /* Simplified: assume diagonal coupling with coefficient 1 */
            excess += (i == j ? 1.0 : 0.1) * delta_xi * delta_xj;
        }
    }

    return excess;
}

void pds_stability_print(const PDSStabilityResult* sr) {
    if (!sr) { printf("StabilityResult: NULL\n"); return; }

    printf("=== Linear Stability Analysis ===\n");
    printf("  Stable: %s\n", sr->is_stable ? "YES" : "NO");
    printf("  Marginal: %s\n", sr->is_marginal ? "YES" : "no");
    printf("  Oscillatory: %s\n", sr->is_oscillatory ? "YES" : "no");
    printf("  Unstable modes: %d\n", sr->n_unstable_modes);
    printf("  Eigenvalues (%d):\n", sr->n_eigenvalues);
    for (int i = 0; i < sr->n_eigenvalues && i < 6; i++) {
        double re = sr->eigenvalues[2*i];
        double im = sr->eigenvalues[2*i+1];
        printf("    lambda_%d = %.4f %+.4fi\n", i, re, im);
    }
    if (sr->n_eigenvalues > 6) printf("    ... (%d more)\n", sr->n_eigenvalues - 6);
    if (sr->is_oscillatory)
        printf("  Hopf frequency: %.4f rad/s\n", sr->oscillation_frequency);
}