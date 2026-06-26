/**
 * @file steady_state.c
 * @brief Implementation of steady-state analysis and MEP optimization.
 *
 * Implements the Minimum Entropy Production Principle (Prigogine 1945):
 * finding steady states that minimize sigma = X^T L X subject to
 * constraints, stability analysis, gradient flow dynamics, and
 * the Glansdorff-Prigogine universal evolution criterion.
 */

#include "steady_state.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* ---------------------------------------------------------------------------
 * L2: Core Concepts -- MEP steady state solvers
 * ------------------------------------------------------------------------- */

int ss_find_mep_state(const OnsagerMatrix *L,
                      const SteadyStateConstraints *constraints,
                      double forces[], int n, double tol)
{
    int n_free, n_fixed;
    int i, j, ii, jj;
    int free_idx[EP_MAX_FORCES];
    int fixed_idx[EP_MAX_FORCES];

    if (L == NULL || forces == NULL || n < 1)
        return -1;

    n_free = 0;
    n_fixed = 0;

    if (constraints != NULL) {
        for (i = 0; i < n; i++) {
            if (constraints->constrained[i]) {
                fixed_idx[n_fixed++] = i;
            } else {
                free_idx[n_free++] = i;
            }
        }
    } else {
        /* No constraints: all forces free. The MEP state is X = 0
         * (equilibrium) since sigma = X^T L X >= 0 with min at X = 0. */
        for (i = 0; i < n; i++) {
            forces[i] = 0.0;
        }
        return 0;
    }

    if (n_fixed == n) {
        /* All forces constrained: steady state is trivially the fixed values. */
        for (i = 0; i < n; i++) {
            forces[i] = constraints->fixed_values[i];
        }
        return 0;
    }

    if (n_free == 0) {
        return 0; /* No free forces */
    }

    /* Build reduced system for free forces:
     * L_uu * X_u = -L_uc * X_c
     * where X_u are free forces, X_c are constrained forces. */

    double L_uu[EP_MAX_FORCES * EP_MAX_FORCES];
    double b[EP_MAX_FORCES];
    double X_u[EP_MAX_FORCES];

    /* Build L_uu (n_free x n_free) */
    for (ii = 0; ii < n_free; ii++) {
        int i_real = free_idx[ii];
        for (jj = 0; jj < n_free; jj++) {
            int j_real = free_idx[jj];
            L_uu[ii * n_free + jj] = L->data[i_real * L->dim + j_real];
        }
    }

    /* Build right-hand side: b_i = -sum_j L_{i, j_fixed} * X_j_fixed */
    for (ii = 0; ii < n_free; ii++) {
        int i_real = free_idx[ii];
        b[ii] = 0.0;
        for (j = 0; j < n_fixed; j++) {
            int j_real = fixed_idx[j];
            b[ii] -= L->data[i_real * L->dim + j_real]
                   * constraints->fixed_values[j_real];
        }
    }

    /* Solve L_uu * X_u = b */
    if (onsager_solve(L_uu, b, X_u, n_free) != 0) {
        return -1;
    }

    /* Reconstruct full force vector */
    for (i = 0; i < n; i++) {
        forces[i] = 0.0;
    }
    for (ii = 0; ii < n_free; ii++) {
        forces[free_idx[ii]] = X_u[ii];
    }
    for (j = 0; j < n_fixed; j++) {
        forces[fixed_idx[j]] = constraints->fixed_values[fixed_idx[j]];
    }

    return 0;
}

int ss_compute_steady_fluxes(const OnsagerMatrix *L,
                             const double forces[], double fluxes[], int n)
{
    if (L == NULL || forces == NULL || fluxes == NULL) return -1;
    return onsager_compute_fluxes(L, forces, fluxes);
}

double ss_minimum_sigma(const OnsagerMatrix *L, const double forces[], int n)
{
    int i, j;
    double sigma = 0.0;

    if (L == NULL || forces == NULL || n < 1) return 0.0;

    for (i = 0; i < n; i++) {
        double row_sum = 0.0;
        for (j = 0; j < n; j++) {
            row_sum += L->data[i * L->dim + j] * forces[j];
        }
        sigma += forces[i] * row_sum;
    }

    return sigma;
}

/* ---------------------------------------------------------------------------
 * L3: Mathematical Structures -- Gradient flow dynamics
 * ------------------------------------------------------------------------- */

void ss_force_evolution(const double forces[], double dforces_dt[],
                        const OnsagerMatrix *L, double gamma,
                        const double fixed_forces[], int n)
{
    int i, j;

    if (forces == NULL || dforces_dt == NULL || L == NULL) return;

    for (i = 0; i < n; i++) {
        if (fixed_forces != NULL && !isnan(fixed_forces[i])) {
            /* This force is fixed by boundary conditions */
            dforces_dt[i] = 0.0;
            continue;
        }

        /* Gradient flow: dX_i/dt = -gamma * partial(sigma)/partial(X_i)
         * partial(sigma)/partial(X_i) = 2 sum_j L_{ij} X_j */
        double partial_sigma = 0.0;
        for (j = 0; j < n; j++) {
            partial_sigma += 2.0 * L->data[i * L->dim + j] * forces[j];
        }
        dforces_dt[i] = -gamma * partial_sigma;
    }
}

void ss_rk4_step(double forces[], const OnsagerMatrix *L, double gamma,
                 const double fixed_forces[], int n, double dt)
{
    double k1[EP_MAX_FORCES], k2[EP_MAX_FORCES];
    double k3[EP_MAX_FORCES], k4[EP_MAX_FORCES];
    double X_temp[EP_MAX_FORCES];
    int i;

    if (forces == NULL || L == NULL || n > EP_MAX_FORCES) return;

    /* RK4 for gradient flow:
     * X_{n+1} = X_n + (dt/6) * (k1 + 2*k2 + 2*k3 + k4) */

    /* k1 = f(X_n) */
    ss_force_evolution(forces, k1, L, gamma, fixed_forces, n);

    /* k2 = f(X_n + dt/2 * k1) */
    for (i = 0; i < n; i++) {
        X_temp[i] = forces[i] + 0.5 * dt * k1[i];
        /* Enforce fixed forces */
        if (fixed_forces != NULL && !isnan(fixed_forces[i])) {
            X_temp[i] = fixed_forces[i];
        }
    }
    ss_force_evolution(X_temp, k2, L, gamma, fixed_forces, n);

    /* k3 = f(X_n + dt/2 * k2) */
    for (i = 0; i < n; i++) {
        X_temp[i] = forces[i] + 0.5 * dt * k2[i];
        if (fixed_forces != NULL && !isnan(fixed_forces[i])) {
            X_temp[i] = fixed_forces[i];
        }
    }
    ss_force_evolution(X_temp, k3, L, gamma, fixed_forces, n);

    /* k4 = f(X_n + dt * k3) */
    for (i = 0; i < n; i++) {
        X_temp[i] = forces[i] + dt * k3[i];
        if (fixed_forces != NULL && !isnan(fixed_forces[i])) {
            X_temp[i] = fixed_forces[i];
        }
    }
    ss_force_evolution(X_temp, k4, L, gamma, fixed_forces, n);

    /* Update forces */
    for (i = 0; i < n; i++) {
        if (fixed_forces != NULL && !isnan(fixed_forces[i])) {
            /* Fixed forces remain unchanged */
            continue;
        }
        forces[i] += (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
}

int ss_simulate_relaxation(EvolutionTrajectory *trajectory,
                           const OnsagerMatrix *L, double gamma,
                           const double fixed_forces[],
                           const double forces_init[], int n,
                           double dt, int max_steps, double tol)
{
    double X[EP_MAX_FORCES];
    double dXdt[EP_MAX_FORCES];
    int step, i;

    if (trajectory == NULL || L == NULL || forces_init == NULL) return -1;
    if (n < 1 || n > EP_MAX_FORCES) return -1;

    /* Initialize forces */
    for (i = 0; i < n; i++) {
        X[i] = forces_init[i];
    }

    trajectory->is_monotonic = 1;
    double prev_sigma = ss_minimum_sigma(L, X, n);

    for (step = 0; step < max_steps; step++) {
        /* Record current state */
        if (step < trajectory->n_steps) {
            trajectory->time[step] = step * dt;
            for (i = 0; i < n; i++) {
                trajectory->forces_traj[step * n + i] = X[i];
            }
            trajectory->sigma_traj[step] = prev_sigma;

            /* Compute dsigma/dt */
            ss_force_evolution(X, dXdt, L, gamma, fixed_forces, n);
            double dsigma_dt = 0.0;
            for (i = 0; i < n; i++) {
                double partial_sigma = 0.0;
                int j;
                for (j = 0; j < n; j++) {
                    partial_sigma += 2.0 * L->data[i * L->dim + j] * X[j];
                }
                dsigma_dt += partial_sigma * dXdt[i];
            }
            trajectory->dsigma_dt_traj[step] = dsigma_dt;
        }

        /* Check convergence */
        if (fabs(trajectory->dsigma_dt_traj[step]) < tol) {
            trajectory->relaxation_time = step * dt;
            return step;
        }

        /* Integrate one step */
        ss_rk4_step(X, L, gamma, fixed_forces, n, dt);

        /* Compute new sigma and check monotonicity */
        double new_sigma = ss_minimum_sigma(L, X, n);
        if (new_sigma > prev_sigma + 1e-14) {
            trajectory->is_monotonic = 0;
        }
        prev_sigma = new_sigma;
    }

    trajectory->relaxation_time = max_steps * dt;
    return max_steps;
}

/* ---------------------------------------------------------------------------
 * L4: Fundamental Laws -- Stability theorems
 * ------------------------------------------------------------------------- */

int ss_linear_stability(StabilityResult *result, const double forces[],
                        const OnsagerMatrix *L, double gamma, int n)
{
    double J[EP_MAX_FORCES * EP_MAX_FORCES];
    double eigenvalues[EP_MAX_FORCES];
    int i, j, ret;

    if (result == NULL || L == NULL || n < 1 || n > EP_MAX_FORCES)
        return -1;

    /* Jacobian for gradient flow dX/dt = -2 gamma L X:
     * J_{ij} = -2 * gamma * L_{ij} */
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            J[i * n + j] = -2.0 * gamma * L->data[i * L->dim + j];
        }
    }

    /* Compute eigenvalues of Jacobian */
    double L_copy[EP_MAX_FORCES * EP_MAX_FORCES];
    double D[EP_MAX_FORCES];

    /* Use eigenvalues of L to compute eigenvalues of J:
     * lambda_i(J) = -2 * gamma * lambda_i(L) */
    OnsagerMatrix Ltemp;
    double Ldata[EP_MAX_FORCES * EP_MAX_FORCES];
    Ltemp.data = Ldata;
    Ltemp.dim = n;

    memcpy(Ldata, L->data, (size_t)(n * n) * sizeof(double));

    ret = onsager_eigenvalues(&Ltemp, eigenvalues, 100);
    if (ret != 0) return -1;

    result->dim = n;
    result->max_real_part = -2.0 * gamma * eigenvalues[n - 1]; /* smallest eigenvalue of L */
    result->is_asymptotically_stable = (result->max_real_part < 0.0);
    result->is_linearly_stable = (result->max_real_part <= 1e-14);
    result->spectral_abscissa = result->max_real_part;

    /* Damping ratio: min |Re(lambda)| / max |Im(lambda)|
     * For purely real eigenvalues, Im = 0, damping_ratio = infinity (stable node). */
    double min_abs_re = 1e100;
    double max_im = 0.0;
    for (i = 0; i < n; i++) {
        double re = -2.0 * gamma * eigenvalues[i];
        double abs_re = fabs(re);
        if (abs_re < min_abs_re && abs_re > 1e-14) min_abs_re = abs_re;
    }
    result->damping_ratio = (max_im > 1e-14) ? min_abs_re / max_im : 1e100;

    /* Store eigenvalues */
    for (i = 0; i < n; i++) {
        result->eigenvalues[i] = -2.0 * gamma * eigenvalues[i];
    }

    return 0;
}

double ss_excess_entropy_production(const double forces[],
                                    const double forces_steady[],
                                    const OnsagerMatrix *L, int n)
{
    double delta_X[EP_MAX_FORCES];
    double delta2_sigma = 0.0;
    int i, j;

    if (forces == NULL || forces_steady == NULL || L == NULL || n < 1)
        return 0.0;

    /* delta^2 sigma = sum_i (J_i - J_i^*) (X_i - X_i^*) = delta_X^T L delta_X */
    for (i = 0; i < n; i++) {
        if (forces_steady != NULL) {
            delta_X[i] = forces[i] - forces_steady[i];
        } else {
            delta_X[i] = forces[i];
        }
    }

    for (i = 0; i < n; i++) {
        double row_sum = 0.0;
        for (j = 0; j < n; j++) {
            row_sum += L->data[i * L->dim + j] * delta_X[j];
        }
        delta2_sigma += delta_X[i] * row_sum;
    }

    return delta2_sigma;
}

double ss_glansdorff_prigogine_criterion(const double forces[],
                                          const double dforces_dt[],
                                          const OnsagerMatrix *L, int n)
{
    double dX_sigma = 0.0;
    double fluxes[EP_MAX_FORCES];
    int i;

    if (forces == NULL || dforces_dt == NULL || L == NULL) return 0.0;

    /* Compute fluxes: J_i = sum_j L_{ij} X_j */
    onsager_compute_fluxes(L, forces, fluxes);

    /* d_X sigma = sum_i J_i * dX_i/dt */
    for (i = 0; i < n; i++) {
        dX_sigma += fluxes[i] * dforces_dt[i];
    }

    return dX_sigma;
}

/* ---------------------------------------------------------------------------
 * L5: Algorithms -- Constrained optimization
 * ------------------------------------------------------------------------- */

int ss_lagrange_mep(const double *L, const double *A, const double *b,
                    double *X, int n, int m)
{
    /* Build KKT system:
     * [ 2L   A^T ] [ X     ]   [ 0 ]
     * [ A    0   ] [ lambda ] = [ b ]
     *
     * Size: (n+m) x (n+m)
     * Solve using block elimination:
     *   1. Solve 2L * X_0 = 0 (which gives X_0 = 0 if L is PSD)
     *   2. Compute Schur complement: S = A * (2L)^{-1} * A^T
     *   3. Solve S * lambda = b - A * X_0
     *   4. Compute X = X_0 - (2L)^{-1} * A^T * lambda
     */

    double L_inv[64 * 64];  /* Inverse of 2L */
    double S[64 * 64];      /* Schur complement */
    double lambda[64];
    double temp[64];
    int i, j, k, ret;

    if (L == NULL || A == NULL || b == NULL || X == NULL) return -1;
    if (n < 1 || m < 1) return -1;
    if (n > 64 || m > 64) return -1;

    /* Step 1: Compute (2L)^{-1} */
    double L2[64 * 64];
    for (i = 0; i < n * n; i++) {
        L2[i] = 2.0 * L[i];
    }

    ret = onsager_invert(L2, L_inv, n);
    if (ret != 0) return -1;

    /* Step 2: Schur complement S = A * L_inv * A^T (m x m) */
    for (i = 0; i < m; i++) {
        for (j = 0; j < m; j++) {
            double sum = 0.0;
            for (k = 0; k < n; k++) {
                double A_Linv_ik = 0.0;
                int l;
                for (l = 0; l < n; l++) {
                    A_Linv_ik += A[i * n + l] * L_inv[l * n + k];
                }
                sum += A_Linv_ik * A[j * n + k];
            }
            S[i * m + j] = sum;
        }
    }

    /* Step 3: Solve S * lambda = b */
    ret = onsager_solve(S, b, lambda, m);
    if (ret != 0) return -1;

    /* Step 4: X = -(2L)^{-1} * A^T * lambda */
    for (i = 0; i < n; i++) {
        double AT_lambda_i = 0.0;
        for (j = 0; j < m; j++) {
            AT_lambda_i += A[j * n + i] * lambda[j];
        }
        temp[i] = AT_lambda_i;
    }

    for (i = 0; i < n; i++) {
        X[i] = 0.0;
        for (j = 0; j < n; j++) {
            X[i] -= L_inv[i * n + j] * temp[j];
        }
    }

    return 0;
}

int ss_projected_gradient(double forces[], const OnsagerMatrix *L,
                          const int constrained[], const double fixed_vals[],
                          int n, double step_size, int max_iter, double tol)
{
    double grad[EP_MAX_FORCES];
    int iter, i, j;

    if (forces == NULL || L == NULL || n < 1) return -1;

    for (iter = 0; iter < max_iter; iter++) {
        /* Compute gradient of unconstrained forces */
        double grad_norm = 0.0;
        for (i = 0; i < n; i++) {
            if (constrained != NULL && constrained[i]) {
                grad[i] = 0.0;
                continue;
            }
            grad[i] = 0.0;
            for (j = 0; j < n; j++) {
                grad[i] += 2.0 * L->data[i * L->dim + j] * forces[j];
            }
            grad_norm += grad[i] * grad[i];
        }
        grad_norm = sqrt(grad_norm);

        if (grad_norm < tol) {
            return iter;
        }

        /* Gradient step with projection */
        for (i = 0; i < n; i++) {
            if (constrained != NULL && constrained[i]) {
                if (fixed_vals != NULL) {
                    forces[i] = fixed_vals[i];
                }
                continue;
            }
            forces[i] -= step_size * grad[i];
        }
    }

    return -1;
}

/* ---------------------------------------------------------------------------
 * Utility functions
 * ------------------------------------------------------------------------- */

void ss_init_constraints(SteadyStateConstraints *c, int num_forces)
{
    int i;
    if (c == NULL) return;
    c->num_forces = num_forces;
    c->num_constrained = 0;
    for (i = 0; i < EP_MAX_FORCES; i++) {
        c->constrained[i] = 0;
        c->fixed_values[i] = 0.0;
    }
}

int ss_add_constraint(SteadyStateConstraints *c, int force_idx, double value)
{
    if (c == NULL) return -1;
    if (force_idx < 0 || force_idx >= c->num_forces) return -1;

    c->constrained[force_idx] = 1;
    c->fixed_values[force_idx] = value;
    c->num_constrained++;

    return 0;
}

void ss_init_steady_state(SteadyState *ss, int num_forces)
{
    int i;
    if (ss == NULL) return;
    memset(ss, 0, sizeof(SteadyState));
    ss->num_forces = num_forces;
    for (i = 0; i < EP_MAX_FORCES; i++) {
        ss->forces[i] = 0.0;
        ss->fluxes[i] = 0.0;
    }
}

void ss_free_steady_state(SteadyState *ss)
{
    /* Nothing heap-allocated in SteadyState (uses fixed-size arrays) */
    (void)ss;
}

void ss_init_trajectory(EvolutionTrajectory *traj, int dim, int max_steps)
{
    if (traj == NULL) return;

    traj->n_steps = max_steps;
    traj->dim = dim;
    traj->relaxation_time = 0.0;
    traj->is_monotonic = 1;

    traj->time = (double *)calloc((size_t)max_steps, sizeof(double));
    traj->forces_traj = (double *)calloc((size_t)(max_steps * dim), sizeof(double));
    traj->sigma_traj = (double *)calloc((size_t)max_steps, sizeof(double));
    traj->dsigma_dt_traj = (double *)calloc((size_t)max_steps, sizeof(double));
}

void ss_free_trajectory(EvolutionTrajectory *traj)
{
    if (traj == NULL) return;
    free(traj->time);
    free(traj->forces_traj);
    free(traj->sigma_traj);
    free(traj->dsigma_dt_traj);
}