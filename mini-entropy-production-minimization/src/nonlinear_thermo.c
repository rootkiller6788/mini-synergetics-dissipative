/**
 * @file nonlinear_thermo.c
 * @brief Implementation of nonlinear irreversible thermodynamics.
 *
 * Extends the linear Onsager formalism to nonlinear regimes:
 * nonlinear constitutive relations (polynomial flux-force),
 * generalized entropy production, nonlinear stability analysis,
 * local potential method (Glansdorff-Prigogine), and the
 * Maximum Entropy Production Principle debate.
 */

#include "nonlinear_thermo.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ---------------------------------------------------------------------------
 * L2: Core Concepts -- Nonlinear constitutive relations
 * ------------------------------------------------------------------------- */

int nt_compute_nonlinear_fluxes(const NonlinearOnsager *nl,
                                 const double forces[], double fluxes[])
{
    int i, j, k, l;
    int n;

    if (nl == NULL || forces == NULL || fluxes == NULL) return -1;

    n = nl->dim;

    /* Initialize fluxes to zero */
    for (i = 0; i < n; i++) {
        fluxes[i] = 0.0;
    }

    /* Linear term: J_i^(1) = sum_j L_{ij} X_j */
    if (nl->linear != NULL) {
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                fluxes[i] += nl->linear[i * n + j] * forces[j];
            }
        }
    }

    /* Quadratic term: J_i^(2) = sum_{j,k} M_{ijk} X_j X_k */
    if (nl->order >= 2 && nl->quadratic != NULL) {
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                for (k = 0; k < n; k++) {
                    double M_ijk = nl->quadratic[(i * n + j) * n + k];
                    fluxes[i] += M_ijk * forces[j] * forces[k];
                }
            }
        }
    }

    /* Cubic term: J_i^(3) = sum_{j,k,l} N_{ijkl} X_j X_k X_l */
    if (nl->order >= 3 && nl->cubic != NULL) {
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                for (k = 0; k < n; k++) {
                    for (l = 0; l < n; l++) {
                        double N_ijkl = nl->cubic[((i * n + j) * n + k) * n + l];
                        fluxes[i] += N_ijkl * forces[j] * forces[k] * forces[l];
                    }
                }
            }
        }
    }

    return 0;
}

int nt_nonlinear_sigma(const NonlinearOnsager *nl,
                       const double forces[], NonlinearEntropyProduction *ep)
{
    double fluxes[64];
    double sigma;
    int i;

    if (nl == NULL || forces == NULL || ep == NULL) return -1;
    if (nl->dim > 64) return -1;

    /* Compute fluxes from nonlinear relations */
    nt_compute_nonlinear_fluxes(nl, forces, fluxes);

    /* Total sigma = sum_i J_i X_i */
    sigma = 0.0;
    for (i = 0; i < nl->dim; i++) {
        sigma += fluxes[i] * forces[i];
    }

    ep->sigma_total = sigma;
    ep->is_positive = (sigma >= 0.0);

    /* Decompose by order */
    ep->sigma_linear = 0.0;
    ep->sigma_quadratic = 0.0;
    ep->sigma_cubic = 0.0;

    if (nl->linear != NULL) {
        for (i = 0; i < nl->dim; i++) {
            for (int j = 0; j < nl->dim; j++) {
                ep->sigma_linear += nl->linear[i * nl->dim + j]
                                  * forces[i] * forces[j];
            }
        }
    }

    ep->sigma_quadratic = sigma - ep->sigma_linear;
    ep->sigma_cubic = 0.0;  /* Already included in sigma_quadratic */

    return 0;
}

int nt_check_nonlinear_second_law(const NonlinearOnsager *nl,
                                   double X_max, int n_grid)
{
    int i, j;
    double grid[16][64];  /* Up to 16 grid points per dimension */

    if (nl == NULL || nl->dim < 1) return 0;
    if (n_grid > 16) n_grid = 16;

    /* Sample forces on a grid in the n-dimensional hypercube [-X_max, X_max]^n.
     * For simplicity, scan along each axis and diagonals. */

    double X[64];
    memset(X, 0, sizeof(X));

    /* Check along each axis */
    for (i = 0; i < nl->dim; i++) {
        for (j = 0; j < n_grid; j++) {
            double val = -X_max + 2.0 * X_max * j / (n_grid - 1);
            X[i] = val;

            NonlinearEntropyProduction ep;
            nt_nonlinear_sigma(nl, X, &ep);
            if (!ep.is_positive) return 0;

            X[i] = 0.0;
        }
    }

    /* Check diagonal directions */
    for (j = 0; j < n_grid; j++) {
        double val = -X_max + 2.0 * X_max * j / (n_grid - 1);
        for (i = 0; i < nl->dim; i++) {
            X[i] = val;
        }

        NonlinearEntropyProduction ep;
        nt_nonlinear_sigma(nl, X, &ep);
        if (!ep.is_positive) return 0;

        for (i = 0; i < nl->dim; i++) {
            X[i] = 0.0;
        }
    }

    return 1;
}

/* ---------------------------------------------------------------------------
 * L3: Mathematical Structures -- Nonlinear Jacobian and evolution
 * ------------------------------------------------------------------------- */

int nt_nonlinear_jacobian(const NonlinearOnsager *nl,
                          const double forces[], double *jacobian)
{
    int i, j, k;
    int n;

    if (nl == NULL || forces == NULL || jacobian == NULL) return -1;
    n = nl->dim;

    /* Jacobian of flux with respect to forces:
     * J_{ij} = partial(J_i) / partial(X_j)
     *        = L_{ij} + sum_k (M_{ijk} + M_{ikj}) X_k + O(X^2) */

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            double J_ij = 0.0;

            /* Linear part */
            if (nl->linear != NULL) {
                J_ij += nl->linear[i * n + j];
            }

            /* Quadratic contribution */
            if (nl->order >= 2 && nl->quadratic != NULL) {
                for (k = 0; k < n; k++) {
                    double M_ijk = nl->quadratic[(i * n + j) * n + k];
                    double M_ikj = nl->quadratic[(i * n + k) * n + j];
                    J_ij += (M_ijk + M_ikj) * forces[k];
                }
            }

            /* Cubic contribution */
            if (nl->order >= 3 && nl->cubic != NULL) {
                for (k = 0; k < n; k++) {
                    for (int l = 0; l < n; l++) {
                        double N_ijkl = nl->cubic[((i * n + j) * n + k) * n + l];
                        double N_ikjl = nl->cubic[((i * n + k) * n + j) * n + l];
                        double N_iklj = nl->cubic[((i * n + k) * n + l) * n + j];
                        J_ij += (N_ijkl + N_ikjl + N_iklj)
                              * forces[k] * forces[l];
                    }
                }
            }

            jacobian[i * n + j] = J_ij;
        }
    }

    return 0;
}

void nt_nonlinear_force_evolution(const double forces[], double dforces_dt[],
                                  const NonlinearOnsager *nl, double gamma,
                                  const double fixed[])
{
    int i, j, k;
    int n;
    double fluxes[64];

    if (forces == NULL || dforces_dt == NULL || nl == NULL) return;
    n = nl->dim;
    if (n > 64) return;

    nt_compute_nonlinear_fluxes(nl, forces, fluxes);

    /* Gradient flow with nonlinear sigma:
     * dX_i/dt = -gamma * partial(sigma)/partial(X_i)
     *
     * partial(sigma)/partial(X_i) = J_i + sum_k X_k * partial(J_k)/partial(X_i)
     *
     * For linear relations: partial(sigma)/partial(X_i) = 2 sum_j L_{ij} X_j
     * For nonlinear: need the derivative of the flux-force product.
     */

    for (i = 0; i < n; i++) {
        if (fixed != NULL && !isnan(fixed[i])) {
            dforces_dt[i] = 0.0;
            continue;
        }

        /* partial(sigma)/partial(X_i) = J_i + sum_k X_k * (partial J_k / partial X_i) */
        double partial_sigma = fluxes[i];

        /* Add contribution from derivative of fluxes */
        for (j = 0; j < n; j++) {
            double dJj_dXi = 0.0;

            /* Linear part */
            if (nl->linear != NULL) {
                dJj_dXi += nl->linear[j * n + i];
            }

            /* Quadratic part */
            if (nl->order >= 2 && nl->quadratic != NULL) {
                for (k = 0; k < n; k++) {
                    dJj_dXi += (nl->quadratic[(j * n + i) * n + k]
                              + nl->quadratic[(j * n + k) * n + i])
                              * forces[k];
                }
            }

            partial_sigma += forces[j] * dJj_dXi;
        }

        dforces_dt[i] = -gamma * partial_sigma;
    }
}

int nt_find_nonlinear_steady_states(double forces[], const NonlinearOnsager *nl,
                                    const double fixed[], int max_iter, double tol)
{
    double dXdt[64];
    int i, iter;
    int n;

    if (forces == NULL || nl == NULL) return -1;
    n = nl->dim;
    if (n > 64) return -1;

    /* Simple relaxation method: integrate gradient flow until steady state */
    double dt = 0.01;

    for (iter = 0; iter < max_iter; iter++) {
        nt_nonlinear_force_evolution(forces, dXdt, nl, 1.0, fixed);

        double max_dXdt = 0.0;
        for (i = 0; i < n; i++) {
            if (fixed != NULL && !isnan(fixed[i])) continue;
            if (fabs(dXdt[i]) > max_dXdt) max_dXdt = fabs(dXdt[i]);
        }

        if (max_dXdt < tol) {
            return 1;  /* Converged to a steady state */
        }

        /* Euler step */
        for (i = 0; i < n; i++) {
            if (fixed != NULL && !isnan(fixed[i])) continue;
            forces[i] += dt * dXdt[i];
        }
    }

    return -1;
}

/* ---------------------------------------------------------------------------
 * L4: Fundamental Laws -- Local potential and MEP/MaxEP debate
 * ------------------------------------------------------------------------- */

double nt_local_potential_evolution(const double forces[],
                                     const double dforces_dt[],
                                     const NonlinearOnsager *nl)
{
    /* d_X sigma = sum_i J_i * dX_i/dt (Glansdorff-Prigogine criterion)
     * This should be <= 0 for evolution toward a steady state. */

    double fluxes[64];
    double dX_sigma = 0.0;
    int i;

    if (forces == NULL || dforces_dt == NULL || nl == NULL) return 0.0;

    nt_compute_nonlinear_fluxes(nl, forces, fluxes);

    for (i = 0; i < nl->dim; i++) {
        dX_sigma += fluxes[i] * dforces_dt[i];
    }

    return dX_sigma;
}

int nt_find_max_entropy_state(double forces[], const NonlinearOnsager *nl,
                               const int constrained[], const double fixed_vals[],
                               double step_size, int max_iter, double tol)
{
    /* Gradient ASCENT (maximization instead of minimization):
     * X_new = X + alpha * partial(sigma)/partial(X)
     * with projection for constraints. */

    int i, j, k, iter;
    int n;
    double grad[64];
    double fluxes[64];

    if (forces == NULL || nl == NULL) return -1;
    n = nl->dim;
    if (n > 64) return -1;

    for (iter = 0; iter < max_iter; iter++) {
        nt_compute_nonlinear_fluxes(nl, forces, fluxes);

        /* Gradient of sigma with respect to forces:
         * partial(sigma)/partial(X_i) = J_i + sum_k X_k * partial(J_k)/partial(X_i) */
        double grad_norm = 0.0;
        for (i = 0; i < n; i++) {
            if (constrained != NULL && constrained[i]) {
                grad[i] = 0.0;
                continue;
            }

            grad[i] = fluxes[i];
            double jac[64 * 64];
            nt_nonlinear_jacobian(nl, forces, jac);
            for (j = 0; j < n; j++) {
                grad[i] += forces[j] * jac[j * n + i];
            }
            grad_norm += grad[i] * grad[i];
        }
        grad_norm = sqrt(grad_norm);

        if (grad_norm < tol) {
            return 0;  /* Converged to local extremum */
        }

        /* Gradient ascent step */
        for (i = 0; i < n; i++) {
            if (constrained != NULL && constrained[i]) {
                if (fixed_vals != NULL) forces[i] = fixed_vals[i];
                continue;
            }
            forces[i] += step_size * grad[i];  /* ASCENT (plus sign) */
        }
    }

    return -1;
}

/* ---------------------------------------------------------------------------
 * L5: Algorithms -- Data fitting and regime classification
 * ------------------------------------------------------------------------- */

int nt_fit_nonlinear_onsager(const double *forces, const double *fluxes,
                              int m, int n, int order, NonlinearOnsager *nl)
{
    /* Fit nonlinear constitutive relation J = f(X) using polynomial regression.
     * J_i = sum_j L_{ij} X_j + sum_{j,k} M_{ijk} X_j X_k + ...
     *
     * For simplicity, fit each output component J_i independently.
     */

    if (forces == NULL || fluxes == NULL || nl == NULL || m < 1 || n < 1)
        return -1;

    nt_init_nonlinear_onsager(nl, n, order);

    /* Number of parameters per output (linear + quadratic terms) */
    int n_linear = n;
    int n_quadratic = (order >= 2) ? n * n : 0;

    /* Least squares: solve A * coeff = b for each output i */
    for (int i_out = 0; i_out < n; i_out++) {
        int n_params = n_linear + n_quadratic;
        double *A = (double *)calloc((size_t)(n_params * n_params), sizeof(double));
        double *b = (double *)calloc((size_t)n_params, sizeof(double));

        if (A == NULL || b == NULL) {
            free(A); free(b);
            return -1;
        }

        /* Build normal equations: A = Phi^T Phi, b = Phi^T y */
        for (int sample = 0; sample < m; sample++) {
            double phi[128];  /* Design matrix row */
            int p = 0;

            /* Linear features: X_0, X_1, ..., X_{n-1} */
            for (int j = 0; j < n; j++) {
                phi[p++] = forces[sample * n + j];
            }

            /* Quadratic features: X_0 X_0, X_0 X_1, ..., X_{n-1} X_{n-1} */
            if (order >= 2) {
                for (int j1 = 0; j1 < n; j1++) {
                    for (int j2 = 0; j2 < n; j2++) {
                        phi[p++] = forces[sample * n + j1]
                                 * forces[sample * n + j2];
                    }
                }
            }

            /* Accumulate A = Phi^T Phi */
            for (int p1 = 0; p1 < n_params; p1++) {
                for (int p2 = 0; p2 < n_params; p2++) {
                    A[p1 * n_params + p2] += phi[p1] * phi[p2];
                }
                b[p1] += phi[p1] * fluxes[sample * n + i_out];
            }
        }

        /* Solve A * coeff = b with simple Gaussian elimination */
        /* (For production use, LDL^T would be better) */
        double *coeff = (double *)calloc((size_t)n_params, sizeof(double));

        /* Simple approach: if n is small, use direct inversion */
        if (n <= 4 && n_params <= 32) {
            /* Copy A to augmented matrix */
            double *aug = (double *)malloc((size_t)(n_params * (n_params + 1)) * sizeof(double));
            for (int r = 0; r < n_params; r++) {
                for (int c = 0; c < n_params; c++) {
                    aug[r * (n_params + 1) + c] = A[r * n_params + c];
                }
                aug[r * (n_params + 1) + n_params] = b[r];
            }

            /* Gaussian elimination with partial pivoting */
            for (int col = 0; col < n_params; col++) {
                double max_val = fabs(aug[col * (n_params + 1) + col]);
                int max_row = col;
                for (int row = col + 1; row < n_params; row++) {
                    if (fabs(aug[row * (n_params + 1) + col]) > max_val) {
                        max_val = fabs(aug[row * (n_params + 1) + col]);
                        max_row = row;
                    }
                }
                if (max_val < 1e-15) continue;

                if (max_row != col) {
                    for (int c = 0; c < n_params + 1; c++) {
                        double tmp = aug[col * (n_params + 1) + c];
                        aug[col * (n_params + 1) + c] = aug[max_row * (n_params + 1) + c];
                        aug[max_row * (n_params + 1) + c] = tmp;
                    }
                }

                for (int row = col + 1; row < n_params; row++) {
                    double factor = aug[row * (n_params + 1) + col]
                                  / aug[col * (n_params + 1) + col];
                    for (int c = col; c < n_params + 1; c++) {
                        aug[row * (n_params + 1) + c] -= factor
                            * aug[col * (n_params + 1) + c];
                    }
                }
            }

            /* Back substitution */
            for (int row = n_params - 1; row >= 0; row--) {
                double sum = aug[row * (n_params + 1) + n_params];
                for (int c = row + 1; c < n_params; c++) {
                    sum -= aug[row * (n_params + 1) + c] * coeff[c];
                }
                coeff[row] = (fabs(aug[row * (n_params + 1) + row]) > 1e-15)
                           ? sum / aug[row * (n_params + 1) + row] : 0.0;
            }
            free(aug);
        }

        /* Store fitted coefficients */
        int p = 0;
        for (int j = 0; j < n; j++) {
            nl->linear[i_out * n + j] = coeff[p++];
        }
        if (order >= 2) {
            for (int j1 = 0; j1 < n; j1++) {
                for (int j2 = 0; j2 < n; j2++) {
                    nl->quadratic[(i_out * n + j1) * n + j2] = coeff[p++];
                }
            }
        }

        free(A); free(b); free(coeff);
    }

    return 0;
}

void nt_classify_regime(FarFromEquilibriumState *state,
                        const double forces[], const OnsagerMatrix *L)
{
    double norm_X = 0.0;
    int i, n;

    if (state == NULL || forces == NULL || L == NULL) return;

    n = L->dim;

    /* Force magnitude */
    for (i = 0; i < n; i++) {
        norm_X += forces[i] * forces[i];
    }
    norm_X = sqrt(norm_X);
    state->distance_from_equilibrium = norm_X;

    /* Linearity measure: compute fluxes, compare to linear approximation */
    double fluxes_linear[64];
    onsager_compute_fluxes(L, forces, fluxes_linear);

    double norm_J_linear = 0.0;
    for (i = 0; i < n; i++) {
        norm_J_linear += fluxes_linear[i] * fluxes_linear[i];
    }
    norm_J_linear = sqrt(norm_J_linear);

    /* Without actual nonlinear flux data, we estimate linearity
     * based on the force magnitude: for |X| small, linear is good. */
    if (norm_X < 1e-6) {
        state->regime = 0; /* Equilibrium */
        state->linearity_measure = 0.0;
    } else if (norm_X < 0.1) {
        state->regime = 1; /* Linear regime */
        state->linearity_measure = norm_X / 0.1;
    } else if (norm_X < 1.0) {
        state->regime = 2; /* Nonlinear near-equilibrium */
        state->linearity_measure = 1.0;
    } else {
        state->regime = 3; /* Far from equilibrium */
        state->linearity_measure = 2.0;
    }
}

/* ---------------------------------------------------------------------------
 * Utility functions
 * ------------------------------------------------------------------------- */

void nt_init_nonlinear_onsager(NonlinearOnsager *nl, int dim, int order)
{
    if (nl == NULL) return;

    memset(nl, 0, sizeof(NonlinearOnsager));
    nl->dim = dim;
    nl->order = order;

    /* Allocate linear matrix */
    nl->linear = (double *)calloc((size_t)(dim * dim), sizeof(double));

    /* Allocate quadratic tensor */
    if (order >= 2) {
        nl->quadratic = (double *)calloc((size_t)(dim * dim * dim), sizeof(double));
    }

    /* Allocate cubic tensor */
    if (order >= 3) {
        nl->cubic = (double *)calloc((size_t)(dim * dim * dim * dim), sizeof(double));
    }
}

void nt_free_nonlinear_onsager(NonlinearOnsager *nl)
{
    if (nl == NULL) return;
    free(nl->linear);
    free(nl->quadratic);
    free(nl->cubic);
    memset(nl, 0, sizeof(NonlinearOnsager));
}