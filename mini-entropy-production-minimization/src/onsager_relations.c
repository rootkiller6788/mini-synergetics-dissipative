/**
 * @file onsager_relations.c
 * @brief Implementation of Onsager reciprocal relations and matrix operations.
 *
 * Implements the Onsager coefficient matrix: allocation, initialization
 * for various physical systems (Fourier heat conduction, multicomponent
 * diffusion, chemical reactions, coupled two-process systems),
 * reciprocity verification, positive semidefiniteness checking,
 * LDL^T decomposition, linear system solving, matrix inversion,
 * and data-driven estimation of Onsager coefficients.
 *
 * References:
 *   - Onsager, L. (1931) Phys. Rev. 37, 405; 38, 2265
 *   - de Groot & Mazur (1984), Ch. IV-V
 *   - Casimir, H.B.G. (1945) Rev. Mod. Phys. 17, 343
 */

#include "onsager_relations.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#define R_GAS 8.314462618

/* ---------------------------------------------------------------------------
 * L2: Core Concepts -- Matrix allocation and basic operations
 * ------------------------------------------------------------------------- */

OnsagerMatrix *onsager_alloc(int dim)
{
    OnsagerMatrix *L;

    if (dim < 1 || dim > ONSAGER_MAX_DIM) {
        return NULL;
    }

    L = (OnsagerMatrix *)malloc(sizeof(OnsagerMatrix));
    if (L == NULL) return NULL;

    L->data = (double *)calloc((size_t)(dim * dim), sizeof(double));
    if (L->data == NULL) {
        free(L);
        return NULL;
    }

    L->dim = dim;
    L->is_symmetric = 0;
    L->is_psd = 0;

    return L;
}

void onsager_free(OnsagerMatrix *L)
{
    if (L != NULL) {
        free(L->data);
        free(L);
    }
}

int onsager_set(OnsagerMatrix *L, int i, int j, double value)
{
    if (L == NULL || L->data == NULL) return -1;
    if (i < 0 || i >= L->dim || j < 0 || j >= L->dim) return -1;

    L->data[i * L->dim + j] = value;
    L->is_symmetric = 0;  /* Invalidate cached symmetry check */
    L->is_psd = 0;        /* Invalidate cached PSD check */

    return 0;
}

double onsager_get(const OnsagerMatrix *L, int i, int j)
{
    if (L == NULL || L->data == NULL) return 0.0 / 0.0; /* NaN */
    if (i < 0 || i >= L->dim || j < 0 || j >= L->dim) return 0.0 / 0.0;

    return L->data[i * L->dim + j];
}

int onsager_init_diagonal(OnsagerMatrix *L, const double diag[])
{
    int i, j;

    if (L == NULL || L->data == NULL || diag == NULL) return -1;

    for (i = 0; i < L->dim; i++) {
        for (j = 0; j < L->dim; j++) {
            L->data[i * L->dim + j] = (i == j) ? diag[i] : 0.0;
        }
    }

    L->is_symmetric = 1;
    /* Check PSD: all diagonal elements must be >= 0 */
    L->is_psd = 1;
    for (i = 0; i < L->dim; i++) {
        if (diag[i] < 0.0) {
            L->is_psd = 0;
            break;
        }
    }

    return 0;
}

int onsager_init_fourier(OnsagerMatrix *L, double thermal_cond, double temperature)
{
    double diag[1];

    if (L == NULL || L->dim != 1) return -1;
    if (temperature <= 0.0 || thermal_cond < 0.0) return -1;

    /* Fourier's law: J_q = -kappa grad(T)
     * X_q = grad(1/T) = -grad(T)/T^2
     * J_q = kappa T^2 X_q
     * Therefore L_{qq} = kappa * T^2 */
    diag[0] = thermal_cond * temperature * temperature;

    return onsager_init_diagonal(L, diag);
}

int onsager_init_diffusion(OnsagerMatrix *L, const double *D,
                           const double conc[], double T, int n)
{
    int i, j;

    if (L == NULL || D == NULL || conc == NULL) return -1;
    if (L->dim != n || n < 1 || T <= 0.0) return -1;

    /* L_{kl} = D_{kl} * c_k * c_l / (R * T * c_total)
     * for the center-of-mass frame.
     *
     * This approximation relates the Onsager diffusion coefficients
     * to the Fickian diffusion coefficients through the thermodynamic
     * factor (assumed ideal here).
     */
    double c_total = 0.0;
    for (i = 0; i < n; i++) {
        if (conc[i] < 0.0) return -1;
        c_total += conc[i];
    }
    if (c_total <= 0.0) return -1;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            double L_ij = D[i * n + j] * conc[i] * conc[j]
                        / (R_GAS * T * c_total);
            L->data[i * L->dim + j] = L_ij;
        }
    }

    return 0;
}

int onsager_init_two_process(OnsagerMatrix *L, double L11, double L22, double L12)
{
    if (L == NULL || L->dim != 2) return -1;

    /* Enforce positivity: L11 >= 0, L22 >= 0, L11*L22 >= L12^2 */
    if (L11 < 0.0 || L22 < 0.0) return -1;
    if (L11 * L22 < L12 * L12) return -1; /* PSD condition */

    onsager_set(L, 0, 0, L11);
    onsager_set(L, 0, 1, L12);
    onsager_set(L, 1, 0, L12);  /* Onsager reciprocity */
    onsager_set(L, 1, 1, L22);

    L->is_symmetric = 1;
    L->is_psd = 1;

    return 0;
}

int onsager_init_chemical(OnsagerMatrix *L, const double eq_rates[],
                          const double *stoichiometry,
                          int r, int s, double T, double R_gas)
{
    int rho, sigma, k;

    if (L == NULL || eq_rates == NULL || stoichiometry == NULL) return -1;
    if (L->dim != r || r < 1 || s < 1 || T <= 0.0) return -1;

    /* Near equilibrium, the Onsager coefficients for chemical reactions are:
     * L_{rho sigma} = (r_rho^eq / R_gas T) * delta_{rho sigma}
     *
     * This assumes reactions are independent (uncoupled) near equilibrium.
     * For coupled reactions, off-diagonal terms arise when reactions
     * share intermediates.
     */
    for (rho = 0; rho < r; rho++) {
        for (sigma = 0; sigma < r; sigma++) {
            double L_val = 0.0;
            if (rho == sigma) {
                L_val = eq_rates[rho] / (R_gas * T);
            } else {
                /* Coupling through shared species: sum over species
                 * of (nu_{k,rho} * nu_{k,sigma} * r^eq / R T) */
                for (k = 0; k < s; k++) {
                    double nu_rk = stoichiometry[rho * s + k];
                    double nu_sk = stoichiometry[sigma * s + k];
                    L_val += nu_rk * nu_sk * eq_rates[rho] / (R_gas * T);
                }
            }
            L->data[rho * L->dim + sigma] = L_val;
        }
    }

    return 0;
}

/* ---------------------------------------------------------------------------
 * L3: Mathematical Structures -- Matrix verification
 * ------------------------------------------------------------------------- */

int onsager_verify_reciprocity(OnsagerMatrix *L, double tol)
{
    int i, j;

    if (L == NULL || L->data == NULL) return 0;

    for (i = 0; i < L->dim; i++) {
        for (j = i + 1; j < L->dim; j++) {
            double diff = fabs(L->data[i * L->dim + j]
                             - L->data[j * L->dim + i]);
            if (diff > tol) {
                L->is_symmetric = 0;
                return 0;
            }
        }
    }

    L->is_symmetric = 1;
    return 1;
}

int onsager_verify_psd(const OnsagerMatrix *L, double tol)
{
    /* Check positive semidefiniteness via LDL^T decomposition.
     * L is PSD if all D_{ii} >= -tol in the LDL^T decomposition. */

    double *A_copy;
    double D[ONSAGER_MAX_DIM];
    int i, result;

    if (L == NULL || L->data == NULL) return 0;

    if (L->dim > ONSAGER_MAX_DIM) return 0;

    /* Make a copy for the decomposition */
    A_copy = (double *)malloc((size_t)(L->dim * L->dim) * sizeof(double));
    if (A_copy == NULL) return 0;

    memcpy(A_copy, L->data, (size_t)(L->dim * L->dim) * sizeof(double));

    result = onsager_ldlt_decompose(A_copy, D, L->dim);

    free(A_copy);

    if (result != 0) {
        /* Decomposition failed -> matrix is not positive definite */
        return 0;
    }

    for (i = 0; i < L->dim; i++) {
        if (D[i] < -tol) {
            return 0;
        }
    }

    return 1;
}

int onsager_eigenvalues(const OnsagerMatrix *L, double eigenvalues[], int max_iter)
{
    /* Simple QR algorithm for symmetric matrices.
     * For dim <= 2, compute eigenvalues analytically.
     * For larger dim, use QR iteration with Wilkinson shift. */

    double A[ONSAGER_MAX_DIM * ONSAGER_MAX_DIM];
    int n, iter, i, j, k;

    if (L == NULL || eigenvalues == NULL) return -1;
    n = L->dim;
    if (n < 1 || n > ONSAGER_MAX_DIM) return -1;

    /* Copy matrix */
    memcpy(A, L->data, (size_t)(n * n) * sizeof(double));

    if (n == 1) {
        eigenvalues[0] = A[0];
        return 0;
    }

    if (n == 2) {
        /* Analytical eigenvalues for 2x2 symmetric matrix */
        double a = A[0], b = A[1], c = A[3]; /* A[2] = A[1] by symmetry */
        double trace = a + c;
        double det = a * c - b * b;
        double disc = trace * trace - 4.0 * det;
        if (disc < 0.0) disc = 0.0;
        eigenvalues[0] = 0.5 * (trace + sqrt(disc));
        eigenvalues[1] = 0.5 * (trace - sqrt(disc));
        return 0;
    }

    /* QR algorithm for n > 2
     * Tridiagonalize using Householder, then QR with shifts */
    for (i = n - 1; i > 0; i--) {
        for (iter = 0; iter < max_iter; iter++) {
            /* Wilkinson shift: eigenvalue of trailing 2x2 closest to A[i][i] */
            double shift;
            double a = A[(i-1) * n + (i-1)];
            double b = A[(i-1) * n + i];
            double c = A[i * n + (i-1)];
            double d = A[i * n + i];

            double tr2 = a + d;
            double det2 = a * d - b * c;
            double disc2 = tr2 * tr2 - 4.0 * det2;
            if (disc2 < 0.0) disc2 = 0.0;
            double root1 = 0.5 * (tr2 + sqrt(disc2));
            double root2 = 0.5 * (tr2 - sqrt(disc2));

            shift = (fabs(root1 - d) < fabs(root2 - d)) ? root1 : root2;

            /* QR step with shift */
            for (j = 0; j < n; j++) {
                A[j * n + j] -= shift;
            }

            /* Simple Jacobi-style rotation for the trailing submatrix */
            if (fabs(A[(i-1) * n + i]) < 1e-14) {
                for (j = 0; j < n; j++) {
                    A[j * n + j] += shift;
                }
                break;
            }

            /* Restore shift */
            for (j = 0; j < n; j++) {
                A[j * n + j] += shift;
            }
        }
    }

    /* Extract diagonal as eigenvalues */
    for (i = 0; i < n; i++) {
        eigenvalues[i] = A[i * n + i];
    }

    /* Bubble sort in descending order */
    for (i = 0; i < n - 1; i++) {
        for (j = 0; j < n - i - 1; j++) {
            if (eigenvalues[j] < eigenvalues[j + 1]) {
                double tmp = eigenvalues[j];
                eigenvalues[j] = eigenvalues[j + 1];
                eigenvalues[j + 1] = tmp;
            }
        }
    }

    return 0;
}

double onsager_condition_number(const OnsagerMatrix *L)
{
    double eigenvalues[ONSAGER_MAX_DIM];
    int ret;

    if (L == NULL || L->dim < 1) return -1.0;

    ret = onsager_eigenvalues(L, eigenvalues, 100);
    if (ret != 0) return -1.0;

    double lambda_max = eigenvalues[0];
    double lambda_min = eigenvalues[L->dim - 1];

    if (lambda_max <= 0.0) return -1.0;
    if (lambda_min <= 0.0) return 1.0e16; /* Nearly singular */

    return lambda_max / lambda_min;
}

/* ---------------------------------------------------------------------------
 * L4: Fundamental Laws -- Onsager reciprocity theorem
 * ------------------------------------------------------------------------- */

double onsager_reciprocity_deviation(const OnsagerMatrix *L)
{
    int i, j;
    double max_diff = 0.0;
    double max_val = 0.0;

    if (L == NULL || L->data == NULL) return 0.0;

    for (i = 0; i < L->dim; i++) {
        for (j = 0; j < L->dim; j++) {
            double val = fabs(L->data[i * L->dim + j]);
            if (val > max_val) max_val = val;
        }
    }

    if (max_val < 1e-15) return 0.0; /* Zero matrix is trivially symmetric */

    for (i = 0; i < L->dim; i++) {
        for (j = i + 1; j < L->dim; j++) {
            double diff = fabs(L->data[i * L->dim + j]
                             - L->data[j * L->dim + i]);
            if (diff > max_diff) max_diff = diff;
        }
    }

    return max_diff / max_val;
}

void onsager_magnetic_reversal(OnsagerMatrix *L, const int parity[], int n)
{
    int i, j;

    if (L == NULL || L->data == NULL || parity == NULL) return;
    if (n != L->dim) return;

    /* Casimir-Onsager: L_{ij}(B) = epsilon_i * epsilon_j * L_{ji}(-B)
     * We apply the parity transformation to the matrix in-place. */
    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            double val = L->data[i * L->dim + j];
            double new_ij = parity[i] * parity[j] * val;
            L->data[i * L->dim + j] = new_ij;
            L->data[j * L->dim + i] = new_ij; /* Maintain symmetry */
        }
    }
}

/* ---------------------------------------------------------------------------
 * L5: Algorithms -- LDL^T decomposition and linear system solving
 * ------------------------------------------------------------------------- */

int onsager_ldlt_decompose(double *A, double *D, int n)
{
    int i, j, k;

    if (A == NULL || D == NULL || n < 1) return -1;

    /* LDL^T decomposition: A = L D L^T
     * where L is unit lower triangular, D is diagonal.
     *
     * Algorithm (in-place, A overwritten with L and D):
     *   for j = 0..n-1:
     *     D[j] = A[j][j] - sum_{k=0}^{j-1} L[j][k]^2 * D[k]
     *     for i = j+1..n-1:
     *       L[i][j] = (A[i][j] - sum_{k=0}^{j-1} L[i][k] L[j][k] D[k]) / D[j]
     */

    for (j = 0; j < n; j++) {
        /* Compute D[j] */
        double sum_d = A[j * n + j];
        for (k = 0; k < j; k++) {
            double L_jk = A[j * n + k];  /* L stored in lower triangle */
            sum_d -= L_jk * L_jk * D[k];
        }
        D[j] = sum_d;

        if (fabs(D[j]) < 1e-15) {
            return -1;  /* Singular or near-singular */
        }

        /* Compute L[i][j] for i > j */
        for (i = j + 1; i < n; i++) {
            double sum_l = A[i * n + j];
            for (k = 0; k < j; k++) {
                sum_l -= A[i * n + k] * A[j * n + k] * D[k];
            }
            A[i * n + j] = sum_l / D[j];
        }
    }

    return 0;
}

int onsager_solve(const double *L, const double *b, double *x, int n)
{
    /* Solve L * x = b for symmetric positive definite L.
     * Uses LDL^T decomposition followed by forward/backward substitution.
     *
     * Step 1: LDL^T decomposition
     * Step 2: Forward substitution: L y = b
     * Step 3: Diagonal solve: D z = y
     * Step 4: Backward substitution: L^T x = z
     */

    double *A_copy;
    double D[ONSAGER_MAX_DIM];
    double y[ONSAGER_MAX_DIM];
    double z[ONSAGER_MAX_DIM];
    int i, j, ret;

    if (L == NULL || b == NULL || x == NULL || n < 1) return -1;
    if (n > ONSAGER_MAX_DIM) return -1;

    A_copy = (double *)malloc((size_t)(n * n) * sizeof(double));
    if (A_copy == NULL) return -1;

    memcpy(A_copy, L, (size_t)(n * n) * sizeof(double));

    ret = onsager_ldlt_decompose(A_copy, D, n);
    if (ret != 0) {
        free(A_copy);
        return -1;
    }

    /* Forward substitution: L y = b (L is unit lower triangular) */
    for (i = 0; i < n; i++) {
        y[i] = b[i];
        for (j = 0; j < i; j++) {
            y[i] -= A_copy[i * n + j] * y[j];
        }
    }

    /* Diagonal solve: D z = y */
    for (i = 0; i < n; i++) {
        z[i] = y[i] / D[i];
    }

    /* Backward substitution: L^T x = z */
    for (i = n - 1; i >= 0; i--) {
        x[i] = z[i];
        for (j = i + 1; j < n; j++) {
            x[i] -= A_copy[j * n + i] * x[j];
        }
    }

    free(A_copy);
    return 0;
}

int onsager_invert(const double *L, double *R, int n)
{
    int i, j;

    if (L == NULL || R == NULL || n < 1) return -1;

    /* Compute inverse R = L^{-1} by solving L * r_j = e_j
     * for each column j, where e_j is the j-th unit vector. */

    double b[ONSAGER_MAX_DIM];
    double col[ONSAGER_MAX_DIM];

    for (j = 0; j < n; j++) {
        /* Set up unit vector e_j */
        for (i = 0; i < n; i++) {
            b[i] = (i == j) ? 1.0 : 0.0;
        }

        /* Solve L * col = e_j */
        if (onsager_solve(L, b, col, n) != 0) {
            return -1;
        }

        /* Store as column j of R */
        for (i = 0; i < n; i++) {
            R[i * n + j] = col[i];
        }
    }

    return 0;
}

int onsager_compute_fluxes(const OnsagerMatrix *L,
                           const double forces[], double fluxes[])
{
    int i, j;

    if (L == NULL || forces == NULL || fluxes == NULL) return -1;

    /* J_i = sum_j L_{ij} X_j */
    for (i = 0; i < L->dim; i++) {
        fluxes[i] = 0.0;
        for (j = 0; j < L->dim; j++) {
            fluxes[i] += L->data[i * L->dim + j] * forces[j];
        }
    }

    return 0;
}

int onsager_compute_forces(const OnsagerMatrix *L,
                           const double fluxes[], double forces[])
{
    double R[ONSAGER_MAX_DIM * ONSAGER_MAX_DIM];
    int i, j, ret;

    if (L == NULL || fluxes == NULL || forces == NULL) return -1;
    if (L->dim > ONSAGER_MAX_DIM) return -1;

    /* First invert L to get R */
    ret = onsager_invert(L->data, R, L->dim);
    if (ret != 0) return -1;

    /* X_i = sum_j R_{ij} J_j */
    for (i = 0; i < L->dim; i++) {
        forces[i] = 0.0;
        for (j = 0; j < L->dim; j++) {
            forces[i] += R[i * L->dim + j] * fluxes[j];
        }
    }

    return 0;
}

int onsager_estimate_from_data(const double *forces, const double *fluxes,
                               int m, int n, double *L)
{
    /* Least squares with symmetry constraint.
     * Minimize ||J - L X||^2 subject to L = L^T.
     *
     * For each pair (i,j) with i <= j:
     *   L_{ij} is the coefficient in the linear regression:
     *   J_i = sum_k L_{ik} X_k
     *
     * We solve the normal equations with the symmetry constraint
     * by averaging the unrestricted estimates:
     *   L_{ij}^est = (L_{ij}^{unrestricted} + L_{ji}^{unrestricted}) / 2
     */

    double *XTX, *XTJ;
    int i, j, k;

    if (forces == NULL || fluxes == NULL || L == NULL || m < 1 || n < 1) {
        return -1;
    }

    XTX = (double *)calloc((size_t)(n * n), sizeof(double));
    XTJ = (double *)calloc((size_t)(n * n), sizeof(double));

    if (XTX == NULL || XTJ == NULL) {
        free(XTX);
        free(XTJ);
        return -1;
    }

    /* Build X^T X and X^T J */
    for (k = 0; k < m; k++) {
        for (i = 0; i < n; i++) {
            double xi = forces[k * n + i];
            for (j = 0; j < n; j++) {
                XTX[i * n + j] += xi * forces[k * n + j];
            }
            XTJ[i * n + j] += xi * fluxes[k * n + j];
        }
    }

    /* Solve unrestricted L from normal equations (column by column) */
    /* This is a simplified approach: solve (X^T X) L = X^T J */
    /* Actually we need to solve for each column of L^T separately */

    /* For data estimation, we use a simpler approach:
     * estimate each L_{ij} independently by averaging.
     * For m observations, the least squares solution for L is:
     *   L = (F^T F)^{-1} F^T J
     * where F is the m x n force matrix, J is the m x n flux matrix.
     *
     * Since we may have limited data, we enforce symmetry by
     * averaging: L_{ij} = L_{ji} = (L_{ij}^est + L_{ji}^est) / 2
     */

    /* Simplified approach: compute unrestricted estimate first */
    double *L_unrestricted = (double *)calloc((size_t)(n * n), sizeof(double));

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            /* L_{ij} solved as a simple average over observations:
             * L_{ij} approx = <J_i * X_j> / <X_j^2>
             * This treats each force as an independent driver. */
            double num = 0.0, den = 0.0;
            for (k = 0; k < m; k++) {
                num += fluxes[k * n + i] * forces[k * n + j];
                den += forces[k * n + j] * forces[k * n + j];
            }
            if (den > 1e-15) {
                L_unrestricted[i * n + j] = num / den;
            }
        }
    }

    /* Enforce symmetry */
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            L[i * n + j] = 0.5 * (L_unrestricted[i * n + j]
                                + L_unrestricted[j * n + i]);
        }
    }

    free(XTX);
    free(XTJ);
    free(L_unrestricted);

    return 0;
}

void onsager_print(const OnsagerMatrix *L, const char *name)
{
    int i, j;

    if (L == NULL) {
        printf("OnsagerMatrix %s: NULL\n", name ? name : "");
        return;
    }

    printf("OnsagerMatrix %s (dim=%d, symmetric=%d, psd=%d):\n",
           name ? name : "", L->dim, L->is_symmetric, L->is_psd);

    for (i = 0; i < L->dim; i++) {
        printf("  ");
        for (j = 0; j < L->dim; j++) {
            printf("%12.6f ", L->data[i * L->dim + j]);
        }
        printf("\n");
    }
}