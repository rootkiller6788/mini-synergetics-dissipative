/* ==============================================================
 * haken_mode_decomp.c — Mode Decomposition and Spectral Analysis
 *
 * Implements the complete spectral decomposition pipeline essential
 * for applying Haken's slaving principle: matrix balancing,
 * Hessenberg reduction, QR decomposition, eigenvalue computation,
 * eigenvector refinement, spectral gap analysis, and participation
 * factors.
 *
 * References:
 *   Golub & Van Loan (2013) Matrix Computations, 4th ed. §7.4-7.6
 *   Wilkinson (1965) The Algebraic Eigenvalue Problem
 *   Osborne (1960) J. ACM 7, 338-345 (balancing)
 *   Trefethen & Embree (2005) Spectra and Pseudospectra
 * ============================================================== */

#include "haken_mode_decomp.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* ==============================================================
 * L1: Matrix balancing (Osborne 1960)
 *
 * Balances a matrix by diagonal similarity transformation to reduce
 * the norm and improve eigenvalue accuracy. Finds D such that
 * D^{-1} A D has approximately equal row and column norms.
 * ============================================================== */

void haken_balance_matrix(const double* A, int n, double* balanced, double* D) {
    if (!A || n <= 0 || !balanced) return;

    /* Initialize D = I, copy A to balanced */
    for (int i = 0; i < n; i++) { D[i] = 1.0; }
    memcpy(balanced, A, (size_t)(n * n) * sizeof(double));

    const double radix = 2.0;
    const double sq_radix = radix * radix;
    int converged = 0;

    while (!converged) {
        converged = 1;
        for (int i = 0; i < n; i++) {
            double row_norm = 0.0, col_norm = 0.0;
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    row_norm += fabs(balanced[i * n + j]);
                    col_norm += fabs(balanced[j * n + i]);
                }
            }
            if (row_norm == 0.0 || col_norm == 0.0) continue;

            double f = 1.0;
            double ratio = col_norm / row_norm;
            while (ratio >= sq_radix) { f *= radix; ratio /= sq_radix; }
            while (ratio < 1.0 / sq_radix) { f /= radix; ratio *= sq_radix; }

            if (fabs(f - 1.0) > 1e-10) {
                converged = 0;
                D[i] *= f;
                for (int j = 0; j < n; j++) balanced[i * n + j] *= (1.0 / f);
                for (int j = 0; j < n; j++) balanced[j * n + i] *= f;
            }
        }
    }
}

/* ==============================================================
 * L2: Hessenberg reduction via Householder reflections
 * ============================================================== */

void haken_hessenberg(const double* A, int n, double* H, double* Q) {
    if (!A || n <= 0 || !H || !Q) return;

    /* Copy A to H */
    memcpy(H, A, (size_t)(n * n) * sizeof(double));

    /* Initialize Q = I */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) Q[i * n + j] = 0.0;
        Q[i * n + i] = 1.0;
    }

    double* v = (double*)calloc((size_t)n, sizeof(double));

    for (int k = 0; k < n - 2; k++) {
        /* Compute norm of column k below subdiagonal (from k+1 to n-1) */
        double sigma = 0.0;
        for (int i = k + 1; i < n; i++) sigma += H[i * n + k] * H[i * n + k];

        if (sigma < 1e-30) continue;

        double xk = H[(k + 1) * n + k];
        double mu = sqrt(xk * xk + sigma);
        double alpha = (xk > 0) ? -mu : mu;

        /* Build Householder vector: v[k+1] = xk - alpha, v[i] = H[i,k] for i > k+1 */
        v[k + 1] = xk - alpha;
        for (int i = k + 2; i < n; i++) v[i] = H[i * n + k];

        /* β = 2 / (v^T v) */
        double vtv = v[k + 1] * v[k + 1];
        for (int i = k + 2; i < n; i++) vtv += v[i] * v[i];
        double beta = 2.0 / vtv;

        /* Apply from left: H = (I - β·v·v^T)·H */
        for (int j = k; j < n; j++) {
            double s = 0.0;
            for (int i = k + 1; i < n; i++) s += v[i] * H[i * n + j];
            s *= beta;
            for (int i = k + 1; i < n; i++) H[i * n + j] -= s * v[i];
        }

        /* Apply from right: H = H·(I - β·v·v^T) */
        for (int i = 0; i < n; i++) {
            double s = 0.0;
            for (int j = k + 1; j < n; j++) s += H[i * n + j] * v[j];
            s *= beta;
            for (int j = k + 1; j < n; j++) H[i * n + j] -= s * v[j];
        }

        /* Accumulate Q = Q·(I - β·v·v^T) */
        for (int i = 0; i < n; i++) {
            double s = 0.0;
            for (int j = k + 1; j < n; j++) s += Q[i * n + j] * v[j];
            s *= beta;
            for (int j = k + 1; j < n; j++) Q[i * n + j] -= s * v[j];
        }

        /* Store alpha in H[k+1,k] */
        H[(k + 1) * n + k] = alpha;
    }
    free(v);
}

/* ==============================================================
 * L3: QR decomposition via Householder reflections
 * ============================================================== */

void haken_qr_decompose(const double* A, int rows, int cols,
                         double* Q, double* R) {
    if (!A || rows <= 0 || cols <= 0 || !Q || !R) return;

    /* Initialize Q = I (rows × rows) */
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < rows; j++) Q[i * rows + j] = 0.0;
        Q[i * rows + i] = 1.0;
    }
    /* Copy A to R (rows × cols) */
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            R[i * cols + j] = A[i * cols + j];

    double* v = (double*)calloc((size_t)rows, sizeof(double));
    double* w = (double*)calloc((size_t)cols, sizeof(double));
    int min_mn = (rows < cols) ? rows : cols;

    for (int k = 0; k < min_mn; k++) {
        /* Norm of the k-th column below and including diagonal */
        double sigma = 0.0;
        for (int i = k; i < rows; i++) sigma += R[i * cols + k] * R[i * cols + k];
        if (sigma < 1e-30) continue;

        double xk = R[k * cols + k];
        double mu = sqrt(sigma);
        double alpha = (xk > 0) ? -mu : mu;

        /* Householder vector v: normalized so v[k] = xk - alpha, v[i] = x[i] for i > k */
        v[k] = xk - alpha;
        for (int i = k + 1; i < rows; i++) v[i] = R[i * cols + k];

        /* β = 2 / (v^T v) */
        double vtv = v[k] * v[k];
        for (int i = k + 1; i < rows; i++) vtv += v[i] * v[i];
        double beta = 2.0 / vtv;

        /* Apply H = I - β·v·v^T to R from the left */
        for (int j = k; j < cols; j++) {
            double s = 0.0;
            for (int i = k; i < rows; i++) s += v[i] * R[i * cols + j];
            s *= beta;
            for (int i = k; i < rows; i++) R[i * cols + j] -= s * v[i];
        }

        /* Apply H^T = H to Q from the right: Q = Q·H */
        for (int i = 0; i < rows; i++) {
            double s = 0.0;
            for (int j = k; j < rows; j++) s += Q[i * rows + j] * v[j];
            s *= beta;
            for (int j = k; j < rows; j++) Q[i * rows + j] -= s * v[j];
        }

        /* Store the reflection vector in R below diagonal */
        R[k * cols + k] = alpha;
    }
    free(v); free(w);
}

/* ==============================================================
 * Matrix multiplication
 * ============================================================== */

void haken_matrix_multiply(const double* A, const double* B,
                            int n, double* C) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int k = 0; k < n; k++) sum += A[i * n + k] * B[k * n + j];
            C[i * n + j] = sum;
        }
    }
}

void haken_matrix_transpose(const double* A, int n, double* At) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            At[j * n + i] = A[i * n + j];
}

/* ==============================================================
 * L4: Eigenvalue computation with QR algorithm
 * ============================================================== */

int haken_qr_eigenvalues(double* A, int n,
                          double* eigenvalues_re,
                          double* eigenvalues_im,
                          const Haken_QRConfig* config) {
    if (!A || n <= 0 || !eigenvalues_re || !eigenvalues_im) return -1;

    int max_iter = config ? config->max_iterations : 100;
    double tol   = config ? config->tolerance : 1e-12;

    /* Balance */
    double* work = (double*)malloc((size_t)(n * n) * sizeof(double));
    double* D    = (double*)malloc((size_t)n * sizeof(double));
    haken_balance_matrix(A, n, work, D);

    /* Hessenberg */
    double* Q = (double*)calloc((size_t)(n * n), sizeof(double));
    double* H = (double*)malloc((size_t)(n * n) * sizeof(double));
    haken_hessenberg(work, n, H, Q);

    /* QR iterations */
    int m = n;
    for (int iter = 0; iter < max_iter * n; iter++) {
        /* Find the largest m such that subdiagonal is negligible */
        while (m > 1) {
            if (fabs(H[(m - 1) * n + (m - 2)]) < tol) m--;
            else if (m > 2 && fabs(H[(m - 2) * n + (m - 3)]) < tol) m--;
            else break;
        }
        if (m <= 1) break; /* Fully converged */

        /* Compute shift from trailing 2×2 */
        int p = m - 1;
        double a11 = H[(p - 1) * n + (p - 1)];
        double a12 = H[(p - 1) * n + p];
        double a21 = H[p * n + (p - 1)];
        double a22 = H[p * n + p];

        /* Choose shift closer to a22 */
        double tr = a11 + a22;
        double det = a11 * a22 - a12 * a21;
        double disc = tr * tr - 4.0 * det;
        double shift;
        if (disc >= 0) {
            double sd = sqrt(disc);
            double r1 = (tr + sd) / 2.0, r2 = (tr - sd) / 2.0;
            shift = (fabs(r1 - a22) < fabs(r2 - a22)) ? r1 : r2;
        } else {
            shift = tr / 2.0;
        }

        /* Apply single shift QR step to H(0:m-1, 0:m-1) */
        for (int j = 0; j < m; j++) H[j * n + j] -= shift;
        for (int k = 0; k < m - 1; k++) {
            double x = H[k * n + k], y = H[(k + 1) * n + k];
            double r = sqrt(x * x + y * y);
            if (r < 1e-30) continue;
            double c = x / r, s = y / r;
            for (int j = k; j < m; j++) {
                double t1 =  c * H[k * n + j] + s * H[(k + 1) * n + j];
                double t2 = -s * H[k * n + j] + c * H[(k + 1) * n + j];
                H[k * n + j] = t1; H[(k + 1) * n + j] = t2;
            }
            for (int i = 0; i <= k + 1 && i < m; i++) {
                double t1 =  c * H[i * n + k] + s * H[i * n + (k + 1)];
                double t2 = -s * H[i * n + k] + c * H[i * n + (k + 1)];
                H[i * n + k] = t1; H[i * n + (k + 1)] = t2;
            }
        }
        for (int j = 0; j < m; j++) H[j * n + j] += shift;
    }

    /* Extract eigenvalues */
    int idx = 0;
    int i = 0;
    while (i < n) {
        if (i == n - 1 || fabs(H[(i + 1) * n + i]) < tol) {
            eigenvalues_re[idx] = H[i * n + i];
            eigenvalues_im[idx] = 0.0;
            idx++; i++;
        } else {
            double a = H[i * n + i], b = H[i * n + (i + 1)];
            double c = H[(i + 1) * n + i], d = H[(i + 1) * n + (i + 1)];
            double trace = a + d, determinant = a * d - b * c;
            double discriminant = trace * trace - 4.0 * determinant;
            if (discriminant >= 0) {
                double sd = sqrt(discriminant);
                eigenvalues_re[idx] = (trace + sd) / 2.0; eigenvalues_im[idx] = 0.0; idx++;
                eigenvalues_re[idx] = (trace - sd) / 2.0; eigenvalues_im[idx] = 0.0; idx++;
            } else {
                eigenvalues_re[idx] = trace / 2.0;
                eigenvalues_im[idx] = sqrt(-discriminant) / 2.0; idx++;
                eigenvalues_re[idx] = trace / 2.0;
                eigenvalues_im[idx] = -sqrt(-discriminant) / 2.0; idx++;
            }
            i += 2;
        }
    }

    free(work); free(D); free(Q); free(H);
    return 0;
}

/* ==============================================================
 * Inverse iteration for eigenvector refinement
 * (Wilkinson 1965, §9)
 * ============================================================== */

void haken_inverse_iteration(const double* A, int n,
                              double lambda_re, double lambda_im,
                              double* eigenvector_re,
                              double* eigenvector_im_,
                              int max_iter, double tol) {
    (void)lambda_im; (void)eigenvector_im_;
    if (!A || n <= 0 || !eigenvector_re) return;

    /* Build shifted matrix B = A - λI */
    double* B = (double*)malloc((size_t)(n * n) * sizeof(double));
    memcpy(B, A, (size_t)(n * n) * sizeof(double));
    for (int i = 0; i < n; i++) B[i * n + i] -= lambda_re;

    /* LU decomposition of B */
    for (int k = 0; k < n; k++) {
        if (fabs(B[k * n + k]) < 1e-14) B[k * n + k] = 1e-14;
        for (int i = k + 1; i < n; i++) {
            B[i * n + k] /= B[k * n + k];
            for (int j = k + 1; j < n; j++)
                B[i * n + j] -= B[i * n + k] * B[k * n + j];
        }
    }

    /* Normalize initial guess */
    double norm = 0.0;
    for (int i = 0; i < n; i++) norm += eigenvector_re[i] * eigenvector_re[i];
    norm = sqrt(norm);
    if (norm < tol) {
        for (int i = 0; i < n; i++) eigenvector_re[i] = 1.0;
        norm = sqrt((double)n);
    }
    for (int i = 0; i < n; i++) eigenvector_re[i] /= norm;

    for (int iter = 0; iter < max_iter; iter++) {
        /* Solve B y = x (forward substitution) */
        double* y = (double*)malloc((size_t)n * sizeof(double));
        memcpy(y, eigenvector_re, (size_t)n * sizeof(double));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) y[i] -= B[i * n + j] * y[j];
        }
        /* Back substitution */
        for (int i = n - 1; i >= 0; i--) {
            for (int j = i + 1; j < n; j++) y[i] -= B[i * n + j] * y[j];
            y[i] /= B[i * n + i];
        }
        /* Normalize */
        norm = 0.0;
        for (int i = 0; i < n; i++) norm += y[i] * y[i];
        norm = sqrt(norm);
        if (norm < 1e-30) { free(y); break; }
        double diff = 0.0;
        for (int i = 0; i < n; i++) {
            double old = eigenvector_re[i];
            eigenvector_re[i] = y[i] / norm;
            diff += (eigenvector_re[i] - old) * (eigenvector_re[i] - old);
        }
        free(y);
        if (sqrt(diff) < tol) break;
    }
    free(B);
}

void haken_left_eigenvectors(const double* A, int n,
                              const Haken_Eigenpair* right_pairs,
                              Haken_Eigenpair* left_pairs) {
    if (!A || n <= 0 || !right_pairs || !left_pairs) return;

    /* Left eigenvectors of A = right eigenvectors of A^T */
    double* At = (double*)malloc((size_t)(n * n) * sizeof(double));
    haken_matrix_transpose(A, n, At);

    for (int k = 0; k < n; k++) {
        left_pairs[k].lambda_re = right_pairs[k].lambda_re;
        left_pairs[k].lambda_im = right_pairs[k].lambda_im;
        /* Compute via inverse iteration on A^T */
        for (int i = 0; i < n; i++) left_pairs[k].vector_re[i] = 1.0;
        haken_inverse_iteration(At, n,
            right_pairs[k].lambda_re, right_pairs[k].lambda_im,
            left_pairs[k].vector_re, left_pairs[k].vector_im, 20, 1e-10);
    }
    free(At);
}

/* ==============================================================
 * L5: Spectral gap analysis
 * ============================================================== */

double haken_spectral_gap_decomp(const Haken_ModeDecomp* decomp,
                                   double eps_critical) {
    if (!decomp || decomp->n_modes <= 0) return 0.0;

    double max_crit = -INFINITY, max_stable = -INFINITY;
    int has_crit = 0, has_stable = 0;

    for (int k = 0; k < decomp->n_modes; k++) {
        double re = decomp->modes[k].lambda_re;
        if (fabs(re) <= eps_critical) {
            if (re > max_crit) max_crit = re;
            has_crit = 1;
        } else if (re < 0) {
            if (re > max_stable) max_stable = re;
            has_stable = 1;
        }
    }
    if (!has_crit || !has_stable) return 0.0;
    return max_crit - max_stable;
}

bool haken_timescale_separation(const Haken_ModeDecomp* decomp,
                                 double min_gap_ratio) {
    if (!decomp) return false;
    double gap = decomp->spectral_gap;
    if (gap <= 0) return false;

    /* Compute timescale ratio: τ_slow / τ_fast */
    double max_stable_re = -INFINITY;
    for (int k = 0; k < decomp->n_modes; k++) {
        if (decomp->modes[k].type == HAKEN_MODE_STABLE) {
            if (decomp->modes[k].lambda_re > max_stable_re)
                max_stable_re = decomp->modes[k].lambda_re;
        }
    }
    if (max_stable_re >= 0) return false;
    return (gap / fabs(max_stable_re)) > min_gap_ratio;
}

void haken_participation_factors(const Haken_Eigenpair* right,
                                  const Haken_Eigenpair* left,
                                  int n, double* p_matrix) {
    if (!right || !left || n <= 0 || !p_matrix) return;

    for (int k = 0; k < n; k++) {
        double sum = 0.0;
        for (int i = 0; i < n; i++) {
            double pk = right[k].vector_re[i] * left[k].vector_re[i];
            p_matrix[i * n + k] = pk;
            sum += pk;
        }
        if (sum > 1e-15) {
            for (int i = 0; i < n; i++) p_matrix[i * n + k] /= sum;
        }
    }
}

void haken_sort_modes_by_timescale(Haken_Mode* modes, int n) {
    if (!modes || n <= 1) return;
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            if (modes[j].relaxation_time > modes[i].relaxation_time) {
                Haken_Mode tmp = modes[i];
                modes[i] = modes[j];
                modes[j] = tmp;
            }
        }
    }
}

int haken_find_most_critical_mode(const Haken_Mode* modes, int n) {
    if (!modes || n <= 0) return -1;
    int best = 0;
    double min_abs = fabs(modes[0].lambda_re);
    for (int k = 1; k < n; k++) {
        double abs_re = fabs(modes[k].lambda_re);
        if (abs_re < min_abs) { min_abs = abs_re; best = k; }
    }
    return best;
}

double haken_eigenvector_condition(const double* eigenvectors, int n) {
    if (!eigenvectors || n <= 0) return INFINITY;
    /* Compute Frobenius norm of V and V^{-1} */
    double norm_V = 0.0;
    for (int i = 0; i < n * n; i++) norm_V += eigenvectors[i] * eigenvectors[i];
    norm_V = sqrt(norm_V);
    /* Estimate norm of inverse via condition of V^T V */
    double norm_inv = 0.0;
    double* VtV = (double*)calloc((size_t)(n * n), sizeof(double));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++)
                VtV[i * n + j] += eigenvectors[k * n + i] * eigenvectors[k * n + j];
    for (int i = 0; i < n; i++) {
        double diag = VtV[i * n + i];
        if (diag > 1e-15) norm_inv += 1.0 / diag;
    }
    norm_inv = sqrt(norm_inv);
    free(VtV);
    return (norm_V > 0 && norm_inv < INFINITY) ? norm_V * norm_inv : INFINITY;
}

double haken_resolvent_norm(const double* A, int n,
                             double lambda_re, double lambda_im_) {
    (void)lambda_im_;
    if (!A || n <= 0) return INFINITY;
    /* Compute ||(A - λI)^{-1}||_F approximately */
    double* B = (double*)malloc((size_t)(n * n) * sizeof(double));
    memcpy(B, A, (size_t)(n * n) * sizeof(double));
    for (int i = 0; i < n; i++) B[i * n + i] -= lambda_re;

    /* Estimate via QR decomposition and back-substitution */
    double norm = 0.0;
    double* y = (double*)calloc((size_t)n, sizeof(double));
    for (int trial = 0; trial < 10; trial++) {
        for (int j = 0; j < n; j++) y[j] = (double)(rand()) / RAND_MAX;
        /* Solve B x = y by Gaussian elimination with partial pivoting */
        double* Bcopy = (double*)malloc((size_t)(n * n) * sizeof(double));
        memcpy(Bcopy, B, (size_t)(n * n) * sizeof(double));
        double* x = (double*)malloc((size_t)n * sizeof(double));
        memcpy(x, y, (size_t)n * sizeof(double));
        for (int k = 0; k < n; k++) {
            int max_row = k;
            double max_val = fabs(Bcopy[k * n + k]);
            for (int r = k + 1; r < n; r++) {
                if (fabs(Bcopy[r * n + k]) > max_val) {
                    max_val = fabs(Bcopy[r * n + k]); max_row = r;
                }
            }
            if (max_val < 1e-15) { free(x); free(Bcopy); continue; }
            if (max_row != k) {
                for (int j = 0; j < n; j++) {
                    double t = Bcopy[k * n + j]; Bcopy[k * n + j] = Bcopy[max_row * n + j]; Bcopy[max_row * n + j] = t; }
                double t = x[k]; x[k] = x[max_row]; x[max_row] = t;
            }
            for (int r = k + 1; r < n; r++) {
                double factor = Bcopy[r * n + k] / Bcopy[k * n + k];
                for (int j = k; j < n; j++) Bcopy[r * n + j] -= factor * Bcopy[k * n + j];
                x[r] -= factor * x[k];
            }
        }
        for (int k = n - 1; k >= 0; k--) {
            for (int j = k + 1; j < n; j++) x[k] -= Bcopy[k * n + j] * x[j];
            if (fabs(Bcopy[k * n + k]) > 1e-15) x[k] /= Bcopy[k * n + k];
        }
        double xn = 0.0;
        for (int j = 0; j < n; j++) xn += x[j] * x[j];
        if (sqrt(xn) > norm) norm = sqrt(xn);
        free(x); free(Bcopy);
    }
    free(B); free(y);
    return norm;
}

int haken_pseudospectrum(const double* A, int n, double epsilon,
                          int grid_size, double range,
                          double* pseudo_grid) {
    if (!A || n <= 0 || !pseudo_grid || grid_size <= 0) return -1;
    double step = 2.0 * range / (grid_size + 1);
    int count = 0;
    for (int ix = 0; ix < grid_size; ix++) {
        double rx = -range + step * (ix + 1);
        for (int iy = 0; iy < grid_size; iy++) {
            double ry = -range + step * (iy + 1);
            double rnorm = haken_resolvent_norm(A, n, rx, ry);
            if (rnorm > 1.0 / epsilon) {
                if (count < grid_size * grid_size * 3) {
                    pseudo_grid[count * 3 + 0] = rx;
                    pseudo_grid[count * 3 + 1] = ry;
                    pseudo_grid[count * 3 + 2] = rnorm;
                    count++;
                }
            }
        }
    }
    return count;
}

int haken_compute_eigenpairs(const double* A, int n,
                              Haken_Eigenpair* pairs,
                              const Haken_QRConfig* config) {
    if (!A || n <= 0 || !pairs) return -1;
    double* work = (double*)malloc((size_t)(n * n) * sizeof(double));
    memcpy(work, A, (size_t)(n * n) * sizeof(double));
    double* re = (double*)malloc((size_t)n * sizeof(double));
    double* im = (double*)malloc((size_t)n * sizeof(double));
    int ret = haken_qr_eigenvalues(work, n, re, im, config);
    if (ret != 0) { free(work); free(re); free(im); return ret; }
    for (int k = 0; k < n; k++) {
        pairs[k].lambda_re = re[k];
        pairs[k].lambda_im = im[k];
        if (!pairs[k].vector_re) pairs[k].vector_re = (double*)calloc((size_t)n, sizeof(double));
        if (!pairs[k].vector_im) pairs[k].vector_im = (double*)calloc((size_t)n, sizeof(double));
        for (int i = 0; i < n; i++) pairs[k].vector_re[i] = (i == k) ? 1.0 : 0.0;
        haken_inverse_iteration(A, n, re[k], im[k], pairs[k].vector_re,
                                 pairs[k].vector_im,
                                 config ? config->max_iterations : 20, 1e-10);
        /* Compute residual */
        double res = 0.0;
        for (int i = 0; i < n; i++) {
            double Av = 0.0;
            for (int j = 0; j < n; j++) Av += A[i * n + j] * pairs[k].vector_re[j];
            double diff = Av - re[k] * pairs[k].vector_re[i];
            res += diff * diff;
        }
        pairs[k].residual = sqrt(res);
    }
    free(work); free(re); free(im);
    return 0;
}
