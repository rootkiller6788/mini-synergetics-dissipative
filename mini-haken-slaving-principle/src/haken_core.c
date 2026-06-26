/* ==============================================================
 * haken_core.c — Core Implementation: System Lifecycle & Mode Analysis
 *
 * Implements the fundamental data structures and operations for
 * Haken's slaving principle in synergetics.
 *
 * References:
 *   Haken (1977, 1983) Synergetics
 *   Golub & Van Loan (2013) Matrix Computations, 4th ed.
 *   Press, Teukolsky, Vetterling, Flannery (2007) Numerical Recipes 3rd ed.
 * ============================================================== */

#include "haken_core.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==============================================================
 * Memory management helpers
 * ============================================================== */

static double* haken_alloc_vector(int n) {
    double* v = (double*)calloc((size_t)n, sizeof(double));
    if (!v) { fprintf(stderr, "haken: allocation failed for vector[%d]\n", n); exit(1); }
    return v;
}

static double* haken_alloc_matrix(int rows, int cols) {
    double* m = (double*)calloc((size_t)(rows * cols), sizeof(double));
    if (!m) { fprintf(stderr, "haken: allocation failed for matrix[%dx%d]\n", rows, cols); exit(1); }
    return m;
}

static void haken_free_vector(double* v) { if (v) free(v); }
static void haken_free_matrix(double* m) { if (m) free(m); }

#define MAT_AT(M, n, i, j) ((M)[(i) + (n)*(j)])

/* ==============================================================
 * L1: System lifecycle
 * ============================================================== */

Haken_System* haken_system_create(int n_dim,
    void (*L_func)(double alpha, int n, double* L_out),
    void (*N_func)(const double* q, int n, double* Nq))
{
    if (n_dim <= 0) return NULL;

    Haken_System* sys = (Haken_System*)calloc(1, sizeof(Haken_System));
    if (!sys) return NULL;

    sys->n_dim          = n_dim;
    sys->n_critical     = 0;
    sys->n_stable       = 0;
    sys->n_unstable     = 0;
    sys->control_param  = 0.0;
    sys->critical_param = 0.0;

    sys->state          = haken_alloc_vector(n_dim);
    sys->linear_matrix  = haken_alloc_matrix(n_dim, n_dim);
    sys->modes          = (Haken_Mode*)calloc((size_t)n_dim, sizeof(Haken_Mode));
    if (!sys->modes) { free(sys->state); free(sys->linear_matrix); free(sys); return NULL; }
    for (int k = 0; k < n_dim; k++) {
        sys->modes[k].eigenvector_re = haken_alloc_vector(n_dim);
        sys->modes[k].eigenvector_im = haken_alloc_vector(n_dim);
        sys->modes[k].left_eigenv_re = haken_alloc_vector(n_dim);
        sys->modes[k].left_eigenv_im = haken_alloc_vector(n_dim);
        sys->modes[k].index = k;
        sys->modes[k].type = HAKEN_MODE_STABLE;
    }

    sys->critical_indices = (int*)calloc((size_t)n_dim, sizeof(int));
    sys->stable_indices   = (int*)calloc((size_t)n_dim, sizeof(int));
    sys->unstable_indices = (int*)calloc((size_t)n_dim, sizeof(int));

    sys->linear_matrix_func = L_func;
    sys->nonlinear_coupling = N_func;
    sys->coupling_type = HAKEN_COUPLING_GENERAL;
    sys->context = NULL;

    return sys;
}

Haken_System* haken_system_clone(const Haken_System* src) {
    if (!src) return NULL;
    Haken_System* dst = haken_system_create(src->n_dim,
        src->linear_matrix_func, src->nonlinear_coupling);
    if (!dst) return NULL;

    dst->control_param  = src->control_param;
    dst->critical_param = src->critical_param;
    dst->coupling_type  = src->coupling_type;
    dst->context        = src->context;

    memcpy(dst->state, src->state, (size_t)(src->n_dim) * sizeof(double));
    memcpy(dst->linear_matrix, src->linear_matrix,
           (size_t)(src->n_dim * src->n_dim) * sizeof(double));

    for (int k = 0; k < src->n_dim; k++) {
        dst->modes[k].lambda_re  = src->modes[k].lambda_re;
        dst->modes[k].lambda_im  = src->modes[k].lambda_im;
        dst->modes[k].type       = src->modes[k].type;
        dst->modes[k].relaxation_time = src->modes[k].relaxation_time;
        dst->modes[k].oscillation_period = src->modes[k].oscillation_period;
        memcpy(dst->modes[k].eigenvector_re, src->modes[k].eigenvector_re,
               (size_t)(src->n_dim) * sizeof(double));
        memcpy(dst->modes[k].eigenvector_im, src->modes[k].eigenvector_im,
               (size_t)(src->n_dim) * sizeof(double));
    }
    return dst;
}

void haken_system_free(Haken_System* sys) {
    if (!sys) return;
    if (sys->modes) {
        for (int k = 0; k < sys->n_dim; k++) {
            haken_free_vector(sys->modes[k].eigenvector_re);
            haken_free_vector(sys->modes[k].eigenvector_im);
            haken_free_vector(sys->modes[k].left_eigenv_re);
            haken_free_vector(sys->modes[k].left_eigenv_im);
        }
        free(sys->modes);
    }
    haken_free_vector(sys->state);
    haken_free_matrix(sys->linear_matrix);
    free(sys->critical_indices);
    free(sys->stable_indices);
    free(sys->unstable_indices);
    free(sys);
}

void haken_set_control_param(Haken_System* sys, double alpha) {
    if (!sys) return;
    sys->control_param = alpha;
}

void haken_set_critical_param(Haken_System* sys, double alpha_c) {
    if (!sys) return;
    sys->critical_param = alpha_c;
}

void haken_set_state(Haken_System* sys, const double* q) {
    if (!sys || !q) return;
    memcpy(sys->state, q, (size_t)(sys->n_dim) * sizeof(double));
}

void haken_update_linear_matrix(Haken_System* sys) {
    if (!sys) return;
    if (sys->linear_matrix_func) {
        sys->linear_matrix_func(sys->control_param, sys->n_dim, sys->linear_matrix);
    }
}

/* ==============================================================
 * L3: QR algorithm for eigenvalue computation
 *
 * Implements the implicit double-shift QR algorithm for real
 * matrices, as described in Golub & Van Loan §7.5.
 *
 * The algorithm reduces A to upper Hessenberg form, then iterates
 * QR steps with Wilkinson shifts until convergence to quasi-triangular
 * (real Schur) form, from which eigenvalues are extracted.
 * ============================================================== */

/* Householder reflection: P = I - β·v·v^T */
static double householder_vec(double* x, int n) {
    double sigma = 0.0;
    for (int i = 1; i < n; i++) sigma += x[i] * x[i];
    double x0 = x[0];
    double beta;
    if (sigma == 0.0 && x0 >= 0) { beta = 0.0; return beta; }
    double mu = sqrt(x0 * x0 + sigma);
    double v0;
    if (x0 <= 0) { v0 = x0 - mu; } else { v0 = -sigma / (x0 + mu); }
    beta = 2.0 * v0 * v0 / (sigma + v0 * v0);
    double inv_v0 = 1.0 / v0;
    x[0] = v0;
    for (int i = 1; i < n; i++) x[i] *= inv_v0;
    return beta;
}

static void apply_householder_left(double* A, int n, const double* v, double beta, int col) {
    if (beta == 0.0) return;
    for (int j = col; j < n; j++) {
        double s = 0.0;
        for (int i = col; i < n; i++) s += v[i - col] * A[i * n + j];
        s *= beta;
        for (int i = col; i < n; i++) A[i * n + j] -= s * v[i - col];
    }
}

static void apply_householder_right(double* A, int n, const double* v, double beta, int row_end) {
    if (beta == 0.0) return;
    for (int i = 0; i <= row_end; i++) {
        double s = 0.0;
        for (int j = 0; j < n; j++) s += A[i * n + j] * v[j];
        s *= beta;
        for (int j = 0; j < n; j++) A[i * n + j] -= s * v[j];
    }
}

/* Reduce to upper Hessenberg form (destroying A, building Q) */
static void to_hessenberg(double* A, int n, double* Q) {
    double* v = (double*)calloc((size_t)n, sizeof(double));
    for (int j = 0; j < n - 2; j++) {
        for (int i = j + 1; i < n; i++) v[i - j - 1] = A[i * n + j];
        double beta = householder_vec(v, n - j - 1);
        if (beta > 0.0) {
            apply_householder_left(A + j * n, n, v, beta, j);
            apply_householder_right(A, n, v, beta, n - 1);
        }
        /* Accumulate Q */
        for (int i = 0; i < n; i++) {
            double s = 0.0;
            for (int k = j + 1; k < n; k++) s += Q[i * n + k] * v[k - j - 1];
            s *= beta;
            for (int k = j + 1; k < n; k++) Q[i * n + k] -= s * v[k - j - 1];
        }
    }
    free(v);
}

/* Wilkinson shift for a 2×2 trailing submatrix */
static void wilkinson_shift(double a, double b, double c, double d,
                             double* sr, double* si) {
    double trace = a + d;
    double det   = a * d - b * c;
    double disc  = trace * trace - 4.0 * det;
    if (disc >= 0) {
        double sd = sqrt(disc);
        if (trace >= 0) { *sr = (trace + sd) / 2.0; *si = 0.0; }
        else            { *sr = (trace - sd) / 2.0; *si = 0.0; }
    } else {
        *sr = trace / 2.0;
        *si = sqrt(-disc) / 2.0;
    }
}

/* Givens rotation: eliminate A[i-1, i] */
static void givens_rotation(double x, double y, double* c, double* s) {
    if (fabs(y) < 1e-15) { *c = 1.0; *s = 0.0; return; }
    if (fabs(x) > fabs(y)) {
        double t = y / x;
        double r = sqrt(1.0 + t * t);
        *c = 1.0 / r;
        *s = t / r;
    } else {
        double t = x / y;
        double r = sqrt(1.0 + t * t);
        *s = 1.0 / r;
        *c = t / r;
    }
}

/* Extract eigenvalues from quasi-triangular matrix.
 * For 1×1 diagonal blocks: eigenvalue is the diagonal entry.
 * For 2×2 diagonal blocks: eigenvalues are roots of char poly. */
static int extract_eigenvalues(double* H, int n, double* re, double* im) {
    int count = 0;
    int i = 0;
    while (i < n) {
        if (i == n - 1 || fabs(H[(i + 1) * n + i]) < 1e-14) {
            re[count] = H[i * n + i];
            im[count] = 0.0;
            count++; i++;
        } else {
            double a = H[i * n + i],    b = H[i * n + (i + 1)];
            double c_val = H[(i + 1) * n + i], d = H[(i + 1) * n + (i + 1)];
            double tr = a + d, det = a * d - b * c_val;
            double disc = tr * tr - 4.0 * det;
            if (disc >= 0) {
                double sd = sqrt(disc);
                re[count] = (tr + sd) / 2.0; im[count] = 0.0; count++;
                re[count] = (tr - sd) / 2.0; im[count] = 0.0; count++;
            } else {
                re[count] = tr / 2.0; im[count] = sqrt(-disc) / 2.0; count++;
                re[count] = tr / 2.0; im[count] = -sqrt(-disc) / 2.0; count++;
            }
            i += 2;
        }
    }
    return count;
}

/* -------------------------------------------------------------
 * Public eigenvalue API
 * ------------------------------------------------------------- */

int haken_compute_eigenmodes(Haken_System* sys) {
    if (!sys) return -1;
    int n = sys->n_dim;
    if (n <= 0) return -1;

    haken_update_linear_matrix(sys);

    /* Allocate working copy */
    double* A = haken_alloc_matrix(n, n);
    memcpy(A, sys->linear_matrix, (size_t)(n * n) * sizeof(double));

    /* Initialize Q as identity */
    double* Q = haken_alloc_matrix(n, n);
    for (int i = 0; i < n; i++) Q[i * n + i] = 1.0;

    /* Reduce to Hessenberg */
    to_hessenberg(A, n, Q);

    /* QR iterations */
    double* v = (double*)calloc((size_t)n, sizeof(double));
    int max_iter = 100 * n;
    for (int iter = 0; iter < max_iter; iter++) {
        /* Check deflation: look for negligible subdiagonal entries */
        int m = n;
        while (m > 1) {
            if (fabs(A[(m - 1) * n + (m - 2)]) < 1e-15) {
                m--;
            } else if (m > 2 && fabs(A[(m - 2) * n + (m - 3)]) < 1e-15) {
                m--;
            } else {
                break;
            }
        }
        if (m <= 1) break;

        /* Wilkinson shift from 2×2 trailing block */
        int k = m - 1;
        double sr, si;
        wilkinson_shift(A[(k - 1) * n + (k - 1)], A[(k - 1) * n + k],
                        A[k * n + (k - 1)], A[k * n + k], &sr, &si);

        /* QR step with implicit shift */
        double x = A[0] - sr;
        double y = A[1 * n];
        for (int j = 0; j < m - 1; j++) {
            double c, s;
            givens_rotation(x, y, &c, &s);
            /* Apply Givens to rows j and j+1 */
            for (int col = 0; col < m; col++) {
                double t1 = c * A[j * n + col] + s * A[(j + 1) * n + col];
                double t2 = -s * A[j * n + col] + c * A[(j + 1) * n + col];
                A[j * n + col] = t1;
                A[(j + 1) * n + col] = t2;
            }
            /* Apply Givens transpose to columns j and j+1 */
            for (int row = 0; row < m; row++) {
                double t1 = c * A[row * n + j] + s * A[row * n + (j + 1)];
                double t2 = -s * A[row * n + j] + c * A[row * n + (j + 1)];
                A[row * n + j] = t1;
                A[row * n + (j + 1)] = t2;
            }
            /* Accumulate Q */
            for (int row = 0; row < n; row++) {
                double t1 = c * Q[row * n + j] + s * Q[row * n + (j + 1)];
                double t2 = -s * Q[row * n + j] + c * Q[row * n + (j + 1)];
                Q[row * n + j] = t1;
                Q[row * n + (j + 1)] = t2;
            }
            if (j < m - 2) {
                x = A[(j + 1) * n + j];
                y = A[(j + 2) * n + j];
            }
        }
    }

    /* Extract eigenvalues */
    double* ev_re = haken_alloc_vector(n);
    double* ev_im = haken_alloc_vector(n);
    extract_eigenvalues(A, n, ev_re, ev_im);

    /* Sort by real part descending (slowest first) */
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            if (ev_re[j] > ev_re[i]) {
                double tmp = ev_re[i]; ev_re[i] = ev_re[j]; ev_re[j] = tmp;
                tmp = ev_im[i]; ev_im[i] = ev_im[j]; ev_im[j] = tmp;
                /* Also swap eigenvectors in Q */
                for (int r = 0; r < n; r++) {
                    tmp = Q[r * n + i]; Q[r * n + i] = Q[r * n + j]; Q[r * n + j] = tmp;
                }
            }
        }
    }

    /* Populate modes */
    for (int k = 0; k < n; k++) {
        sys->modes[k].lambda_re = ev_re[k];
        sys->modes[k].lambda_im = ev_im[k];
        sys->modes[k].relaxation_time = (ev_re[k] < 0) ? (-1.0 / ev_re[k]) : INFINITY;
        sys->modes[k].oscillation_period = (fabs(ev_im[k]) > 1e-15)
            ? (2.0 * M_PI / fabs(ev_im[k])) : INFINITY;
        for (int i = 0; i < n; i++) {
            sys->modes[k].eigenvector_re[i] = Q[i * n + k];
            sys->modes[k].eigenvector_im[i] = 0.0;  /* Real eigenvectors from QR */
        }
    }

    haken_free_matrix(A);
    haken_free_matrix(Q);
    haken_free_vector(ev_re);
    haken_free_vector(ev_im);
    free(v);
    return 0;
}

/* ==============================================================
 * Mode classification
 * ============================================================== */

void haken_classify_modes(Haken_System* sys, double eps_crit, double eps_stable) {
    if (!sys) return;
    int n = sys->n_dim;
    sys->n_critical = 0;
    sys->n_stable   = 0;
    sys->n_unstable = 0;

    for (int k = 0; k < n; k++) {
        double r = sys->modes[k].lambda_re;
        if (r > eps_crit) {
            /* Growing mode → unstable */
            sys->modes[k].type = HAKEN_MODE_UNSTABLE;
            sys->unstable_indices[sys->n_unstable++] = k;
        } else if (fabs(r) <= eps_crit) {
            /* Near zero → critical (order parameter candidate) */
            sys->modes[k].type = HAKEN_MODE_CRITICAL;
            sys->critical_indices[sys->n_critical++] = k;
        } else if (r < -eps_stable) {
            /* Strongly damped → stable / slaved */
            sys->modes[k].type = HAKEN_MODE_STABLE;
            sys->stable_indices[sys->n_stable++] = k;
        } else {
            /* Between critical and stable → marginal */
            sys->modes[k].type = HAKEN_MODE_MARGINAL;
            sys->stable_indices[sys->n_stable++] = k;
        }
    }
}

/* ==============================================================
 * Subspace projection
 * ============================================================== */

void haken_project_critical(const Haken_System* sys, const double* q,
                             double* q_s) {
    if (!sys || !q || !q_s) return;
    int nc = sys->n_critical;
    for (int a = 0; a < nc; a++) {
        int k = sys->critical_indices[a];
        double proj = 0.0;
        for (int i = 0; i < sys->n_dim; i++) {
            /* Project using left eigenvector · state (scalar product) */
            proj += sys->modes[k].left_eigenv_re[i] * q[i];
        }
        q_s[a] = proj;
    }
}

void haken_project_stable(const Haken_System* sys, const double* q,
                           double* q_f) {
    if (!sys || !q || !q_f) return;
    int ns = sys->n_stable;
    for (int a = 0; a < ns; a++) {
        int k = sys->stable_indices[a];
        double proj = 0.0;
        for (int i = 0; i < sys->n_dim; i++) {
            proj += sys->modes[k].left_eigenv_re[i] * q[i];
        }
        q_f[a] = proj;
    }
}

/* ==============================================================
 * Utility: compute the spectral gap of the system
 * ============================================================== */

double haken_spectral_gap(const Haken_System* sys) {
    if (!sys || sys->n_critical == 0 || sys->n_stable == 0) return 0.0;
    double max_crit_re = -INFINITY;
    for (int a = 0; a < sys->n_critical; a++) {
        double r = sys->modes[sys->critical_indices[a]].lambda_re;
        if (r > max_crit_re) max_crit_re = r;
    }
    double max_stable_re = -INFINITY;
    for (int a = 0; a < sys->n_stable; a++) {
        double r = sys->modes[sys->stable_indices[a]].lambda_re;
        if (r > max_stable_re) max_stable_re = r;
    }
    return max_crit_re - max_stable_re;
}
