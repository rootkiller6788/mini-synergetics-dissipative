#include "oep_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==============================================================
 * Core Linear Algebra and Utility Implementation
 *
 * Reference: Golub & Van Loan (2013) "Matrix Computations"
 * All matrix operations are O(n³) worst case, optimized for
 * the small-to-medium matrices (n ≤ OEP_MAX_DIM) typical in
 * order parameter analysis.
 * ============================================================== */

/* ── RNG State ───────────────────────────────────────── */
static unsigned int oep_rng_state = 123456789;

/* ── Vector Lifecycle ────────────────────────────────── */

OEPVector* oep_vec_create(int size) {
    if (size <= 0 || size > OEP_MAX_DIM) return NULL;
    OEPVector* v = calloc(1, sizeof(OEPVector));
    if (!v) return NULL;
    v->size = size;
    v->capacity = size;
    v->data = calloc((size_t)size, sizeof(double));
    if (!v->data) { free(v); return NULL; }
    return v;
}

OEPVector* oep_vec_create_from(const double* data, int size) {
    OEPVector* v = oep_vec_create(size);
    if (v && data) memcpy(v->data, data, (size_t)size * sizeof(double));
    return v;
}

OEPVector* oep_vec_clone(const OEPVector* v) {
    if (!v) return NULL;
    return oep_vec_create_from(v->data, v->size);
}

void oep_vec_free(OEPVector* v) {
    if (!v) return;
    free(v->data);
    free(v);
}

void oep_vec_set(OEPVector* v, int i, double val) {
    if (v && i >= 0 && i < v->size) v->data[i] = val;
}

double oep_vec_get(const OEPVector* v, int i) {
    if (!v || i < 0 || i >= v->size) return 0.0;
    return v->data[i];
}

void oep_vec_set_all(OEPVector* v, const double* data) {
    if (v && data) memcpy(v->data, data, (size_t)v->size * sizeof(double));
}

void oep_vec_zero(OEPVector* v) {
    if (v) memset(v->data, 0, (size_t)v->size * sizeof(double));
}

void oep_vec_copy(OEPVector* dst, const OEPVector* src) {
    if (!dst || !src) return;
    int n = (dst->size < src->size) ? dst->size : src->size;
    memcpy(dst->data, src->data, (size_t)n * sizeof(double));
}

void oep_vec_print(const OEPVector* v, const char* name) {
    if (!v) return;
    printf("%s [%d]:", name ? name : "vec", v->size);
    for (int i = 0; i < v->size; i++) {
        if (i % 8 == 0) printf("\n  ");
        printf("%10.6f ", v->data[i]);
    }
    printf("\n");
}

/* ── Vector Operations ───────────────────────────────── */

double oep_vec_dot(const OEPVector* a, const OEPVector* b) {
    if (!a || !b) return 0.0;
    int n = (a->size < b->size) ? a->size : b->size;
    double s = 0.0;
    for (int i = 0; i < n; i++) s += a->data[i] * b->data[i];
    return s;
}

double oep_vec_norm(const OEPVector* v) {
    return v ? sqrt(oep_vec_dot(v, v)) : 0.0;
}

double oep_vec_norm_max(const OEPVector* v) {
    if (!v) return 0.0;
    double m = 0.0;
    for (int i = 0; i < v->size; i++) {
        double a = fabs(v->data[i]);
        if (a > m) m = a;
    }
    return m;
}

OEPVector* oep_vec_add(const OEPVector* a, const OEPVector* b) {
    if (!a || !b || a->size != b->size) return NULL;
    OEPVector* r = oep_vec_create(a->size);
    if (!r) return NULL;
    for (int i = 0; i < a->size; i++) r->data[i] = a->data[i] + b->data[i];
    return r;
}

OEPVector* oep_vec_sub(const OEPVector* a, const OEPVector* b) {
    if (!a || !b || a->size != b->size) return NULL;
    OEPVector* r = oep_vec_create(a->size);
    if (!r) return NULL;
    for (int i = 0; i < a->size; i++) r->data[i] = a->data[i] - b->data[i];
    return r;
}

OEPVector* oep_vec_scale(const OEPVector* v, double s) {
    if (!v) return NULL;
    OEPVector* r = oep_vec_create(v->size);
    if (!r) return NULL;
    for (int i = 0; i < v->size; i++) r->data[i] = v->data[i] * s;
    return r;
}

OEPVector* oep_vec_linear_combination(double a, const OEPVector* va,
                                       double b, const OEPVector* vb) {
    if (!va || !vb || va->size != vb->size) return NULL;
    OEPVector* r = oep_vec_create(va->size);
    if (!r) return NULL;
    for (int i = 0; i < va->size; i++)
        r->data[i] = a * va->data[i] + b * vb->data[i];
    return r;
}

double oep_vec_mean(const OEPVector* v) {
    if (!v || v->size == 0) return 0.0;
    double s = 0.0;
    for (int i = 0; i < v->size; i++) s += v->data[i];
    return s / (double)v->size;
}

double oep_vec_variance(const OEPVector* v) {
    if (!v || v->size < 2) return 0.0;
    double mu = oep_vec_mean(v);
    double s = 0.0;
    for (int i = 0; i < v->size; i++) {
        double d = v->data[i] - mu;
        s += d * d;
    }
    return s / (double)(v->size - 1);
}

double oep_vec_correlation(const OEPVector* a, const OEPVector* b) {
    if (!a || !b || a->size < 2 || a->size != b->size) return 0.0;
    int n = a->size;
    double ma = oep_vec_mean(a), mb = oep_vec_mean(b);
    double sa = 0.0, sb = 0.0, sab = 0.0;
    for (int i = 0; i < n; i++) {
        double da = a->data[i] - ma;
        double db = b->data[i] - mb;
        sa += da * da;
        sb += db * db;
        sab += da * db;
    }
    if (sa < OEP_EPSILON || sb < OEP_EPSILON) return 0.0;
    return sab / sqrt(sa * sb);
}

/* ── Matrix Lifecycle ────────────────────────────────── */

OEPMatrix* oep_mat_create(int rows, int cols) {
    if (rows <= 0 || cols <= 0 || rows > OEP_MAX_DIM || cols > OEP_MAX_DIM)
        return NULL;
    OEPMatrix* m = calloc(1, sizeof(OEPMatrix));
    if (!m) return NULL;
    m->rows = rows;
    m->cols = cols;
    m->capacity = rows * cols;
    m->data = calloc((size_t)m->capacity, sizeof(double));
    if (!m->data) { free(m); return NULL; }
    return m;
}

OEPMatrix* oep_mat_create_from(const double* data, int rows, int cols) {
    OEPMatrix* m = oep_mat_create(rows, cols);
    if (m && data) memcpy(m->data, data, (size_t)(rows * cols) * sizeof(double));
    return m;
}

OEPMatrix* oep_mat_clone(const OEPMatrix* m) {
    if (!m) return NULL;
    return oep_mat_create_from(m->data, m->rows, m->cols);
}

void oep_mat_free(OEPMatrix* m) {
    if (!m) return;
    free(m->data);
    free(m);
}

void oep_mat_set(OEPMatrix* m, int i, int j, double val) {
    if (m && i >= 0 && i < m->rows && j >= 0 && j < m->cols)
        m->data[i * m->cols + j] = val;
}

double oep_mat_get(const OEPMatrix* m, int i, int j) {
    if (!m || i < 0 || i >= m->rows || j < 0 || j >= m->cols) return 0.0;
    return m->data[i * m->cols + j];
}

void oep_mat_set_all(OEPMatrix* m, const double* data) {
    if (m && data) memcpy(m->data, data, (size_t)(m->rows * m->cols) * sizeof(double));
}

void oep_mat_identity(OEPMatrix* m) {
    if (!m) return;
    oep_mat_zero(m);
    int n = (m->rows < m->cols) ? m->rows : m->cols;
    for (int i = 0; i < n; i++) m->data[i * m->cols + i] = 1.0;
}

void oep_mat_zero(OEPMatrix* m) {
    if (m) memset(m->data, 0, (size_t)(m->rows * m->cols) * sizeof(double));
}

void oep_mat_copy(OEPMatrix* dst, const OEPMatrix* src) {
    if (!dst || !src) return;
    int n = (dst->rows * dst->cols < src->rows * src->cols)
            ? dst->rows * dst->cols : src->rows * src->cols;
    memcpy(dst->data, src->data, (size_t)n * sizeof(double));
}

void oep_mat_print(const OEPMatrix* m, const char* name) {
    if (!m) return;
    printf("%s [%dx%d]:\n", name ? name : "mat", m->rows, m->cols);
    for (int i = 0; i < m->rows; i++) {
        printf("  ");
        for (int j = 0; j < m->cols; j++)
            printf("%10.4f ", m->data[i * m->cols + j]);
        printf("\n");
    }
}

/* ── Matrix Operations ───────────────────────────────── */

OEPMatrix* oep_mat_add(const OEPMatrix* a, const OEPMatrix* b) {
    if (!a || !b || a->rows != b->rows || a->cols != b->cols) return NULL;
    OEPMatrix* r = oep_mat_create(a->rows, a->cols);
    if (!r) return NULL;
    int n = a->rows * a->cols;
    for (int i = 0; i < n; i++) r->data[i] = a->data[i] + b->data[i];
    return r;
}

OEPMatrix* oep_mat_sub(const OEPMatrix* a, const OEPMatrix* b) {
    if (!a || !b || a->rows != b->rows || a->cols != b->cols) return NULL;
    OEPMatrix* r = oep_mat_create(a->rows, a->cols);
    if (!r) return NULL;
    int n = a->rows * a->cols;
    for (int i = 0; i < n; i++) r->data[i] = a->data[i] - b->data[i];
    return r;
}

OEPMatrix* oep_mat_mul(const OEPMatrix* a, const OEPMatrix* b) {
    if (!a || !b || a->cols != b->rows) return NULL;
    OEPMatrix* r = oep_mat_create(a->rows, b->cols);
    if (!r) return NULL;
    for (int i = 0; i < a->rows; i++) {
        for (int j = 0; j < b->cols; j++) {
            double s = 0.0;
            for (int k = 0; k < a->cols; k++)
                s += a->data[i * a->cols + k] * b->data[k * b->cols + j];
            r->data[i * r->cols + j] = s;
        }
    }
    return r;
}

OEPMatrix* oep_mat_scale(const OEPMatrix* m, double s) {
    if (!m) return NULL;
    OEPMatrix* r = oep_mat_create(m->rows, m->cols);
    if (!r) return NULL;
    int n = m->rows * m->cols;
    for (int i = 0; i < n; i++) r->data[i] = m->data[i] * s;
    return r;
}

OEPMatrix* oep_mat_transpose(const OEPMatrix* m) {
    if (!m) return NULL;
    OEPMatrix* r = oep_mat_create(m->cols, m->rows);
    if (!r) return NULL;
    for (int i = 0; i < m->rows; i++)
        for (int j = 0; j < m->cols; j++)
            r->data[j * r->cols + i] = m->data[i * m->cols + j];
    return r;
}

OEPVector* oep_mat_vec_mul(const OEPMatrix* A, const OEPVector* x) {
    if (!A || !x || A->cols != x->size) return NULL;
    OEPVector* r = oep_vec_create(A->rows);
    if (!r) return NULL;
    for (int i = 0; i < A->rows; i++) {
        double s = 0.0;
        for (int j = 0; j < A->cols; j++)
            s += A->data[i * A->cols + j] * x->data[j];
        r->data[i] = s;
    }
    return r;
}

double oep_mat_determinant(const OEPMatrix* m) {
    if (!m || m->rows != m->cols || m->rows == 0) return 0.0;
    /* Gaussian elimination with partial pivoting for determinant.
     * Complexity: O(n³), Reference: Golub & Van Loan §3.4 */
    int n = m->rows;
    double* a = malloc((size_t)(n * n) * sizeof(double));
    if (!a) return 0.0;
    memcpy(a, m->data, (size_t)(n * n) * sizeof(double));

    double det = 1.0;
    for (int k = 0; k < n; k++) {
        /* Partial pivoting */
        int pivot = k;
        double max_val = fabs(a[k * n + k]);
        for (int i = k + 1; i < n; i++) {
            double v = fabs(a[i * n + k]);
            if (v > max_val) { max_val = v; pivot = i; }
        }
        if (max_val < OEP_EPSILON) { free(a); return 0.0; }
        if (pivot != k) {
            det = -det;
            for (int j = k; j < n; j++) {
                double t = a[k * n + j];
                a[k * n + j] = a[pivot * n + j];
                a[pivot * n + j] = t;
            }
        }
        det *= a[k * n + k];
        for (int i = k + 1; i < n; i++) {
            double factor = a[i * n + k] / a[k * n + k];
            for (int j = k + 1; j < n; j++)
                a[i * n + j] -= factor * a[k * n + j];
        }
    }
    free(a);
    return det;
}

double oep_mat_trace(const OEPMatrix* m) {
    if (!m) return 0.0;
    int n = (m->rows < m->cols) ? m->rows : m->cols;
    double t = 0.0;
    for (int i = 0; i < n; i++) t += m->data[i * m->cols + i];
    return t;
}

double oep_mat_norm_frobenius(const OEPMatrix* m) {
    if (!m) return 0.0;
    double s = 0.0;
    int n = m->rows * m->cols;
    for (int i = 0; i < n; i++) s += m->data[i] * m->data[i];
    return sqrt(s);
}

double oep_mat_norm_inf(const OEPMatrix* m) {
    if (!m) return 0.0;
    double max_sum = 0.0;
    for (int i = 0; i < m->rows; i++) {
        double row_sum = 0.0;
        for (int j = 0; j < m->cols; j++)
            row_sum += fabs(m->data[i * m->cols + j]);
        if (row_sum > max_sum) max_sum = row_sum;
    }
    return max_sum;
}

int oep_mat_rank_gaussian(const OEPMatrix* m) {
    if (!m) return 0;
    /* Compute rank via Gaussian elimination with tolerance.
     * Reference: Golub & Van Loan (2013) §5.4 */
    int r = m->rows, c = m->cols;
    if (r == 0 || c == 0) return 0;
    int n_max = (r * c > OEP_MAX_DIM * OEP_MAX_DIM) ? OEP_MAX_DIM * OEP_MAX_DIM : r * c;
    double* a = malloc((size_t)n_max * sizeof(double));
    if (!a) return 0;
    memcpy(a, m->data, (size_t)(r * c) * sizeof(double));

    int rank = 0;
    int col = 0;
    for (int row = 0; row < r && col < c; col++) {
        /* Find pivot */
        int pivot = -1;
        double max_val = OEP_EPSILON;
        for (int i = row; i < r; i++) {
            double v = fabs(a[i * c + col]);
            if (v > max_val) { max_val = v; pivot = i; }
        }
        if (pivot < 0) continue;
        /* Swap rows */
        if (pivot != row) {
            for (int j = col; j < c; j++) {
                double t = a[row * c + j];
                a[row * c + j] = a[pivot * c + j];
                a[pivot * c + j] = t;
            }
        }
        /* Eliminate below */
        double piv = a[row * c + col];
        for (int i = row + 1; i < r; i++) {
            double factor = a[i * c + col] / piv;
            for (int j = col; j < c; j++)
                a[i * c + j] -= factor * a[row * c + j];
        }
        rank++;
        row++;
    }
    free(a);
    return rank;
}

/* ── ODE Solver (Classical Runge-Kutta 4) ─────────────── */

double oep_ode_rk4_step(double t, const double* y, int n, double dt,
                         void (*f)(double, const double*, double*, void*),
                         void* ctx, double* y_out) {
    /* Classical RK4: k1 = dt*f(t, y), k2 = dt*f(t+dt/2, y+k1/2),
     *   k3 = dt*f(t+dt/2, y+k2/2), k4 = dt*f(t+dt, y+k3).
     * Then y_{n+1} = y_n + (k1 + 2k2 + 2k3 + k4)/6.
     * Complexity: 4 evaluations of f per step. */
    if (n <= 0 || dt <= 0.0 || !y || !f || !y_out) return t;

    double* k1 = malloc((size_t)(4 * n) * sizeof(double));
    if (!k1) { memcpy(y_out, y, (size_t)n * sizeof(double)); return t + dt; }
    double* k2 = k1 + n;
    double* k3 = k2 + n;
    double* k4 = k3 + n;
    double* ytmp = malloc((size_t)n * sizeof(double));

    /* k1 = f(t, y) */
    f(t, y, k1, ctx);

    /* k2 = f(t + dt/2, y + k1*dt/2) */
    for (int i = 0; i < n; i++) ytmp[i] = y[i] + 0.5 * dt * k1[i];
    f(t + 0.5 * dt, ytmp, k2, ctx);

    /* k3 = f(t + dt/2, y + k2*dt/2) */
    for (int i = 0; i < n; i++) ytmp[i] = y[i] + 0.5 * dt * k2[i];
    f(t + 0.5 * dt, ytmp, k3, ctx);

    /* k4 = f(t + dt, y + k3*dt) */
    for (int i = 0; i < n; i++) ytmp[i] = y[i] + dt * k3[i];
    f(t + dt, ytmp, k4, ctx);

    for (int i = 0; i < n; i++)
        y_out[i] = y[i] + dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;

    free(ytmp);
    free(k1);
    return t + dt;
}

void oep_ode_rk4_multi(double t0, const double* y0, int n, int steps,
                        double dt,
                        void (*f)(double, const double*, double*, void*),
                        void* ctx, double** y_trajectory) {
    /* Multi-step RK4 integration. Stores trajectory in y_trajectory[n_steps*n].
     * Each row of y_trajectory is one time point. */
    if (!y0 || !f || !y_trajectory || n <= 0 || steps <= 0) return;

    double* y_cur = malloc((size_t)n * sizeof(double));
    double* y_next = malloc((size_t)n * sizeof(double));
    if (!y_cur || !y_next) { free(y_cur); free(y_next); return; }
    memcpy(y_cur, y0, (size_t)n * sizeof(double));

    /* Store initial state */
    memcpy(y_trajectory[0], y_cur, (size_t)n * sizeof(double));

    double t = t0;
    for (int step = 0; step < steps - 1; step++) {
        t = oep_ode_rk4_step(t, y_cur, n, dt, f, ctx, y_next);
        memcpy(y_trajectory[step + 1], y_next, (size_t)n * sizeof(double));
        double* tmp = y_cur;
        y_cur = y_next;
        y_next = tmp;
    }
    free(y_cur);
    free(y_next);
}

/* ── Utility Functions ─────────────────────────────────── */

double oep_max(double a, double b) { return (a > b) ? a : b; }
double oep_min(double a, double b) { return (a < b) ? a : b; }

double oep_clamp(double x, double lo, double hi) {
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

double oep_sigmoid(double x) { return 1.0 / (1.0 + exp(-x)); }

double oep_tanh_smooth(double x, double scale) {
    return tanh(scale * x);
}

double oep_gaussian(double x, double mu, double sigma) {
    if (sigma <= 0.0) return (fabs(x - mu) < OEP_EPSILON) ? 1.0 : 0.0;
    double z = (x - mu) / sigma;
    return exp(-0.5 * z * z) / (sigma * sqrt(2.0 * OEP_PI));
}

bool oep_is_zero(double x) { return fabs(x) < OEP_EPSILON; }
bool oep_is_equal(double a, double b) { return fabs(a - b) < OEP_EPSILON; }
double oep_sign(double x) { return (x > 0.0) ? 1.0 : ((x < 0.0) ? -1.0 : 0.0); }

double oep_heaviside(double x, double threshold) {
    return (x >= threshold) ? 1.0 : 0.0;
}

double oep_deg2rad(double deg) { return deg * OEP_PI / 180.0; }
double oep_rad2deg(double rad) { return rad * 180.0 / OEP_PI; }

/* ── Random Number Generation ──────────────────────────── */
/* Using a simple LCG for reproducibility (not cryptographic).
 * Reference: Park & Miller (1988) "Minimal Standard" LCG. */

void oep_rng_set_seed(unsigned int seed) { oep_rng_state = seed; }

double oep_rng_uniform(void) {
    oep_rng_state = oep_rng_state * 1103515245U + 12345U;
    return (double)(oep_rng_state & 0x7FFFFFFF) / (double)0x7FFFFFFF;
}

double oep_rng_normal(double mu, double sigma) {
    /* Box-Muller transform */
    double u1 = oep_rng_uniform();
    double u2 = oep_rng_uniform();
    /* Guard against u1 == 0 */
    if (u1 < 1e-10) u1 = 1e-10;
    double z = sqrt(-2.0 * log(u1)) * cos(2.0 * OEP_PI * u2);
    return mu + sigma * z;
}

double oep_rng_exponential(double lambda) {
    if (lambda <= 0.0) return 0.0;
    double u = oep_rng_uniform();
    return -log(u) / lambda;
}

void oep_rng_normal_vector(double* out, int n, double mu, double sigma) {
    if (!out) return;
    for (int i = 0; i < n; i++) out[i] = oep_rng_normal(mu, sigma);
}

/* ── Timeseries Utilities ──────────────────────────────── */

double oep_ts_mean(const double* data, int n) {
    if (!data || n <= 0) return 0.0;
    double s = 0.0;
    for (int i = 0; i < n; i++) s += data[i];
    return s / (double)n;
}

double oep_ts_variance(const double* data, int n) {
    if (!data || n < 2) return 0.0;
    double mu = oep_ts_mean(data, n);
    double s = 0.0;
    for (int i = 0; i < n; i++) { double d = data[i] - mu; s += d * d; }
    return s / (double)(n - 1);
}

double oep_ts_autocorrelation(const double* data, int n, int lag) {
    if (!data || n < 2 || lag < 0 || lag >= n) return 0.0;
    double mu = oep_ts_mean(data, n);
    double num = 0.0, den = 0.0;
    for (int i = 0; i < n; i++) {
        double d = data[i] - mu;
        den += d * d;
    }
    if (den < OEP_EPSILON) return 0.0;
    for (int i = 0; i < n - lag; i++) {
        num += (data[i] - mu) * (data[i + lag] - mu);
    }
    return num / den;
}

double oep_ts_cross_correlation(const double* x, const double* y, int n, int lag) {
    if (!x || !y || n < 2 || lag < 0 || lag >= n) return 0.0;
    double mx = oep_ts_mean(x, n), my = oep_ts_mean(y, n);
    double sxx = 0.0, syy = 0.0;
    for (int i = 0; i < n; i++) {
        sxx += (x[i] - mx) * (x[i] - mx);
        syy += (y[i] - my) * (y[i] - my);
    }
    if (sxx < OEP_EPSILON || syy < OEP_EPSILON) return 0.0;
    double sxy = 0.0;
    for (int i = 0; i < n - lag; i++)
        sxy += (x[i] - mx) * (y[i + lag] - my);
    return sxy / sqrt(sxx * syy);
}

void oep_ts_detrend_linear(double* data, int n) {
    if (!data || n < 3) return;
    /* Fit y = a + b*t via OLS, subtract trend */
    double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
    for (int i = 0; i < n; i++) {
        double t = (double)i;
        sx += t;
        sy += data[i];
        sxx += t * t;
        sxy += t * data[i];
    }
    double b = (n * sxy - sx * sy) / (n * sxx - sx * sx);
    double a = (sy - b * sx) / (double)n;
    for (int i = 0; i < n; i++) data[i] -= (a + b * (double)i);
}

void oep_ts_smooth_moving_average(const double* in, double* out, int n, int window) {
    if (!in || !out || n <= 0 || window <= 0) return;
    if (window > n) window = n;
    int half = window / 2;
    for (int i = 0; i < n; i++) {
        int start = i - half;
        int end = i + half;
        if (start < 0) { end -= start; start = 0; }
        if (end >= n) { start -= (end - (n - 1)); end = n - 1; }
        if (start < 0) start = 0;
        double sum = 0.0;
        int count = 0;
        for (int j = start; j <= end; j++) { sum += in[j]; count++; }
        out[i] = (count > 0) ? (sum / (double)count) : in[i];
    }
}

double oep_ts_entropy_shannon(const double* data, int n, int bins) {
    if (!data || n < 2 || bins < 2) return 0.0;
    /* Find min/max for binning */
    double vmin = data[0], vmax = data[0];
    for (int i = 1; i < n; i++) {
        if (data[i] < vmin) vmin = data[i];
        if (data[i] > vmax) vmax = data[i];
    }
    if (vmax - vmin < OEP_EPSILON) return 0.0;

    int* hist = calloc((size_t)bins, sizeof(int));
    if (!hist) return 0.0;
    double bin_width = (vmax - vmin) / (double)bins;
    for (int i = 0; i < n; i++) {
        int b = (int)((data[i] - vmin) / bin_width);
        if (b >= bins) b = bins - 1;
        if (b < 0) b = 0;
        hist[b]++;
    }
    double entropy = 0.0;
    for (int i = 0; i < bins; i++) {
        if (hist[i] > 0) {
            double p = (double)hist[i] / (double)n;
            entropy -= p * log(p);
        }
    }
    free(hist);
    return entropy;
}

double oep_ts_lyapunov_exponent(const double* data, int n, int embed_dim, int delay) {
    /* Estimate largest Lyapunov exponent via Rosenstein algorithm.
     * Reference: Rosenstein et al. (1993) Physica D 65:117-134.
     * Complexity: O(N²) naive — adequate for moderate N.
     *
     * Algorithm: reconstruct phase space, find nearest neighbors,
     * measure average divergence rate. */
    if (!data || n < embed_dim * delay + 2) return 0.0;
    int N = n - (embed_dim - 1) * delay;
    if (N < 20) return 0.0;

    /* Build embedded vectors */
    double** embedded = malloc((size_t)N * sizeof(double*));
    if (!embedded) return 0.0;
    for (int i = 0; i < N; i++) {
        embedded[i] = malloc((size_t)embed_dim * sizeof(double));
        if (!embedded[i]) {
            for (int j = 0; j < i; j++) free(embedded[j]);
            free(embedded);
            return 0.0;
        }
        for (int d = 0; d < embed_dim; d++)
            embedded[i][d] = data[i + d * delay];
    }

    /* Find nearest neighbor for each point (exclude temporally close) */
    int exclusion = 10;
    double sum_log_div = 0.0;
    int count = 0;
    for (int i = 0; i < N - 1; i++) {
        double min_dist = 1e100;
        int nn = -1;
        for (int j = 0; j < N - 1; j++) {
            if (abs(i - j) < exclusion) continue;
            double dist = 0.0;
            for (int d = 0; d < embed_dim; d++) {
                double diff = embedded[i][d] - embedded[j][d];
                dist += diff * diff;
            }
            if (dist < min_dist) { min_dist = dist; nn = j; }
        }
        if (nn >= 0 && min_dist > OEP_EPSILON) {
            /* Compute divergence after 1 step */
            double div = 0.0;
            for (int d = 0; d < embed_dim; d++) {
                double diff = embedded[i + 1][d] - embedded[nn + 1][d];
                div += diff * diff;
            }
            sum_log_div += log(sqrt(div) / sqrt(min_dist));
            count++;
        }
    }

    for (int i = 0; i < N; i++) free(embedded[i]);
    free(embedded);
    return (count > 0) ? (sum_log_div / (double)count) : 0.0;
}
