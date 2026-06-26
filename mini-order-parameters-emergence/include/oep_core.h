#ifndef OEP_CORE_H
#define OEP_CORE_H
#include <stdbool.h>
#include <stddef.h>
#include <math.h>

/* ==============================================================
 * Order Parameters & Emergence — Core Types and Primitives
 * Based on: Haken (1977, 1983), Prigogine (1980), Nicolis (1989)
 *
 * This module implements the mathematical framework for order
 * parameter identification and emergence quantification in
 * complex synergetic systems.
 *
 * Key concepts:
 *   dξ/dt = αξ - βξ³ + F(t)           — Order parameter equation
 *   τ_fast << τ_slow                    — Time-scale separation
 *   dq/dt = f(q, ξ)                     — Slaved variables
 * ============================================================== */

/* ── Constants ─────────────────────────────────────────── */
#define OEP_EPSILON      1e-14
#define OEP_MAX_DIM      128
#define OEP_MAX_ITER     10000
#define OEP_PI           3.14159265358979323846
#define OEP_E             2.71828182845904523536

/* ── Core Types ────────────────────────────────────────── */

/* Complex number for stability analysis */
typedef struct {
    double re;
    double im;
} OEPComplex;

/* Dense vector */
typedef struct {
    double* data;
    int     size;
    int     capacity;
} OEPVector;

/* Dense matrix, row-major storage */
typedef struct {
    double* data;
    int     rows;
    int     cols;
    int     capacity;
} OEPMatrix;

/* ── Vector API ────────────────────────────────────────── */

OEPVector* oep_vec_create(int size);
OEPVector* oep_vec_create_from(const double* data, int size);
OEPVector* oep_vec_clone(const OEPVector* v);
void       oep_vec_free(OEPVector* v);
void       oep_vec_set(OEPVector* v, int i, double val);
double     oep_vec_get(const OEPVector* v, int i);
void       oep_vec_set_all(OEPVector* v, const double* data);
void       oep_vec_zero(OEPVector* v);
void       oep_vec_copy(OEPVector* dst, const OEPVector* src);
void       oep_vec_print(const OEPVector* v, const char* name);

double     oep_vec_dot(const OEPVector* a, const OEPVector* b);
double     oep_vec_norm(const OEPVector* v);
double     oep_vec_norm_max(const OEPVector* v);
OEPVector* oep_vec_add(const OEPVector* a, const OEPVector* b);
OEPVector* oep_vec_sub(const OEPVector* a, const OEPVector* b);
OEPVector* oep_vec_scale(const OEPVector* v, double s);
OEPVector* oep_vec_linear_combination(double a, const OEPVector* va,
                                      double b, const OEPVector* vb);
double     oep_vec_mean(const OEPVector* v);
double     oep_vec_variance(const OEPVector* v);
double     oep_vec_correlation(const OEPVector* a, const OEPVector* b);

/* ── Matrix API ────────────────────────────────────────── */

OEPMatrix* oep_mat_create(int rows, int cols);
OEPMatrix* oep_mat_create_from(const double* data, int rows, int cols);
OEPMatrix* oep_mat_clone(const OEPMatrix* m);
void       oep_mat_free(OEPMatrix* m);
void       oep_mat_set(OEPMatrix* m, int i, int j, double val);
double     oep_mat_get(const OEPMatrix* m, int i, int j);
void       oep_mat_set_all(OEPMatrix* m, const double* data);
void       oep_mat_identity(OEPMatrix* m);
void       oep_mat_zero(OEPMatrix* m);
void       oep_mat_copy(OEPMatrix* dst, const OEPMatrix* src);
void       oep_mat_print(const OEPMatrix* m, const char* name);

OEPMatrix* oep_mat_add(const OEPMatrix* a, const OEPMatrix* b);
OEPMatrix* oep_mat_sub(const OEPMatrix* a, const OEPMatrix* b);
OEPMatrix* oep_mat_mul(const OEPMatrix* a, const OEPMatrix* b);
OEPMatrix* oep_mat_scale(const OEPMatrix* m, double s);
OEPMatrix* oep_mat_transpose(const OEPMatrix* m);
OEPVector* oep_mat_vec_mul(const OEPMatrix* A, const OEPVector* x);
double     oep_mat_determinant(const OEPMatrix* m);
double     oep_mat_trace(const OEPMatrix* m);
double     oep_mat_norm_frobenius(const OEPMatrix* m);
double     oep_mat_norm_inf(const OEPMatrix* m);
int        oep_mat_rank_gaussian(const OEPMatrix* m);

/* ODE solver utilities (embedded Runge-Kutta) */
double     oep_ode_rk4_step(double t, const double* y, int n, double dt,
                            void (*f)(double, const double*, double*, void*),
                            void* ctx, double* y_out);
void       oep_ode_rk4_multi(double t0, const double* y0, int n, int steps,
                             double dt,
                             void (*f)(double, const double*, double*, void*),
                             void* ctx, double** y_trajectory);

/* ── Utility Functions ─────────────────────────────────── */

double     oep_max(double a, double b);
double     oep_min(double a, double b);
double     oep_clamp(double x, double lo, double hi);
double     oep_sigmoid(double x);
double     oep_tanh_smooth(double x, double scale);
double     oep_gaussian(double x, double mu, double sigma);
bool       oep_is_zero(double x);
bool       oep_is_equal(double a, double b);
double     oep_sign(double x);
double     oep_heaviside(double x, double threshold);
double     oep_deg2rad(double deg);
double     oep_rad2deg(double rad);

/* Random number generation (for stochastic dynamics) */
void       oep_rng_set_seed(unsigned int seed);
double     oep_rng_uniform(void);
double     oep_rng_normal(double mu, double sigma);
double     oep_rng_exponential(double lambda);
void       oep_rng_normal_vector(double* out, int n, double mu, double sigma);

/* ── Timeseries Utilities ──────────────────────────────── */

double     oep_ts_mean(const double* data, int n);
double     oep_ts_variance(const double* data, int n);
double     oep_ts_autocorrelation(const double* data, int n, int lag);
double     oep_ts_cross_correlation(const double* x, const double* y, int n, int lag);
void       oep_ts_detrend_linear(double* data, int n);
void       oep_ts_smooth_moving_average(const double* in, double* out, int n, int window);
double     oep_ts_entropy_shannon(const double* data, int n, int bins);
double     oep_ts_lyapunov_exponent(const double* data, int n, int embed_dim, int delay);

#endif
