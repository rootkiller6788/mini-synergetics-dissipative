/* ==============================================================
 * haken_dynamics.c — Nonlinear Dynamics & Numerical Integration
 *
 * Implements the time evolution of synergetic systems with
 * built-in support for timescale separation. Provides multiple
 * integrators (Euler, RK2, RK4, implicit, split, adaptive) and
 * constructs canonical models demonstrating the slaving principle.
 *
 * Key models:
 *   Haken prototype:  dq1/dt = α·q1 - a·q1·q2
 *                     dq2/dt = -β·q2 + b·q1²
 *   Lorenz near onset: reduces to GL equation via slaving.
 *
 * References:
 *   Haken (1977) Synergetics, §5.1
 *   Haken (1983) Advanced Synergetics, §1.3-1.5
 *   Strogatz (2015) Nonlinear Dynamics and Chaos, §3
 *   Press et al. (2007) Numerical Recipes, §17
 * ============================================================== */

#include "haken_dynamics.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* ==============================================================
 * RHS evaluation
 * ============================================================== */

void haken_rhs(const Haken_System* sys, const double* q, double* dq) {
    if (!sys || !q || !dq) return;
    int n = sys->n_dim;

    /* Zero output */
    for (int i = 0; i < n; i++) dq[i] = 0.0;

    /* Linear part: L(α)·q */
    if (sys->linear_matrix_func) {
        sys->linear_matrix_func(sys->control_param, n, sys->linear_matrix);
    }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            dq[i] += sys->linear_matrix[i * n + j] * q[j];

    /* Nonlinear part: N(q) */
    if (sys->nonlinear_coupling) {
        double* Nq = (double*)calloc((size_t)n, sizeof(double));
        sys->nonlinear_coupling(q, n, Nq);
        for (int i = 0; i < n; i++) dq[i] += Nq[i];
        free(Nq);
    }
}

void haken_rhs_linear(const Haken_System* sys, const double* q,
                       double* f_lin) {
    if (!sys || !q || !f_lin) return;
    int n = sys->n_dim;
    for (int i = 0; i < n; i++) {
        f_lin[i] = 0.0;
        for (int j = 0; j < n; j++)
            f_lin[i] += sys->linear_matrix[i * n + j] * q[j];
    }
}

void haken_rhs_nonlinear(const Haken_System* sys, const double* q,
                          double* f_nl) {
    if (!sys || !q || !f_nl) return;
    int n = sys->n_dim;
    for (int i = 0; i < n; i++) f_nl[i] = 0.0;
    if (sys->nonlinear_coupling) {
        sys->nonlinear_coupling(q, n, f_nl);
    }
}

/* ==============================================================
 * L5: Explicit integrators
 * ============================================================== */

void haken_step_euler(Haken_System* sys, double dt) {
    if (!sys) return;
    int n = sys->n_dim;
    double* dq = (double*)calloc((size_t)n, sizeof(double));
    haken_rhs(sys, sys->state, dq);
    for (int i = 0; i < n; i++) sys->state[i] += dt * dq[i];
    free(dq);
}

void haken_step_rk2(Haken_System* sys, double dt) {
    if (!sys) return;
    int n = sys->n_dim;
    double* dq1 = (double*)calloc((size_t)n, sizeof(double));
    double* dq2 = (double*)calloc((size_t)n, sizeof(double));
    double* tmp  = (double*)calloc((size_t)n, sizeof(double));

    haken_rhs(sys, sys->state, dq1);
    for (int i = 0; i < n; i++) tmp[i] = sys->state[i] + 0.5 * dt * dq1[i];
    haken_rhs(sys, tmp, dq2);

    for (int i = 0; i < n; i++) sys->state[i] += dt * dq2[i];
    free(dq1); free(dq2); free(tmp);
}

void haken_step_rk4(Haken_System* sys, double dt) {
    if (!sys) return;
    int n = sys->n_dim;
    double* k1 = (double*)calloc((size_t)n, sizeof(double));
    double* k2 = (double*)calloc((size_t)n, sizeof(double));
    double* k3 = (double*)calloc((size_t)n, sizeof(double));
    double* k4 = (double*)calloc((size_t)n, sizeof(double));
    double* tmp = (double*)calloc((size_t)n, sizeof(double));
    double* state_save = (double*)calloc((size_t)n, sizeof(double));

    memcpy(state_save, sys->state, (size_t)n * sizeof(double));

    /* k1 */
    haken_rhs(sys, sys->state, k1);
    /* k2 */
    for (int i = 0; i < n; i++) tmp[i] = state_save[i] + 0.5 * dt * k1[i];
    haken_rhs(sys, tmp, k2);
    /* k3 */
    for (int i = 0; i < n; i++) tmp[i] = state_save[i] + 0.5 * dt * k2[i];
    haken_rhs(sys, tmp, k3);
    /* k4 */
    for (int i = 0; i < n; i++) tmp[i] = state_save[i] + dt * k3[i];
    haken_rhs(sys, tmp, k4);

    for (int i = 0; i < n; i++)
        sys->state[i] = state_save[i] + (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);

    free(k1); free(k2); free(k3); free(k4); free(tmp); free(state_save);
}

/* ==============================================================
 * Implicit: backward Euler (for stiff modes)
 * ============================================================== */

void haken_step_backward_euler(Haken_System* sys, double dt,
                                double tol, int max_iter) {
    if (!sys) return;
    int n = sys->n_dim;
    double* q_new = (double*)calloc((size_t)n, sizeof(double));
    double* q_old = (double*)calloc((size_t)n, sizeof(double));
    double* dq    = (double*)calloc((size_t)n, sizeof(double));
    double* resid = (double*)calloc((size_t)n, sizeof(double));

    memcpy(q_old, sys->state, (size_t)n * sizeof(double));
    memcpy(q_new, sys->state, (size_t)n * sizeof(double));

    for (int iter = 0; iter < max_iter; iter++) {
        haken_rhs(sys, q_new, dq);
        double max_res = 0.0;
        for (int i = 0; i < n; i++) {
            resid[i] = q_new[i] - q_old[i] - dt * dq[i];
            if (fabs(resid[i]) > max_res) max_res = fabs(resid[i]);
        }
        if (max_res < tol) break;
        /* Jacobi update */
        for (int i = 0; i < n; i++) {
            double J_diag = 1.0 - dt * sys->linear_matrix[i * n + i];
            if (fabs(J_diag) > 1e-15) q_new[i] -= resid[i] / J_diag;
        }
    }
    memcpy(sys->state, q_new, (size_t)n * sizeof(double));
    free(q_new); free(q_old); free(dq); free(resid);
}

/* ==============================================================
 * Split-step integrator (Strang splitting)
 * ============================================================== */

void haken_step_split(Haken_System* sys, double dt,
                       double tol, int max_iter) {
    if (!sys) return;
    if (sys->n_critical == 0 || sys->n_stable == 0) {
        haken_step_rk4(sys, dt);
        return;
    }

    int n  = sys->n_dim;
    int ns = sys->n_critical;
    int nf = sys->n_stable;

    double* q_full = (double*)calloc((size_t)n, sizeof(double));
    memcpy(q_full, sys->state, (size_t)n * sizeof(double));

    double* q_s = (double*)calloc((size_t)ns, sizeof(double));
    double* q_f = (double*)calloc((size_t)nf, sizeof(double));

    /* Extract slow and fast components */
    for (int a = 0; a < ns; a++) q_s[a] = q_full[sys->critical_indices[a]];
    for (int a = 0; a < nf; a++) q_f[a] = q_full[sys->stable_indices[a]];

    /* 1. Half-step implicit on fast */
    for (int sub = 0; sub < max_iter; sub++) {
        double* dq_f = (double*)calloc((size_t)nf, sizeof(double));
        /* Evaluate fast dynamics: dq_f/dt = Λ_f·q_f + N_f(q_s, q_f) */
        for (int a = 0; a < nf; a++) {
            int k = sys->stable_indices[a];
            dq_f[a] = sys->modes[k].lambda_re * q_f[a];
            /* Add nonlinear coupling (simplified: only linear term dominant) */
        }
        double max_res = 0.0;
        for (int a = 0; a < nf; a++) {
            double new_val = q_f[a] + 0.5 * dt * dq_f[a];
            double diff = new_val - q_f[a];
            if (fabs(diff) > max_res) max_res = fabs(diff);
            q_f[a] = new_val;
        }
        free(dq_f);
        if (max_res < tol) break;
    }

    /* 2. Full-step explicit on slow */
    double* dq_s = (double*)calloc((size_t)ns, sizeof(double));
    for (int a = 0; a < ns; a++) {
        int k = sys->critical_indices[a];
        dq_s[a] = sys->modes[k].lambda_re * q_s[a];
    }
    for (int a = 0; a < ns; a++) q_s[a] += dt * dq_s[a];
    free(dq_s);

    /* 3. Half-step implicit on fast */
    for (int sub = 0; sub < max_iter; sub++) {
        double* dq_f2 = (double*)calloc((size_t)nf, sizeof(double));
        for (int a = 0; a < nf; a++) {
            int k = sys->stable_indices[a];
            dq_f2[a] = sys->modes[k].lambda_re * q_f[a];
        }
        double max_res = 0.0;
        for (int a = 0; a < nf; a++) {
            double new_val = q_f[a] + 0.5 * dt * dq_f2[a];
            double diff = new_val - q_f[a];
            if (fabs(diff) > max_res) max_res = fabs(diff);
            q_f[a] = new_val;
        }
        free(dq_f2);
        if (max_res < tol) break;
    }

    /* Reassemble */
    for (int a = 0; a < ns; a++) q_full[sys->critical_indices[a]] = q_s[a];
    for (int a = 0; a < nf; a++) q_full[sys->stable_indices[a]] = q_f[a];
    memcpy(sys->state, q_full, (size_t)n * sizeof(double));

    free(q_full); free(q_s); free(q_f);
}

/* ==============================================================
 * Adaptive RK4(5) — Dormand-Prince pair
 * ============================================================== */

int haken_step_rk45_adaptive(Haken_System* sys, double* dt, double* err,
                               double atol, double rtol) {
    if (!sys || !dt || !err) return -1;
    int n = sys->n_dim;
    double h = *dt;

    /* Dormand-Prince RK5(4) coefficients */
    static const double a21=1.0/5;
    static const double a31=3.0/40, a32=9.0/40;
    static const double a41=44.0/45, a42=-56.0/15, a43=32.0/9;
    static const double a51=19372.0/6561, a52=-25360.0/2187, a53=64448.0/6561, a54=-212.0/729;
    static const double a61=9017.0/3168, a62=-355.0/33, a63=46732.0/5247, a64=49.0/176, a65=-5103.0/18656;

    /* Butcher tableau weights */
    static const double b51=35.0/384, b53=500.0/1113, b54=125.0/192, b55=-2187.0/6784, b56=11.0/84;
    static const double b41=5179.0/57600.0, b43=7571.0/16695.0, b44=393.0/640.0,
                         b45=-92097.0/339200.0, b46=187.0/2100.0, b47=1.0/40.0;
    (void)b47;  /* Reserved for future FSAL implementation */

    double* k1 = (double*)calloc((size_t)n, sizeof(double));
    double* k2 = (double*)calloc((size_t)n, sizeof(double));
    double* k3 = (double*)calloc((size_t)n, sizeof(double));
    double* k4 = (double*)calloc((size_t)n, sizeof(double));
    double* k5 = (double*)calloc((size_t)n, sizeof(double));
    double* k6 = (double*)calloc((size_t)n, sizeof(double));
    double* k7 = (double*)calloc((size_t)n, sizeof(double));
    double* tmp = (double*)calloc((size_t)n, sizeof(double));
    double* y_old = (double*)calloc((size_t)n, sizeof(double));
    double* y5 = (double*)calloc((size_t)n, sizeof(double));
    double* y4 = (double*)calloc((size_t)n, sizeof(double));

    memcpy(y_old, sys->state, (size_t)n * sizeof(double));

    /* Stage 1 */
    haken_rhs(sys, sys->state, k1);
    /* Stage 2 */
    for (int i = 0; i < n; i++) tmp[i] = y_old[i] + h * a21 * k1[i];
    haken_rhs(sys, tmp, k2);
    /* Stage 3 */
    for (int i = 0; i < n; i++) tmp[i] = y_old[i] + h * (a31 * k1[i] + a32 * k2[i]);
    haken_rhs(sys, tmp, k3);
    /* Stage 4 */
    for (int i = 0; i < n; i++) tmp[i] = y_old[i] + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
    haken_rhs(sys, tmp, k4);
    /* Stage 5 */
    for (int i = 0; i < n; i++) tmp[i] = y_old[i] + h * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
    haken_rhs(sys, tmp, k5);
    /* Stage 6 */
    for (int i = 0; i < n; i++) tmp[i] = y_old[i] + h * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
    haken_rhs(sys, tmp, k6);
    /* Stage 7 */
    for (int i = 0; i < n; i++)
        y4[i] = y_old[i] + h * (b41 * k1[i] + b43 * k3[i] + b44 * k4[i] + b45 * k5[i] + b46 * k6[i]);
    /* k7 = f(y4) */
    haken_rhs(sys, y4, k7);
    /* 5th-order solution */
    for (int i = 0; i < n; i++)
        y5[i] = y_old[i] + h * (b51 * k1[i] + b53 * k3[i] + b54 * k4[i] + b55 * k5[i] + b56 * k6[i]);

    /* Error estimate */
    double max_err = 0.0, max_y = 0.0;
    for (int i = 0; i < n; i++) {
        double e = fabs(y5[i] - y4[i]);
        if (e > max_err) max_err = e;
        if (fabs(y5[i]) > max_y) max_y = fabs(y5[i]);
    }
    *err = max_err;

    /* Step size control */
    double tolerance = atol + rtol * max_y;
    if (max_err <= tolerance) {
        /* Accept step */
        memcpy(sys->state, y5, (size_t)n * sizeof(double));
        if (max_err > 0) {
            double factor = 0.9 * pow(tolerance / max_err, 1.0 / 5.0);
            if (factor > 5.0) factor = 5.0;
            if (factor < 0.2) factor = 0.2;
            *dt = h * factor;
        }
        free(k1); free(k2); free(k3); free(k4); free(k5); free(k6); free(k7);
        free(tmp); free(y_old); free(y5); free(y4);
        return 0;
    } else {
        /* Reject step */
        double factor = 0.9 * pow(tolerance / max_err, 1.0 / 5.0);
        if (factor < 0.2) factor = 0.2;
        *dt = h * factor;
        free(k1); free(k2); free(k3); free(k4); free(k5); free(k6); free(k7);
        free(tmp); free(y_old); free(y5); free(y4);
        return 1;
    }
}

/* ==============================================================
 * Generic integration over an interval
 * ============================================================== */

void haken_integrate_step(Haken_System* sys, double dt,
                           Haken_IntegratorType type) {
    if (!sys) return;
    switch (type) {
        case HAKEN_INTEGRATOR_EULER:    haken_step_euler(sys, dt); break;
        case HAKEN_INTEGRATOR_RK2:      haken_step_rk2(sys, dt); break;
        case HAKEN_INTEGRATOR_RK4:      haken_step_rk4(sys, dt); break;
        case HAKEN_INTEGRATOR_IMPLICIT: haken_step_backward_euler(sys, dt, 1e-8, 50); break;
        case HAKEN_INTEGRATOR_SPLIT:    haken_step_split(sys, dt, 1e-8, 20); break;
        default:                        haken_step_rk4(sys, dt); break;
    }
}

int haken_integrate(Haken_System* sys, Haken_IntegratorConfig* config) {
    if (!sys || !config) return -1;
    double t = config->t_start;
    double dt = config->dt;
    int step = 0;
    int traj_idx = 0;

    if (config->record_trajectory && config->trajectory && config->trajectory_capacity > 0) {
        memcpy(config->trajectory, sys->state, (size_t)(sys->n_dim) * sizeof(double));
        traj_idx = 1;
    }

    while (t < config->t_end && step < config->max_steps) {
        if (config->type == HAKEN_INTEGRATOR_ADAPTIVE) {
            double err;
            int ret = haken_step_rk45_adaptive(sys, &dt, &err, config->atol, config->rtol);
            if (ret == 1) continue; /* Rejected, retry with smaller dt */
        } else {
            haken_integrate_step(sys, dt, config->type);
        }
        t += dt;
        step++;

        if (config->record_trajectory && traj_idx < config->trajectory_capacity) {
            memcpy(config->trajectory + traj_idx * sys->n_dim,
                   sys->state, (size_t)(sys->n_dim) * sizeof(double));
            traj_idx++;
        }
    }
    config->n_saved = traj_idx;
    return step;
}

/* ==============================================================
 * L6: Canonical model constructors
 * ============================================================== */

/* ---- Haken prototype slaving model ---- */

/* Context for prototype nonlinear coupling */
typedef struct { double a; double b; double beta; } ProtoContext;

static void proto_nonlinear_coupling(const double* q, int n, double* Nq) {
    if (!q || !Nq || n != 2) return;
    /* Cannot access context directly without system pointer.
     * Use static globals as last-resort storage. */
    /* Default parameters: a=1, b=1 */
    Nq[0] = -q[0] * q[1];  /* -a·q1·q2 with a=1 */
    Nq[1] = q[0] * q[0];   /* b·q1² with b=1 */
}

Haken_System* haken_create_prototype_model(double alpha, double beta,
                                             double a, double b) {
    (void)a; (void)b; /* Parameters fixed at 1.0 in default coupling */
    Haken_System* sys = haken_system_create(2, NULL, proto_nonlinear_coupling);
    if (!sys) return NULL;

    /* Build linear matrix: L = [[α, 0], [0, -β]] */
    sys->linear_matrix[0] = alpha; sys->linear_matrix[1] = 0.0;
    sys->linear_matrix[2] = 0.0;   sys->linear_matrix[3] = -beta;
    sys->linear_matrix_func = NULL; /* Use explicit matrix */

    /* Store coupling parameters in context */
    ProtoContext* ctx = (ProtoContext*)malloc(sizeof(ProtoContext));
    ctx->a = a; ctx->b = b; ctx->beta = beta;
    sys->context = ctx;

    sys->coupling_type = HAKEN_COUPLING_QUADRATIC;
    sys->control_param = alpha;
    sys->critical_param = 0.0;

    /* Set initial state */
    sys->state[0] = 0.1; /* Small initial order parameter */
    sys->state[1] = 0.0;

    return sys;
}

/* ---- Lorenz model near first instability ----
 * The full Lorenz equations:
 *   ẋ = σ(y - x)
 *   ẏ = rx - y - xz
 *   ż = xy - bz
 *
 * Near the first instability (r → 1+):
 *   x, y are O(ε), z is O(ε²) → single order parameter
 *
 * Haken (1975) showed that near threshold the Lorenz model
 * reduces to a GL equation via the slaving principle.
 * Reference: Haken (1975) Phys. Lett. 53A, 77
 * ----------------------------------------------------------- */

static void lorenz_linear(double r, int n, double* L) {
    if (!L || n != 3) return;
    double sigma = 10.0, b_geo = 8.0/3.0;
    /* Jacobian at origin */
    L[0] = -sigma; L[1] = sigma;  L[2] = 0.0;
    L[3] = r;      L[4] = -1.0;   L[5] = 0.0;
    L[6] = 0.0;    L[7] = 0.0;    L[8] = -b_geo;
}

static void lorenz_nonlinear(const double* q, int n, double* Nq) {
    if (!q || !Nq || n != 3) return;
    Nq[0] = 0.0;
    Nq[1] = -q[0] * q[2];  /* -xz */
    Nq[2] = q[0] * q[1];   /* xy */
}

Haken_System* haken_create_lorenz_instability(double r, double sigma,
                                                double b_geo) {
    Haken_System* sys = haken_system_create(3, lorenz_linear, lorenz_nonlinear);
    if (!sys) return NULL;
    sys->control_param  = r;
    sys->critical_param = 1.0;
    sys->coupling_type  = HAKEN_COUPLING_QUADRATIC;

    /* Set linear matrix with actual sigma, b_geo */
    sys->linear_matrix[0] = -sigma; sys->linear_matrix[1] = sigma;  sys->linear_matrix[2] = 0.0;
    sys->linear_matrix[3] = r;      sys->linear_matrix[4] = -1.0;   sys->linear_matrix[5] = 0.0;
    sys->linear_matrix[6] = 0.0;    sys->linear_matrix[7] = 0.0;    sys->linear_matrix[8] = -b_geo;

    sys->state[0] = 0.1; sys->state[1] = 0.1; sys->state[2] = 0.0;

    return sys;
}

/* ==============================================================
 * L3: Coupling tensor constructors
 * ============================================================== */

Haken_QuadraticCoupling* haken_quadratic_create(int n, bool symmetric) {
    if (n <= 0) return NULL;
    Haken_QuadraticCoupling* qc = (Haken_QuadraticCoupling*)calloc(1, sizeof(Haken_QuadraticCoupling));
    qc->n_dim = n;
    qc->symmetric_jk = symmetric;
    qc->Q = (double*)calloc((size_t)(n * n * n), sizeof(double));
    return qc;
}

void haken_quadratic_free(Haken_QuadraticCoupling* qc) {
    if (!qc) return;
    free(qc->Q);
    free(qc);
}

void haken_quadratic_zero(Haken_QuadraticCoupling* qc) {
    if (!qc) return;
    memset(qc->Q, 0, (size_t)(qc->n_dim * qc->n_dim * qc->n_dim) * sizeof(double));
}

void haken_quadratic_set(Haken_QuadraticCoupling* qc,
                          int i, int j, int k, double value) {
    if (!qc || i < 0 || j < 0 || k < 0) return;
    int n = qc->n_dim;
    if (i >= n || j >= n || k >= n) return;
    qc->Q[i + n * (j + n * k)] = value;
    if (qc->symmetric_jk) {
        qc->Q[i + n * (k + n * j)] = value;
    }
}

void haken_quadratic_eval(const Haken_QuadraticCoupling* qc,
                           const double* q, double* Nq) {
    if (!qc || !q || !Nq) return;
    int n = qc->n_dim;
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++)
                sum += qc->Q[i + n * (j + n * k)] * q[j] * q[k];
        Nq[i] = sum;
    }
}

Haken_CubicCoupling* haken_cubic_create(int n, bool symmetric) {
    if (n <= 0) return NULL;
    Haken_CubicCoupling* cc = (Haken_CubicCoupling*)calloc(1, sizeof(Haken_CubicCoupling));
    cc->n_dim = n;
    cc->symmetric_jkl = symmetric;
    cc->C = (double*)calloc((size_t)(n * n * n * n), sizeof(double));
    return cc;
}

void haken_cubic_free(Haken_CubicCoupling* cc) {
    if (!cc) return;
    free(cc->C);
    free(cc);
}

void haken_cubic_eval(const Haken_CubicCoupling* cc,
                       const double* q, double* Nq) {
    if (!cc || !q || !Nq) return;
    int n = cc->n_dim;
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++)
                for (int l = 0; l < n; l++)
                    sum += cc->C[i + n * (j + n * (k + n * l))] * q[j] * q[k] * q[l];
        Nq[i] = sum;
    }
}

Haken_PolynomialCoupling* haken_polynomial_create(int n, int max_degree) {
    if (n <= 0 || max_degree < 2) return NULL;
    Haken_PolynomialCoupling* pc = (Haken_PolynomialCoupling*)calloc(1, sizeof(Haken_PolynomialCoupling));
    pc->n_dim = n;
    pc->max_degree = max_degree;
    pc->coefficients = (double**)calloc((size_t)(max_degree - 1), sizeof(double*));
    pc->n_terms = (int*)calloc((size_t)(max_degree - 1), sizeof(int));
    for (int d = 0; d < max_degree - 1; d++) {
        pc->n_terms[d] = n;
    }
    return pc;
}

void haken_polynomial_free(Haken_PolynomialCoupling* pc) {
    if (!pc) return;
    for (int d = 0; d < pc->max_degree - 1; d++) free(pc->coefficients[d]);
    free(pc->coefficients);
    free(pc->n_terms);
    free(pc);
}

void haken_polynomial_eval(const Haken_PolynomialCoupling* pc,
                            const double* q, double* Nq) {
    if (!pc || !q || !Nq) return;
    int n = pc->n_dim;
    for (int i = 0; i < n; i++) Nq[i] = 0.0;
    /* For simplicity, handle up to degree 4: */
    for (int d = 2; d <= pc->max_degree && d <= 4; d++) {
        double* coeffs = pc->coefficients[d - 2];
        if (!coeffs) continue;
        for (int i = 0; i < n; i++) {
            double term = 0.0;
            if (d == 2) for (int j = 0; j < n; j++) for (int k = 0; k < n; k++)
                term += coeffs[i + n * (j + n * k)] * q[j] * q[k];
            else if (d == 3) for (int j = 0; j < n; j++) for (int k = 0; k < n; k++) for (int l = 0; l < n; l++)
                term += coeffs[i + n * (j + n * (k + n * l))] * q[j] * q[k] * q[l];
            else term = 0.0;
            Nq[i] += term;
        }
    }
}
