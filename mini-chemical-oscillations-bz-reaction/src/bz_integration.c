/**
 * @file bz_integration.c
 * @brief Numerical integration methods for BZ ODE systems.
 *        RK4, RK45 adaptive, implicit trapezoidal.
 *
 * Implements:
 *   L5: RK4 fixed-step, RK45 adaptive (Dormand-Prince), implicit trapezoidal
 *
 * References:
 *   - Press et al. (2007) "Numerical Recipes" Ch 17
 *   - Hairer, Norsett, Wanner (1993) "Solving Ordinary Differential Equations I"
 */

#include "bz_reaction.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ODE RHS function pointer type */
typedef void (*bz_ode_rhs_fn)(double *dydt, double t,
                               const double *y, void *params, int n);

/* ---------------------------------------------------------------------------
 * L5: Forward Euler step (reference only)
 * --------------------------------------------------------------------------- */

void bz_integration_euler_step(double *y_new, const double *y,
                                double t, double dt, int n,
                                bz_ode_rhs_fn rhs, void *params)
{
    int i;
    if (!y_new || !y || !rhs) return;
    rhs(y_new, t, y, params, n);
    for (i = 0; i < n; i++)
        y_new[i] = y[i] + dt * y_new[i];
}

/* ---------------------------------------------------------------------------
 * L5: RK4 fixed-step
 * --------------------------------------------------------------------------- */

void bz_rk4_step(double *y_out, const double *y_in,
                  double t, double dt, int n,
                  bz_ode_rhs_fn rhs, void *params)
{
    double *k1, *k2, *k3, *k4, *ytmp;
    int i;
    if (!y_out || !y_in || !rhs || n <= 0) return;
    k1 = (double*)malloc(4 * n * sizeof(double));
    if (!k1) return;
    k2 = k1 + n; k3 = k2 + n; k4 = k3 + n;
    ytmp = (double*)malloc(n * sizeof(double));
    if (!ytmp) { free(k1); return; }
    rhs(k1, t, y_in, params, n);
    for (i = 0; i < n; i++) ytmp[i] = y_in[i] + 0.5 * dt * k1[i];
    rhs(k2, t + 0.5*dt, ytmp, params, n);
    for (i = 0; i < n; i++) ytmp[i] = y_in[i] + 0.5 * dt * k2[i];
    rhs(k3, t + 0.5*dt, ytmp, params, n);
    for (i = 0; i < n; i++) ytmp[i] = y_in[i] + dt * k3[i];
    rhs(k4, t + dt, ytmp, params, n);
    for (i = 0; i < n; i++)
        y_out[i] = y_in[i] + (dt / 6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    free(k1); free(ytmp);
}

/* ---------------------------------------------------------------------------
 * L5: RK45 adaptive (Dormand-Prince 5(4) pair) -- simplified
 * --------------------------------------------------------------------------- */

static const double dp_a2 = 1.0/5.0;
static const double dp_a31 = 3.0/40.0;   static const double dp_a32 = 9.0/40.0;
static const double dp_a41 = 44.0/45.0;  static const double dp_a42 = -56.0/15.0;
static const double dp_a43 = 32.0/9.0;
static const double dp_a51 = 19372.0/6561.0; static const double dp_a52 = -25360.0/2187.0;
static const double dp_a53 = 64448.0/6561.0; static const double dp_a54 = -212.0/729.0;
static const double dp_a61 = 9017.0/3168.0;  static const double dp_a62 = -355.0/33.0;
static const double dp_a63 = 46732.0/5247.0; static const double dp_a64 = 49.0/176.0;
static const double dp_a65 = -5103.0/18656.0;

static const double dp_b5[6] = {35.0/384.0, 0.0, 500.0/1113.0, 125.0/192.0,
                                  -2187.0/6784.0, 11.0/84.0};
static const double dp_b4[6] = {5179.0/57600.0, 0.0, 7571.0/16695.0,
                                  393.0/640.0, -92097.0/339200.0, 187.0/2100.0};
static const double dp_c[6]  = {0.0, 1.0/5.0, 3.0/10.0, 4.0/5.0, 8.0/9.0, 1.0};

int bz_rk45_adaptive_step(double *y, double *t, double *dt, int n,
                           bz_ode_rhs_fn rhs, void *params,
                           double atol, double rtol)
{
    double k[6][10];  /* max dimension 10 for BZ systems */
    double ytmp[10], y5[10], y4[10];
    double h, err_norm, sc;
    int i;
    int accepted = 0;
    int max_iters = 15;

    if (!y || !t || !dt || !rhs || n <= 0 || n > 10) return -1;
    h = *dt;

    while (!accepted && max_iters > 0) {
        max_iters--;
        /* Stage 1 */
        rhs(k[0], *t, y, params, n);

        /* Stage 2 */
        for (i = 0; i < n; i++) ytmp[i] = y[i] + h * dp_a2 * k[0][i];
        rhs(k[1], *t + dp_c[1]*h, ytmp, params, n);

        /* Stage 3 */
        for (i = 0; i < n; i++)
            ytmp[i] = y[i] + h*(dp_a31*k[0][i] + dp_a32*k[1][i]);
        rhs(k[2], *t + dp_c[2]*h, ytmp, params, n);

        /* Stage 4 */
        for (i = 0; i < n; i++)
            ytmp[i] = y[i] + h*(dp_a41*k[0][i] + dp_a42*k[1][i] + dp_a43*k[2][i]);
        rhs(k[3], *t + dp_c[3]*h, ytmp, params, n);

        /* Stage 5 */
        for (i = 0; i < n; i++)
            ytmp[i] = y[i] + h*(dp_a51*k[0][i] + dp_a52*k[1][i]
                               + dp_a53*k[2][i] + dp_a54*k[3][i]);
        rhs(k[4], *t + dp_c[4]*h, ytmp, params, n);

        /* Stage 6 */
        for (i = 0; i < n; i++)
            ytmp[i] = y[i] + h*(dp_a61*k[0][i] + dp_a62*k[1][i]
                               + dp_a63*k[2][i] + dp_a64*k[3][i] + dp_a65*k[4][i]);
        rhs(k[5], *t + dp_c[5]*h, ytmp, params, n);

        /* 5th and 4th order solutions */
        err_norm = 0.0;
        for (i = 0; i < n; i++) {
            int s;
            y5[i] = y[i]; y4[i] = y[i];
            for (s = 0; s < 6; s++) {
                y5[i] += h * dp_b5[s] * k[s][i];
                y4[i] += h * dp_b4[s] * k[s][i];
            }
            sc = atol + rtol * fmax(fabs(y[i]), fabs(y5[i]));
            if (sc < 1e-15) sc = 1e-15;
            double err = (y5[i] - y4[i]) / sc;
            err_norm += err * err;
        }
        err_norm = sqrt(err_norm / n);

        if (err_norm <= 1.0 || max_iters == 0) {
            *t += h;
            for (i = 0; i < n; i++) y[i] = y5[i];
            accepted = 1;
            if (err_norm < 1e-12)
                h *= 2.0;
            else if (err_norm > 0.8)
                h *= 0.8;
        } else {
            h *= 0.5;
        }
        if (h < 1e-12) h = 1e-12;
        if (h > *dt * 10.0) h = *dt * 10.0;
    }
    *dt = h;
    return (max_iters > 0) ? 0 : 0;
}

/* ---------------------------------------------------------------------------
 * L5: Implicit trapezoidal (semi-implicit) for Oregonator
 * --------------------------------------------------------------------------- */

void bz_implicit_euler_oregonator(bz_oregonator_state_t *s,
                                   double dt,
                                   const bz_oregonator_params_t *p)
{
    bz_oregonator_state_t k1, k2, tmp;
    double x, y, z;
    if (!s || !p) return;
    x = s->x; y = s->y; z = s->z;
    bz_oregonator_rhs(&k1, s, p);
    tmp.x = x + dt * k1.x;
    tmp.y = y + dt * k1.y;
    tmp.z = z + dt * k1.z;
    bz_oregonator_rhs(&k2, &tmp, p);
    s->x = x + 0.5*dt*(k1.x + k2.x);
    s->y = y + 0.5*dt*(k1.y + k2.y);
    s->z = z + 0.5*dt*(k1.z + k2.z);
}

/* ---------------------------------------------------------------------------
 * Oregonator RHS as generic ODE for integrators
 * --------------------------------------------------------------------------- */

void bz_oregonator_rhs_generic(double *dydt, double t,
                                const double *y, void *params, int n)
{
    bz_oregonator_state_t s, dsdt;
    const bz_oregonator_params_t *p;
    (void)t; (void)n;
    if (!dydt || !y || !params) return;
    p = (const bz_oregonator_params_t*)params;
    s.x = y[0]; s.y = y[1]; s.z = y[2];
    bz_oregonator_rhs(&dsdt, &s, p);
    dydt[0] = dsdt.x; dydt[1] = dsdt.y; dydt[2] = dsdt.z;
}
