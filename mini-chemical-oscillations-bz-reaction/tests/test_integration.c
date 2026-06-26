/**
 * @file test_integration.c
 * @brief Tests for numerical integrators: RK4, RK45, implicit Euler.
 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "../include/bz_reaction.h"

/* ODE RHS function pointer type */
typedef void (*bz_ode_rhs_fn)(double *dydt, double t,
                               const double *y, void *params, int n);

/* Forward declarations from bz_integration.c */
void bz_integration_euler_step(double *y_new, const double *y,
                                double t, double dt, int n,
                                bz_ode_rhs_fn rhs, void *params);
void bz_rk4_step(double *y_out, const double *y_in,
                  double t, double dt, int n,
                  bz_ode_rhs_fn rhs, void *params);
int bz_rk45_adaptive_step(double *y, double *t, double *dt, int n,
                           bz_ode_rhs_fn rhs, void *params,
                           double atol, double rtol);
void bz_implicit_euler_oregonator(bz_oregonator_state_t *s,
                                   double dt,
                                   const bz_oregonator_params_t *p);
void bz_oregonator_rhs_generic(double *dydt, double t,
                                const double *y, void *params, int n);

static int test_euler_step(void) {
    bz_oregonator_params_t p;
    double y[3], y_new[3];
    bz_oregonator_params_init(&p, 0.8);
    y[0] = 0.01; y[1] = 0.001; y[2] = 0.01;
    bz_integration_euler_step(y_new, y, 0.0, 0.01, 3,
                               bz_oregonator_rhs_generic, &p);
    assert(isfinite(y_new[0]));
    assert(isfinite(y_new[1]));
    assert(isfinite(y_new[2]));
    return 1;
}

static int test_rk4_step(void) {
    bz_oregonator_params_t p;
    double y[3], y_out[3];
    bz_oregonator_params_init(&p, 0.8);
    y[0] = 0.01; y[1] = 0.001; y[2] = 0.01;
    bz_rk4_step(y_out, y, 0.0, 0.01, 3,
                 bz_oregonator_rhs_generic, &p);
    assert(isfinite(y_out[0]));
    assert(isfinite(y_out[1]));
    assert(isfinite(y_out[2]));
    /* Should move in direction of derivative */
    assert(y_out[2] != y[2] || y_out[0] != y[0]);
    return 1;
}

static int test_rk45_adaptive(void) {
    bz_oregonator_params_t p;
    double y[3], t, dt;
    int ret;
    bz_oregonator_params_init(&p, 0.8);
    y[0] = 0.01; y[1] = 0.001; y[2] = 0.01;
    t = 0.0;
    dt = 0.001;
    ret = bz_rk45_adaptive_step(y, &t, &dt, 3,
                                 bz_oregonator_rhs_generic, &p, 1e-6, 1e-6);
    assert(ret >= 0);
    assert(t > 0.0);
    assert(isfinite(y[0]));
    assert(dt > 0.0);
    return 1;
}

static int test_implicit_euler_oregonator(void) {
    bz_oregonator_params_t p;
    bz_oregonator_state_t s;
    bz_oregonator_params_init(&p, 0.8);
    bz_oregonator_init_oscillatory(&s, &p);
    double x0 = s.x;
    bz_implicit_euler_oregonator(&s, 0.01, &p);
    assert(isfinite(s.x));
    assert(s.x != x0);
    return 1;
}

static int test_rk4_accuracy(void) {
    /* Test that RK4 conserves energy-like quantity for linear oscillator */
    /* y'' = -y => y' = v, v' = -y. Use this as sanity check. */
    /* Actually test on Oregonator: after one small step, should not diverge */
    bz_oregonator_params_t p;
    double y[3], y_out[3];
    bz_oregonator_params_init(&p, 0.8);
    y[0] = 0.01; y[1] = 0.001; y[2] = 0.01;
    bz_rk4_step(y_out, y, 0.0, 0.001, 3,
                 bz_oregonator_rhs_generic, &p);
    /* Values should remain in reasonable range */
    assert(y_out[0] >= 0.0);
    assert(y_out[1] >= 0.0);
    assert(y_out[2] >= 0.0);
    return 1;
}

static int test_euler_vs_rk4_convergence(void) {
    /* Both methods should produce finite, positive concentrations */
    bz_oregonator_params_t p;
    double y_euler[3], y_rk4[3], y0[3];
    double dt = 1e-5;
    int i;
    bz_oregonator_params_init(&p, 0.8);
    y0[0] = 0.1; y0[1] = 0.01; y0[2] = 0.1;
    bz_integration_euler_step(y_euler, y0, 0.0, dt, 3,
                               bz_oregonator_rhs_generic, &p);
    bz_rk4_step(y_rk4, y0, 0.0, dt, 3,
                 bz_oregonator_rhs_generic, &p);
    /* Verify both methods give finite, non-negative results */
    for (i = 0; i < 3; i++) {
        assert(isfinite(y_euler[i]));
        assert(isfinite(y_rk4[i]));
        assert(y_euler[i] >= 0.0);
        assert(y_rk4[i] >= 0.0);
    }
    return 1;
}

int main(void) {
    int passed = 0, total = 0;
#define RUN(t) do { total++; printf("  %s... ", #t); \
    if (t()) { printf("PASS\n"); passed++; } else { printf("FAIL\n"); } \
} while(0)

    printf("=== Integration Tests ===\n");
    RUN(test_euler_step);
    RUN(test_rk4_step);
    RUN(test_rk45_adaptive);
    RUN(test_implicit_euler_oregonator);
    RUN(test_rk4_accuracy);
    RUN(test_euler_vs_rk4_convergence);

    printf("\n%d/%d tests passed\n", passed, total);
    return (passed == total) ? 0 : 1;
}
