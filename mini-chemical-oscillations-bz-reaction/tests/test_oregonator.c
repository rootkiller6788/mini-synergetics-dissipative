/**
 * @file test_oregonator.c
 * @brief Tests for Oregonator core: initialization, RHS, steady states,
 *        nullclines, oscillation detection.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "../include/bz_reaction.h"
#include "../include/bz_oregonator.h"

static int test_state_init(void) {
    bz_state_t s;
    bz_state_init_typical(&s);
    assert(fabs(s.bro3 - 0.06) < 1e-6);
    assert(fabs(s.ma   - 0.25) < 1e-6);
    assert(fabs(s.ce3  - 0.001) < 1e-6);
    assert(fabs(s.h    - 0.8) < 1e-6);
    assert(fabs(s.hbro2) < 1e-15);
    return 1;
}

static int test_params_init(void) {
    bz_params_t p;
    bz_params_init_default(&p);
    assert(fabs(p.k_a1 - 2.1) < 1e-3);
    assert(fabs(p.k_a2 - 2.0e9) < 1e5);
    assert(fabs(p.k_b1 - 1.4) < 1e-3);
    assert(fabs(p.k_b3 - 3.0e3) < 1.0);
    assert(fabs(p.f     - 0.8) < 1e-6);
    assert(fabs(p.epsilon   - 9.9e-3) < 1e-5);
    return 1;
}

static int test_oregonator_params_init(void) {
    bz_oregonator_params_t p;
    bz_oregonator_params_init(&p, 1.2);
    assert(fabs(p.f - 1.2) < 1e-6);
    assert(fabs(p.q - 8.0e-5) < 1e-10);
    assert(fabs(p.epsilon - 9.9e-3) < 1e-5);
    assert(fabs(p.epsilon_p - 2.0e-5) < 1e-8);
    return 1;
}

static int test_oregonator_init_oscillatory(void) {
    bz_oregonator_params_t p;
    bz_oregonator_state_t s;
    bz_oregonator_params_init(&p, 0.8);
    bz_oregonator_init_oscillatory(&s, &p);
    assert(s.x > 0.0);
    assert(s.y > 0.0);
    assert(s.z > 0.0);
    return 1;
}

static int test_oregonator_rhs(void) {
    bz_oregonator_params_t p;
    bz_oregonator_state_t s, dsdt;
    bz_oregonator_params_init(&p, 0.8);
    bz_oregonator_init_oscillatory(&s, &p);
    bz_oregonator_rhs(&dsdt, &s, &p);
    /* Derivatives should be finite */
    assert(isfinite(dsdt.x));
    assert(isfinite(dsdt.y));
    assert(isfinite(dsdt.z));
    return 1;
}

static int test_fkn_rhs(void) {
    bz_state_t s, dsdt;
    bz_params_t p;
    double norm;
    bz_state_init_typical(&s);
    bz_params_init_default(&p);
    norm = bz_fkn_rhs(&dsdt, &s, &p);
    assert(norm >= 0.0);
    assert(isfinite(norm));
    /* Check bromide is consumed (negative derivative in process A) */
    assert(dsdt.br < 0.0);
    return 1;
}

static int test_tyson_2var(void) {
    bz_oregonator_params_t p;
    double dudt, dvdt;
    bz_oregonator_params_init(&p, 0.8);
    bz_tyson_2var_rhs(&dudt, &dvdt, 0.1, 0.1, &p);
    assert(isfinite(dudt));
    assert(isfinite(dvdt));
    return 1;
}

static int test_steady_states(void) {
    bz_oregonator_params_t p;
    bz_oregonator_state_t ss[4];
    int n;
    bz_oregonator_params_init(&p, 0.8);
    n = bz_oregonator_steady_states(ss, &p);
    assert(n >= 1);
    assert(n <= 3);
    /* Trivial steady state should be first */
    assert(fabs(ss[0].x) < 1e-15);
    return 1;
}

static int test_nullclines(void) {
    bz_oregonator_params_t p;
    double ynx, yny, zn;
    bz_oregonator_params_init(&p, 0.8);
    bz_oregonator_nullclines(&ynx, &yny, &zn, 0.2, 0.15, &p);
    assert(isfinite(ynx));
    assert(isfinite(yny));
    assert(fabs(zn - 0.2) < 1e-10);
    return 1;
}

static int test_is_oscillatory(void) {
    bz_oregonator_params_t p;
    bz_oregonator_params_init(&p, 0.8);
    /* f=0.8 should be in the oscillatory regime for q=8e-5 */
    assert(bz_oregonator_is_oscillatory(&p) == 1);
    /* f=0.1 should be steady state */
    bz_oregonator_params_init(&p, 0.1);
    assert(bz_oregonator_is_oscillatory(&p) == 0);
    return 1;
}

static int test_period_approx(void) {
    bz_oregonator_params_t p;
    double T;
    bz_oregonator_params_init(&p, 0.8);
    T = bz_oregonator_period_approx(&p);
    assert(T > 0.0);
    assert(isfinite(T));
    return 1;
}

static int test_hopf_critical_f(void) {
    double fc;
    fc = bz_oregonator_hopf_critical_f(8.0e-5);
    assert(fc > 0.4);
    assert(fc < 0.7);
    fc = bz_oregonator_hopf_critical_f(0.5);
    assert(fc > 1.0);
    return 1;
}

static int test_jacobian(void) {
    bz_oregonator_params_t p;
    bz_oregonator_state_t s;
    double J[9];
    bz_oregonator_params_init(&p, 0.8);
    bz_oregonator_init_oscillatory(&s, &p);
    bz_oregonator_jacobian(J, &s, &p);
    /* All entries should be finite */
    int i;
    for (i = 0; i < 9; i++) assert(isfinite(J[i]));
    /* J[8] = -1 (from dz/dt = x - z) */
    assert(fabs(J[8] + 1.0) < 1e-10);
    return 1;
}

static int test_eigenvalues(void) {
    bz_oregonator_params_t p;
    bz_oregonator_state_t s;
    double J[9], ev_re[3], ev_im[3];
    int n_unstable;
    bz_oregonator_params_init(&p, 0.8);
    bz_oregonator_init_oscillatory(&s, &p);
    bz_oregonator_jacobian(J, &s, &p);
    n_unstable = bz_eigenvalues_3x3(ev_re, ev_im, J);
    assert(n_unstable >= 0 && n_unstable <= 3);
    return 1;
}

static int test_hopf_test(void) {
    bz_oregonator_params_t p;
    bz_oregonator_state_t ss[4];
    int n;
    bz_oregonator_params_init(&p, 0.6);
    n = bz_oregonator_steady_states(ss, &p);
    assert(n >= 1);
    /* Hopf test runs without crashing */
    bz_oregonator_hopf_test(&p, &ss[1]); /* non-trivial SS */
    return 1;
}

static int test_sweep_f(void) {
    bz_oregonator_params_t p;
    int results[10];
    bz_oregonator_params_init(&p, 1.0);
    bz_oregonator_sweep_f(0.1, 2.0, 10, results, &p);
    /* Should find both non-oscillatory and oscillatory regions */
    assert(results[0] == 0); /* f=0.1: non-osc */
    assert(results[9] == 1); /* f=2.0: osc */
    return 1;
}

static int test_integrate_and_analyze(void) {
    bz_oregonator_params_t p;
    bz_oregonator_state_t s;
    int n_peaks;
    double period, amplitude;
    bz_oregonator_params_init(&p, 0.8);
    bz_oregonator_init_oscillatory(&s, &p);
    bz_oregonator_integrate_and_analyze(&s, &p, 50.0, 5000,
                                         &n_peaks, &period, &amplitude);
    assert(amplitude > 0.0);
    assert(period > 0.0 || n_peaks == 0);
    return 1;
}

int main(void) {
    int passed = 0, total = 0;
#define RUN(t) do { total++; printf("  %s... ", #t); \
    if (t()) { printf("PASS\n"); passed++; } else { printf("FAIL\n"); } \
} while(0)

    printf("=== Oregonator Tests ===\n");
    RUN(test_state_init);
    RUN(test_params_init);
    RUN(test_oregonator_params_init);
    RUN(test_oregonator_init_oscillatory);
    RUN(test_oregonator_rhs);
    RUN(test_fkn_rhs);
    RUN(test_tyson_2var);
    RUN(test_steady_states);
    RUN(test_nullclines);
    RUN(test_is_oscillatory);
    RUN(test_period_approx);
    RUN(test_hopf_critical_f);
    RUN(test_jacobian);
    RUN(test_eigenvalues);
    RUN(test_hopf_test);
    RUN(test_sweep_f);
    RUN(test_integrate_and_analyze);

    printf("\n%d/%d tests passed\n", passed, total);
    return (passed == total) ? 0 : 1;
}
