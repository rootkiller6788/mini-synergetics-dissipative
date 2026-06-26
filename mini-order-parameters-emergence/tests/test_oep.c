#include "oep_core.h"
#include "oep_order_parameter.h"
#include "oep_slaving_principle.h"
#include "oep_phase_transition.h"
#include "oep_emergence_metrics.h"
#include "oep_collective_dynamics.h"
#include "oep_landau_ginzburg.h"
#include "oep_applications.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

/* ==============================================================
 * Comprehensive Test Suite for mini-order-parameters-emergence
 *
 * Covers all L1-L6 core APIs with mathematical assertions.
 * Each test verifies a specific knowledge point.
 * ============================================================== */

static int tests_passed = 0;
static int tests_total = 0;

#define TEST(name) do { \
    tests_total++; \
    printf("  TEST %s ... ", name); \
} while(0)

#define PASS() do { \
    printf("PASS\n"); \
    tests_passed++; \
} while(0)

/* ── L1: Core Definitions Tests ────────────────────────── */

void test_vector_create(void) {
    TEST("vector_create");
    OEPVector* v = oep_vec_create(5);
    assert(v != NULL);
    assert(v->size == 5);
    oep_vec_free(v);
    PASS();
}

void test_matrix_create(void) {
    TEST("matrix_create");
    OEPMatrix* m = oep_mat_create(3, 4);
    assert(m != NULL);
    assert(m->rows == 3);
    assert(m->cols == 4);
    oep_mat_free(m);
    PASS();
}

void test_vector_operations(void) {
    TEST("vector_operations");
    double d1[] = {1, 2, 3}, d2[] = {4, 5, 6};
    OEPVector* v1 = oep_vec_create_from(d1, 3);
    OEPVector* v2 = oep_vec_create_from(d2, 3);
    assert(v1 && v2);

    double dot = oep_vec_dot(v1, v2);
    assert(fabs(dot - 32.0) < 1e-10); /* 1*4 + 2*5 + 3*6 = 32 */

    double norm = oep_vec_norm(v1);
    assert(fabs(norm - sqrt(14.0)) < 1e-10);

    OEPVector* v3 = oep_vec_add(v1, v2);
    assert(fabs(oep_vec_get(v3, 0) - 5.0) < 1e-10);
    assert(fabs(oep_vec_get(v3, 1) - 7.0) < 1e-10);
    assert(fabs(oep_vec_get(v3, 2) - 9.0) < 1e-10);

    oep_vec_free(v1); oep_vec_free(v2); oep_vec_free(v3);
    PASS();
}

void test_matrix_multiplication(void) {
    TEST("matrix_multiplication");
    double dA[] = {1, 2, 3, 4};  /* 2x2 */
    double dB[] = {5, 6, 7, 8};  /* 2x2 */
    OEPMatrix* A = oep_mat_create_from(dA, 2, 2);
    OEPMatrix* B = oep_mat_create_from(dB, 2, 2);
    OEPMatrix* C = oep_mat_mul(A, B);
    assert(C != NULL);
    /* [1 2; 3 4] * [5 6; 7 8] = [19 22; 43 50] */
    assert(fabs(oep_mat_get(C, 0, 0) - 19.0) < 1e-10);
    assert(fabs(oep_mat_get(C, 0, 1) - 22.0) < 1e-10);
    assert(fabs(oep_mat_get(C, 1, 0) - 43.0) < 1e-10);
    assert(fabs(oep_mat_get(C, 1, 1) - 50.0) < 1e-10);
    oep_mat_free(A); oep_mat_free(B); oep_mat_free(C);
    PASS();
}

void test_determinant(void) {
    TEST("determinant");
    double dA[] = {4, 3, 6, 3};  /* det = 4*3 - 3*6 = -6 */
    OEPMatrix* A = oep_mat_create_from(dA, 2, 2);
    double det = oep_mat_determinant(A);
    assert(fabs(det - (-6.0)) < 1e-6);
    oep_mat_free(A);
    PASS();
}

void test_rank(void) {
    TEST("matrix_rank");
    double dA[] = {1, 2, 2, 4};  /* rank 1 (second row = 2*first row) */
    OEPMatrix* A = oep_mat_create_from(dA, 2, 2);
    int r = oep_mat_rank_gaussian(A);
    assert(r == 1);
    oep_mat_free(A);
    PASS();
}

/* ── L2: Order Parameter Tests ──────────────────────────── */

void test_scalar_op_create(void) {
    TEST("scalar_op_create");
    OEPScalarOP* op = oep_op_scalar_create(1.0, 2.0, 0.1);
    assert(op != NULL);
    assert(fabs(op->alpha - 1.0) < 1e-10);
    assert(fabs(op->beta - 2.0) < 1e-10);
    oep_op_scalar_free(op);
    PASS();
}

void test_order_parameter_equilibria_supercritical(void) {
    TEST("op_equilibria_supercritical");
    double eq0, eq1, eq2;
    int n;
    oep_op_scalar_equilibria(1.0, 1.0, &eq0, &eq1, &eq2, &n);
    assert(n == 3);  /* α > 0 → 3 equilibria */
    assert(fabs(eq0) < 1e-10);   /* unstable origin */
    assert(fabs(eq1 - 1.0) < 1e-10);
    assert(fabs(eq2 + 1.0) < 1e-10);
    PASS();
}

void test_order_parameter_equilibria_subcritical(void) {
    TEST("op_equilibria_subcritical");
    double eq0, eq1, eq2;
    int n;
    oep_op_scalar_equilibria(-1.0, 1.0, &eq0, &eq1, &eq2, &n);
    assert(n == 1);  /* α < 0 → only ξ = 0 (stable) */
    assert(fabs(eq0) < 1e-10);
    PASS();
}

void test_stability_check(void) {
    TEST("stability_check");
    /* At ξ=0, α=-1 → should be stable (λ = -1) */
    assert(oep_op_scalar_is_stable(0.0, -1.0, 1.0));
    /* At ξ=0, α=1 → should be unstable (λ = 1) */
    assert(!oep_op_scalar_is_stable(0.0, 1.0, 1.0));
    /* At ξ=1 (nonzero equilibrium for α=1, β=1), λ = 1-3 = -2 → stable */
    assert(oep_op_scalar_is_stable(1.0, 1.0, 1.0));
    PASS();
}

void test_potential_symmetry(void) {
    TEST("potential_symmetry");
    double V_pos = oep_op_scalar_potential(1.5, 1.0, 1.0);
    double V_neg = oep_op_scalar_potential(-1.5, 1.0, 1.0);
    assert(fabs(V_pos - V_neg) < 1e-10);  /* V(ξ) = V(-ξ) */
    PASS();
}

void test_op_evolution(void) {
    TEST("op_evolution_deterministic");
    OEPScalarOP* op = oep_op_scalar_create(1.0, 1.0, 0.0);
    op->value = 0.1;
    int n_steps = 500;
    double* traj = malloc((size_t)n_steps * sizeof(double));
    assert(traj != NULL);
    oep_op_scalar_evolve(op, 0.02, n_steps, traj);
    /* Should approach ξ_eq = √(1/1) = 1.0 after sufficient time */
    assert(traj[n_steps - 1] > 0.5);  /* Growing toward equilibrium */
    assert(fabs(oep_op_scalar_drift(traj[n_steps - 1], 1.0, 1.0)) < 0.1);  /* Near fixed point */
    free(traj);
    oep_op_scalar_free(op);
    PASS();
}

void test_complex_op_create(void) {
    TEST("complex_op_create");
    OEPComplexOP* op = oep_op_complex_create(0.5, 2.0, 1.0, 0.0, 1.0);
    assert(op != NULL);
    assert(op->amplitude > 0.0);
    oep_op_complex_free(op);
    PASS();
}

/* ── L3: Slaving Principle Tests ───────────────────────── */

void test_slaving_create(void) {
    TEST("slaving_create");
    OEPSlavingSystem* sys = oep_slaving_create(2, 5);
    assert(sys != NULL);
    assert(sys->n_order_params == 2);
    assert(sys->slaved->n_vars == 5);
    oep_slaving_free(sys);
    PASS();
}

void test_slaving_timescale_check(void) {
    TEST("slaving_timescale_check");
    OEPSlavingSystem* sys = oep_slaving_create(1, 2);
    assert(sys != NULL);
    /* Set slow mode: α = 0.1, fast modes: γ = 10, 20 */
    oep_slaving_set_linear(sys, 0, 0.1);
    oep_slaving_set_damping(sys, 0, 10.0);
    oep_slaving_set_damping(sys, 1, 20.0);
    double ratio;
    bool ok = oep_slaving_check_timescale_separation(sys, &ratio);
    assert(ok);  /* γ_min / |α_max| = 10 / 0.1 = 100 >> 10 */
    assert(ratio > 10.0);
    oep_slaving_free(sys);
    PASS();
}

void test_quasistationary_q(void) {
    TEST("quasistationary_q");
    OEPSlavingSystem* sys = oep_slaving_create(1, 2);
    assert(sys != NULL);
    oep_slaving_set_linear(sys, 0, 0.5);
    oep_slaving_set_damping(sys, 0, 5.0);
    oep_slaving_set_damping(sys, 1, 5.0);
    oep_slaving_set_coupling(sys, 0, 0, 1.0);  /* order=0, slaved=0: K[0][0]=1 */
    oep_slaving_set_coupling(sys, 0, 1, 2.0);  /* order=0, slaved=1: K[1][0]=2 */

    double xi_data[] = {1.0};
    OEPVector* xi = oep_vec_create_from(xi_data, 1);
    OEPVector* q = oep_slaving_quasistationary_q(sys, xi);
    assert(q != NULL);
    /* q[0] = K[0,0]*(1+0.1*1)/γ₀ = 1*1.1/5 = 0.22 */
    /* q[1] = K[1,0]*(1+0.1*1)/γ₁ = 2*1.1/5 = 0.44 */
    assert(fabs(oep_vec_get(q, 0) - 0.22) < 0.01);
    assert(fabs(oep_vec_get(q, 1) - 0.44) < 0.01);

    oep_vec_free(xi); oep_vec_free(q);
    oep_slaving_free(sys);
    PASS();
}

void test_slaving_adiabatic_elimination(void) {
    TEST("slaving_adiabatic_elimination");
    OEPSlavingSystem* sys = oep_slaving_create(2, 1);
    assert(sys != NULL);
    oep_slaving_set_linear(sys, 0, 0.2);
    oep_slaving_set_linear(sys, 1, 0.3);
    oep_slaving_set_damping(sys, 0, 5.0);
    oep_slaving_set_coupling(sys, 0, 0, 1.0);
    oep_slaving_set_coupling(sys, 0, 1, 0.5);
    oep_slaving_set_reverse_coupling(sys, 0, 0, 0.1);

    OEPVector* eff_alpha = oep_vec_create(2);
    OEPVector* eff_beta = oep_vec_create(2);
    OEPMatrix* eff_coup = oep_mat_create(2, 2);
    oep_slaving_adiabatic_elimination(sys, eff_alpha, eff_beta, eff_coup);
    /* Effective α should include corrections from slaved vars */
    assert(oep_vec_get(eff_alpha, 0) > 0.2);  /* should have positive correction */
    oep_vec_free(eff_alpha); oep_vec_free(eff_beta); oep_mat_free(eff_coup);
    oep_slaving_free(sys);
    PASS();
}

/* ── L4: Phase Transition Tests ────────────────────────── */

void test_pitchfork_bifurcation(void) {
    TEST("pitchfork_bifurcation");
    OEPBifDiag* bd = oep_bifdiag_pitchfork(-2.0, 2.0, 100, 0.0, 1.0, 1.0);
    assert(bd != NULL);
    assert(bd->n_branches == 3);
    /* At p=1 (p > pc, α > 0): nonzero branches should be nonzero */
    /* Check a few points */
    assert(fabs(bd->branches[0][75].equilibrium - 0.0) < 1e-6); /* p=1, branch 0 = 0 */
    /* Branch 1 at p=1 should be sqrt(1) = 1 */
    double xi_expected = sqrt(1.0);  /* α = 1*(1-0) = 1, β = 1, sqrt(1/1) = 1 */
    assert(fabs(bd->branches[1][75].equilibrium - xi_expected) < 0.1);
    oep_bifdiag_free(bd);
    PASS();
}

void test_phase_transition_detection(void) {
    TEST("phase_transition_detection");
    int n = 50;
    double* cp = malloc((size_t)n * sizeof(double));
    double* op = malloc((size_t)n * sizeof(double));
    for (int i = 0; i < n; i++) {
        cp[i] = -1.0 + (double)i * 2.0 / (double)(n - 1);
        double alpha = 1.0 * cp[i];
        op[i] = (alpha > 0.0) ? sqrt(alpha) : 0.0;  /* Simulated pitchfork */
    }
    OEPPhaseTransition* pt = oep_pt_detect(cp, op, n);
    assert(pt != NULL);
    assert(fabs(pt->critical_param) < 0.2);  /* Should be near 0 */
    oep_pt_free(pt);
    free(cp); free(op);
    PASS();
}

void test_landau_free_energy(void) {
    TEST("landau_free_energy");
    /* At p = pc, a = 0, F should be purely quartic */
    double F = oep_pt_landau_free_energy(1.0, 0.0, 0.0, 1.0, 1.0);
    assert(fabs(F - 0.25) < 1e-10);  /* 0 + 1*1⁴/4 = 0.25 */
    PASS();
}

void test_rushbrooke_inequality(void) {
    TEST("rushbrooke_inequality");
    /* Mean-field: α̅=0, β=0.5, γ=1 → 0 + 2*0.5 + 1 = 2 ≥ 2 */
    assert(oep_pt_check_rushbrooke(0.0, 0.5, 1.0));
    PASS();
}

void test_critical_slowing(void) {
    TEST("critical_slowing");
    /* τ ∝ 1/|p-pc|^νz */
    double tau = oep_pt_correlation_time(1.0, 0.1, 0.0, 0.5, 2.0);
    assert(tau > 1.0);  /* Should be larger than τ₀ */
    double tau_near = oep_pt_correlation_time(1.0, 0.01, 0.0, 0.5, 2.0);
    assert(tau_near > tau);  /* Closer to critical → slower */
    PASS();
}

/* ── L5: Emergence Metrics Tests ───────────────────────── */

void test_shannon_entropy(void) {
    TEST("shannon_entropy");
    double unif[] = {0.5, 0.5};
    double H = oep_em_entropy(unif, 2);
    assert(fabs(H - log(2.0)) < 1e-6);  /* H = ln 2 */
    double det[] = {1.0, 0.0};
    H = oep_em_entropy(det, 2);
    assert(fabs(H) < 1e-10);  /* Deterministic → H = 0 */
    PASS();
}

void test_mutual_information(void) {
    TEST("mutual_information");
    double joint[] = {0.25, 0.25, 0.25, 0.25};  /* Independent uniform */
    OEPMatrix* J = oep_mat_create_from(joint, 2, 2);
    double I = oep_em_mutual_information(J, 2, 2);
    assert(fabs(I) < 1e-10);  /* No mutual information */
    oep_mat_free(J);
    PASS();
}

void test_effective_information(void) {
    TEST("effective_information");
    /* Identity transition (deterministic) */
    double tid[] = {1.0, 0.0, 0.0, 1.0};
    OEPMatrix* T = oep_mat_create_from(tid, 2, 2);
    double EI = oep_em_effective_information(T, 2);
    assert(EI > 0.0);  /* Deterministic → high EI */
    oep_mat_free(T);
    PASS();
}

void test_timeseries_autocorrelation(void) {
    TEST("timeseries_autocorrelation");
    double ts[] = {1, 2, 1, 2, 1, 2, 1, 2, 1, 2};  /* Period 2 */
    double ac1 = oep_ts_autocorrelation(ts, 10, 1);
    assert(ac1 < 0.0);  /* Anti-correlated at lag 1 for period-2 */
    PASS();
}

/* ── L6: Collective Dynamics Tests ─────────────────────── */

void test_kuramoto_create(void) {
    TEST("kuramoto_create");
    OEPKuramotoSystem* sys = oep_kura_create(50);
    assert(sys != NULL);
    assert(sys->n_oscillators == 50);
    oep_kura_free(sys);
    PASS();
}

void test_kuramoto_order_parameter(void) {
    TEST("kuramoto_order_parameter");
    OEPKuramotoSystem* sys = oep_kura_create(10);
    assert(sys != NULL);
    /* Set all phases to 0 → perfect synchronization */
    for (int i = 0; i < 10; i++) sys->oscillators[i].phase = 0.0;
    oep_kura_compute_order_parameter(sys);
    assert(fabs(sys->mean_field_coherence - 1.0) < 1e-10);
    oep_kura_free(sys);
    PASS();
}

void test_vicsek_create(void) {
    TEST("vicsek_create");
    OEPVicsekSystem* sys = oep_vicsek_create(100, 10.0);
    assert(sys != NULL);
    assert(sys->n_agents == 100);
    oep_vicsek_free(sys);
    PASS();
}

void test_ising_create(void) {
    TEST("ising_create");
    OEPIsingSystem* sys = oep_ising_create(10, 10);
    assert(sys != NULL);
    assert(sys->n_spins == 100);
    oep_ising_free(sys);
    PASS();
}

void test_wilson_cowan(void) {
    TEST("wilson_cowan");
    OEPWilsonCowan* model = oep_wc_create();
    assert(model != NULL);
    oep_wc_set_default_params(model);
    double dE, dI;
    oep_wc_rhs(model, &dE, &dI);
    /* Derivatives should be finite */
    assert(isfinite(dE) && isfinite(dI));
    oep_wc_free(model);
    PASS();
}

/* ── L7: Landau-Ginzburg Tests ──────────────────────────── */

void test_field_create(void) {
    TEST("field_create");
    OEPSpatialField* field = oep_field_create(32, 32, 1.0, 1.0);
    assert(field != NULL);
    assert(field->nx == 32 && field->ny == 32);
    oep_field_free(field);
    PASS();
}

void test_field_laplacian(void) {
    TEST("field_laplacian");
    OEPSpatialField* field = oep_field_create(16, 16, 1.0, 1.0);
    assert(field != NULL);
    /* Set a peak at center */
    oep_field_set(field, 8, 8, 1.0);
    double lap = oep_field_laplacian(field, 8, 8);
    /* Laplacian of a peak should be negative (curvature < 0) */
    assert(lap < 0.0);
    oep_field_free(field);
    PASS();
}

void test_lg_create(void) {
    TEST("landau_ginzburg_create");
    OEPLandauGinzburg* lg = oep_lg_create(32, 32, 1.0, 1.0, -1.0, 1.0, 0.1, 1.0);
    assert(lg != NULL);
    oep_lg_free(lg);
    PASS();
}

void test_lg_correlation_length(void) {
    TEST("lg_correlation_length");
    /* ξ = √(c/|a|) = √(1.0/1.0) = 1.0 */
    double xi = oep_lg_correlation_length(-1.0, 1.0);
    assert(fabs(xi - 1.0) < 1e-10);
    /* Near critical: ξ = √(1.0/0.01) = 10 */
    double xi_crit = oep_lg_correlation_length(-0.01, 1.0);
    assert(fabs(xi_crit - 10.0) < 1e-6);
    PASS();
}

void test_superconductor_type(void) {
    TEST("superconductor_type");
    /* κ = 0.5 < 1/√2 → Type I */
    assert(!oep_lg_is_type_II(0.5));
    /* κ = 10 > 1/√2 → Type II */
    assert(oep_lg_is_type_II(10.0));
    PASS();
}

void test_domain_wall(void) {
    TEST("domain_wall_energy");
    /* σ = (2√2/3) * |a|^{3/2} * √c / b */
    double sigma = oep_lg_domain_wall_energy(-1.0, 1.0, 0.5);
    assert(sigma > 0.0);
    PASS();
}

/* ── Edge Cases / Robustness Tests ─────────────────────── */

void test_null_pointer_safety(void) {
    TEST("null_pointer_safety");
    oep_vec_free(NULL);
    oep_mat_free(NULL);
    assert(oep_vec_create(-1) == NULL);
    assert(oep_vec_create(0) == NULL);
    assert(oep_mat_create(0, 0) == NULL);
    oep_vec_print(NULL, "test");
    oep_mat_print(NULL, "test");
    PASS();
}

void test_zero_size_operations(void) {
    TEST("zero_size_operations");
    double det = oep_mat_determinant(NULL);
    assert(fabs(det) < 1e-10);
    int rank = oep_mat_rank_gaussian(NULL);
    assert(rank == 0);
    PASS();
}

void test_rng(void) {
    TEST("rng_basic");
    oep_rng_set_seed(42);
    double u1 = oep_rng_uniform();
    assert(u1 >= 0.0 && u1 <= 1.0);
    double n = oep_rng_normal(0.0, 1.0);
    assert(isfinite(n));
    PASS();
}

/* ── Application Tests ──────────────────────────────────── */

void test_enso_application(void) {
    TEST("enso_application");
    double sst[240] = {0}; /* 20 years of monthly data */
    for (int i = 0; i < 240; i++)
        sst[i] = 0.5 * sin(2.0 * OEP_PI * (double)i / 48.0) + oep_rng_normal(0.0, 0.2);
    OEPENSOAnalysis* enso = oep_app_enso_analyze(sst, 240);
    assert(enso != NULL);
    oep_app_enso_free(enso);
    PASS();
}

void test_ecosystem_resilience(void) {
    TEST("ecosystem_resilience");
    /* Generate auto-correlated data: AR(1) with φ=0.7 */
    double biomass[200];
    biomass[0] = 10.0;
    for (int i = 1; i < 200; i++)
        biomass[i] = 0.7 * biomass[i-1] + 3.0 + oep_rng_normal(0.0, 0.1);
    OEPEcosystemResilience* eco = oep_app_ecosystem_resilience_analyze(biomass, 200);
    assert(eco != NULL);
    /* With AR(1) φ=0.7, autocorrelation at lag 1 should be positive after detrending */
    assert(eco->autocorr_lag1 > 0.3);
    oep_app_ecosystem_free(eco);
    PASS();
}

/* ── Main ──────────────────────────────────────────────── */

int main(void) {
    printf("=== mini-order-parameters-emergence Test Suite ===\n\n");

    /* L1: Core Definitions */
    test_vector_create();
    test_matrix_create();
    test_vector_operations();
    test_matrix_multiplication();
    test_determinant();
    test_rank();
    test_null_pointer_safety();
    test_zero_size_operations();
    test_rng();

    /* L2: Order Parameter */
    test_scalar_op_create();
    test_order_parameter_equilibria_supercritical();
    test_order_parameter_equilibria_subcritical();
    test_stability_check();
    test_potential_symmetry();
    test_op_evolution();
    test_complex_op_create();

    /* L3: Slaving Principle */
    test_slaving_create();
    test_slaving_timescale_check();
    test_quasistationary_q();
    test_slaving_adiabatic_elimination();

    /* L4: Phase Transition */
    test_pitchfork_bifurcation();
    test_phase_transition_detection();
    test_landau_free_energy();
    test_rushbrooke_inequality();
    test_critical_slowing();

    /* L5: Emergence Metrics */
    test_shannon_entropy();
    test_mutual_information();
    test_effective_information();
    test_timeseries_autocorrelation();

    /* L6: Collective Dynamics */
    test_kuramoto_create();
    test_kuramoto_order_parameter();
    test_vicsek_create();
    test_ising_create();
    test_wilson_cowan();

    /* L7: Landau-Ginzburg */
    test_field_create();
    test_field_laplacian();
    test_lg_create();
    test_lg_correlation_length();
    test_superconductor_type();
    test_domain_wall();

    /* Applications */
    test_enso_application();
    test_ecosystem_resilience();

    printf("\n=== Results: %d/%d tests passed ===\n", tests_passed, tests_total);
    return (tests_passed == tests_total) ? 0 : 1;
}
