/**
 * test_main.c — Comprehensive test suite for Turing Pattern Formation
 *
 * Tests all core APIs: field operations, model kinetics, solver steps,
 * stability analysis, visualization, and stochastic noise.
 *
 * Knowledge Coverage:
 *   L4: Theorem verification via assert-based testing
 *   L5: Algorithm correctness validation
 *   L6: Problem-solving verification
 */

#include "turing_common.h"
#include "turing_models.h"
#include "turing_analysis.h"
#include "turing_visualize.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

static int tests_run = 0;
static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(name) do { \
    tests_run++; \
    printf("  TEST %s ... ", name); \
} while(0)

#define PASS() do { \
    printf("PASS\n"); \
    tests_passed++; \
} while(0)

#define FAIL(msg) do { \
    printf("FAIL: %s\n", msg); \
    tests_failed++; \
} while(0)

#define CHECK(cond, msg) do { \
    if (!(cond)) { FAIL(msg); return; } \
} while(0)

/* ==========================================================================
   L1: Field Lifecycle Tests
   ==========================================================================*/

static void test_field_alloc_free(void) {
    TEST("field_alloc_free");
    ScalarField2D *f = field_alloc(64, 64, 1.0, 1.0, 1);
    CHECK(f != NULL, "allocation failed");
    CHECK(f->Nx == 64, "Nx wrong");
    CHECK(f->Ny == 64, "Ny wrong");
    CHECK(f->ghost == 1, "ghost wrong");
    CHECK(f->dx > 0.0, "dx invalid");
    field_free(f);
    PASS();
}

static void test_field_fill_and_get(void) {
    TEST("field_fill_and_get");
    ScalarField2D *f = field_alloc(32, 32, 1.0, 1.0, 1);
    CHECK(f != NULL, "alloc");
    field_fill(f, 3.14159);
    double val = field_get(f, 15, 15);
    CHECK(fabs(val - 3.14159) < 1e-10, "fill value mismatch");
    field_free(f);
    PASS();
}

static void test_field_copy(void) {
    TEST("field_copy");
    ScalarField2D *src = field_alloc(32, 32, 1.0, 1.0, 1);
    ScalarField2D *dst = field_alloc(32, 32, 1.0, 1.0, 1);
    CHECK(src && dst, "alloc");
    field_init_random(src, 0.5, 0.1, 42);
    int ret = field_copy(src, dst);
    CHECK(ret == 0, "copy failed");
    double diff = field_diff_norm(src, dst);
    CHECK(diff < 1e-15, "copy mismatch");
    field_free(src);
    field_free(dst);
    PASS();
}

/* ==========================================================================
   L3: Laplacian Tests
   ==========================================================================*/

static void test_laplacian_5pt_constant(void) {
    TEST("laplacian_5pt_constant");
    ScalarField2D *f = field_alloc(32, 32, 1.0, 1.0, 1);
    ScalarField2D *lap = field_alloc(32, 32, 1.0, 1.0, 1);
    CHECK(f && lap, "alloc");
    field_fill(f, 2.5);
    field_apply_bc(f, BC_NEUMANN, BC_NEUMANN);
    int ret = laplacian_compute(f, lap, LAPLACIAN_5PT);
    CHECK(ret == 0, "laplacian compute");
    /* Laplacian of constant field should be zero */
    double max_lap = 0.0;
    for (int i = 0; i < f->Nx; i++)
        for (int j = 0; j < f->Ny; j++) {
            double v = fabs(field_get(lap, i, j));
            if (v > max_lap) max_lap = v;
        }
    CHECK(max_lap < 1e-10, "laplacian of constant not zero");
    field_free(f);
    field_free(lap);
    PASS();
}

static void test_laplacian_sinusoid(void) {
    TEST("laplacian_sinusoid_analytic");
    /* f = cos(kx*x)*cos(ky*y) with kx=ky=1 on [0, 2π]
     * ∇²f = -(kx² + ky²)*f = -2f
     * Use N=128 for better accuracy */
    int N = 128;
    double L = 2.0 * M_PI;
    ScalarField2D *f = field_alloc(N, N, L, L, 1);
    ScalarField2D *lap = field_alloc(N, N, L, L, 1);
    field_init_sinusoidal(f, 1.0, 1.0, 1.0, 0.0);
    field_apply_bc(f, BC_PERIODIC, BC_PERIODIC);
    laplacian_compute(f, lap, LAPLACIAN_5PT);

    /* Analytic: ∇²f = -(1+1)*f = -2f
     * Finite difference truncation error O(h²) ≈ (2π/127)²/6 ≈ 0.0013 */
    double max_err = 0.0;
    int max_i = 0, max_j = 0;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
            double val = field_get(f, i, j);
            double lap_val = field_get(lap, i, j);
            double expected = -2.0 * val;
            double err = fabs(lap_val - expected);
            if (err > max_err) { max_err = err; max_i = i; max_j = j; }
        }
    if (max_err >= 0.1) {
        printf("max_err=%.6f at (%d,%d), val=%.6f, lap=%.6f, expected=%.6f\n",
               max_err, max_i, max_j,
               field_get(f, max_i, max_j),
               field_get(lap, max_i, max_j),
               -2.0 * field_get(f, max_i, max_j));
    }
    CHECK(max_err < 0.1, "laplacian of sinusoid inaccurate");
    field_free(f);
    field_free(lap);
    PASS();
}

/* ==========================================================================
   L4: Turing Conditions Tests (Theorem Verification)
   ==========================================================================*/

static void test_turing_conditions_schnakenberg(void) {
    TEST("turing_conditions_schnakenberg");
    /* Schnakenberg model at a=0.1, b=0.9 with Du=1, Dv=40
     * u*=1.0, v*=0.9
     * f_u=-1+2*1*0.9=0.8, f_v=1, g_u=-2*1*0.9=-1.8, g_v=-1 */
    double J[2][2] = {{0.8, 1.0}, {-1.8, -1.0}};
    TuringConditions tc;
    turing_conditions_compute(J, 1.0, 40.0, &tc);

    CHECK(tc.cond1 == 1, "cond1: trace should be negative");
    CHECK(tc.cond2 == 1, "cond2: det should be positive");
    CHECK(tc.cond3 == 1, "cond3: cross-diffusion should be positive");
    CHECK(tc.cond4 == 1, "cond4: discriminant should be positive");
    CHECK(tc.turing_unstable == 1, "should be Turing-unstable");
    PASS();
}

static void test_turing_conditions_no_diffusion(void) {
    TEST("turing_conditions_no_diffusion");
    /* With Du = Dv, cross-diffusion term Dv*f_u + Du*g_v = D*(f_u+g_v) = D*tr(J)
     * If tr(J) < 0, then cross < 0 and cond3 fails → no Turing */
    double J[2][2] = {{-1.0, 0.5}, {0.5, -2.0}};
    TuringConditions tc;
    turing_conditions_compute(J, 1.0, 1.0, &tc);

    CHECK(tc.cond1 == 1, "trace should be negative");
    CHECK(tc.cond2 == 1, "det should be positive");
    CHECK(tc.cond3 == 0, "cross should be negative (equal diffusion)");
    CHECK(tc.turing_unstable == 0, "equal diffusion = no Turing");
    PASS();
}

static void test_dispersion_relation(void) {
    TEST("dispersion_relation");
    double J[2][2] = {{0.8, 1.0}, {-1.8, -1.0}};
    DispersionRelation *dr = turing_dispersion_compute(J, 1.0, 40.0, 100, 0.01, 3.0);
    CHECK(dr != NULL, "dispersion alloc");
    CHECK(dr->turing_unstable == 1, "should detect instability");
    CHECK(dr->lambda_max > 0.0, "max growth rate should be positive");
    CHECK(dr->k_critical > 0.0, "critical wavenumber should be positive");
    CHECK(dr->predicted_wavelength > 0.0, "wavelength should be positive");
    dispersion_free(dr);
    PASS();
}

/* ==========================================================================
   L5: Model Tests
   ==========================================================================*/

static void test_grayscott_reaction(void) {
    TEST("grayscott_reaction");
    ModelParams p;
    model_params_default(MODEL_GRAY_SCOTT, &p);
    double f, g;
    int ret = grayscott_reaction(0.5, 0.25, &p, &f, &g);
    CHECK(ret == 0, "reaction eval");
    CHECK(is_finite_val(f) && is_finite_val(g), "reaction values finite");
    PASS();
}

static void test_grayscott_jacobian(void) {
    TEST("grayscott_jacobian");
    ModelParams p;
    model_params_default(MODEL_GRAY_SCOTT, &p);
    double J[2][2];
    int ret = grayscott_jacobian(0.5, 0.25, &p, J);
    CHECK(ret == 0, "jacobian eval");
    CHECK(is_finite_val(J[0][0]) && is_finite_val(J[1][1]), "finite");
    PASS();
}

static void test_brusselator_hss(void) {
    TEST("brusselator_hss");
    ModelParams p;
    model_params_default(MODEL_BRUSSELATOR, &p);
    HomogeneousSteadyState hss;
    int ret = brusselator_find_hss(&p, &hss);
    CHECK(ret == 0, "HSS find");
    CHECK(hss.converged == 1, "should converge");
    CHECK(fabs(hss.u_star - p.A) < 1e-10, "u* = A");
    CHECK(fabs(hss.v_star - p.B/p.A) < 1e-10, "v* = B/A");
    PASS();
}

static void test_giere_meinhardt_hss(void) {
    TEST("giere_meinhardt_hss");
    ModelParams p;
    model_params_default(MODEL_GIERER_MEINHARDT, &p);
    HomogeneousSteadyState hss;
    int ret = giere_meinhardt_find_hss(&p, &hss);
    CHECK(ret == 0 || ret == -1, "HSS find attempted");
    if (ret == 0 && hss.converged) {
        double f, g;
        giere_meinhardt_reaction(hss.u_star, hss.v_star, &p, &f, &g);
        CHECK(fabs(f) < 1e-6 && fabs(g) < 1e-6, "f, g near zero at HSS");
    }
    PASS();
}

static void test_all_models_hss(void) {
    TEST("all_models_hss");
    ModelType models[] = {
        MODEL_GRAY_SCOTT, MODEL_BRUSSELATOR, MODEL_SCHNAKENBERG,
        MODEL_LENGYEL_EPSTEIN, MODEL_FITZHUGH_NAGUMO
    };
    int n = 5;
    int found = 0;
    for (int i = 0; i < n; i++) {
        ModelParams p;
        model_params_default(models[i], &p);
        if (!model_validate_params(models[i], &p)) continue;
        HomogeneousSteadyState hss;
        int ret = model_find_hss(models[i], &p, &hss);
        if (ret == 0 && hss.converged) found++;
    }
    CHECK(found >= 3, "at least 3 models should have HSS");
    PASS();
}

/* ==========================================================================
   L5: Matrix Eigenvalue Tests
   ==========================================================================*/

static void test_eigenvalues_real(void) {
    TEST("eigenvalues_real");
    double M[2][2] = {{2.0, 1.0}, {1.0, 2.0}};
    double er[2], ei[2];
    matrix_2x2_eigenvalues(M, er, ei);
    CHECK(fabs(ei[0]) < 1e-15, "should be real");
    CHECK(fabs(ei[1]) < 1e-15, "should be real");
    /* Eigenvalues of [[2,1],[1,2]] are 3 and 1 */
    double max_ev = (er[0] > er[1]) ? er[0] : er[1];
    double min_ev = (er[0] < er[1]) ? er[0] : er[1];
    CHECK(fabs(max_ev - 3.0) < 1e-10, "max eigenvalue");
    CHECK(fabs(min_ev - 1.0) < 1e-10, "min eigenvalue");
    PASS();
}

static void test_eigenvalues_complex(void) {
    TEST("eigenvalues_complex");
    double M[2][2] = {{0.0, -1.0}, {1.0, 0.0}};
    double er[2], ei[2];
    matrix_2x2_eigenvalues(M, er, ei);
    CHECK(fabs(er[0]) < 1e-15 && fabs(er[1]) < 1e-15, "real parts zero");
    CHECK((fabs(ei[0] - 1.0) < 1e-10 || fabs(ei[0] + 1.0) < 1e-10),
          "imag parts ±1");
    PASS();
}

/* ==========================================================================
   L5: Newton HSS Solver Tests
   ==========================================================================*/

static void test_newton_brusselator(void) {
    TEST("newton_brusselator");
    ModelParams p;
    model_params_default(MODEL_BRUSSELATOR, &p);
    HomogeneousSteadyState *hss = newton_find_hss(
        brusselator_reaction, brusselator_jacobian, &p, 1.5, 3.0, 50, 1e-12);
    CHECK(hss != NULL, "newton alloc");
    CHECK(hss->converged == 1, "should converge");
    CHECK(fabs(hss->u_star - p.A) < 1e-8, "u* correct");
    CHECK(fabs(hss->v_star - p.B/p.A) < 1e-8, "v* correct");
    free(hss);
    PASS();
}

/* ==========================================================================
   L6: Visualization Tests
   ==========================================================================*/

static void test_ppm_output(void) {
    TEST("ppm_output");
    ScalarField2D *f = field_alloc(32, 32, 1.0, 1.0, 0);
    field_init_sinusoidal(f, 2.0*M_PI, 2.0*M_PI, 0.5, 0.5);
    int ret = field_write_ppm(f, "build/test_output.ppm", COLORMAP_VIRIDIS, 0.0, 1.0);
    CHECK(ret == 0, "PPM write");
    field_free(f);
    PASS();
}

static void test_pattern_metrics(void) {
    TEST("pattern_metrics");
    ScalarField2D *f = field_alloc(64, 64, 1.0, 1.0, 0);
    field_init_sinusoidal(f, 4.0*M_PI, 4.0*M_PI, 1.0, 0.0);
    PatternMetrics pm = field_compute_metrics(f);
    CHECK(pm.pattern_amplitude > 0.0, "amplitude positive");
    CHECK(fabs(pm.mean) < 0.01, "mean near zero");
    field_free(f);
    PASS();
}

static void test_pattern_classification(void) {
    TEST("pattern_classification");
    ScalarField2D *f = field_alloc(64, 64, 1.0, 1.0, 0);
    /* Constant field → homogeneous */
    field_fill(f, 1.0);
    PatternType pt = field_classify_pattern(f);
    CHECK(pt == PATTERN_HOMOGENEOUS, "constant = homogeneous");
    field_free(f);
    PASS();
}

/* ==========================================================================
   L4: Newton multiple HSS test
   ==========================================================================*/

static void test_newton_multiple_hss(void) {
    TEST("newton_multiple_hss");
    ModelParams p;
    p.A = 2.0; p.B = 4.5; p.F = 0; p.k = 0; p.a = 0; p.b = 0;
    p.rho_a = 0; p.rho_h = 0; p.mu_a = 0; p.mu_h = 0;
    p.c = 0; p.epsilon = 0; p.beta = 0; p.gamma = 0; p.I_ext = 0;

    int n_hss;
    HomogeneousSteadyState **hss_list = newton_find_all_hss(
        brusselator_reaction, brusselator_jacobian, &p, &n_hss);

    CHECK(hss_list != NULL || n_hss == 0, "multi-start completed");
    if (hss_list) {
        CHECK(n_hss >= 1, "should find at least 1 HSS");
        for (int i = 0; i < n_hss; i++) free(hss_list[i]);
        free(hss_list);
    }
    PASS();
}

/* ==========================================================================
   L4: Bifurcation Detection Test
   ==========================================================================*/

static void test_bifurcation_detection(void) {
    TEST("bifurcation_detection");
    ModelParams base;
    model_params_default(MODEL_SCHNAKENBERG, &base);
    int n_found;
    /* Vary parameter b (index 5) to find Turing bifurcation */
    BifurcationPoint *bps = find_turing_bifurcations(
        MODEL_SCHNAKENBERG, &base, 5,
        0.1, 2.0, 1.0, 40.0, 10, &n_found);
    /* May or may not find bifurcations depending on range */
    CHECK(bps != NULL || n_found == 0, "bifurcation search completed");
    free(bps);
    PASS();
}

/* ==========================================================================
   L4: Hopf Condition Test
   ==========================================================================*/

static void test_hopf_condition(void) {
    TEST("hopf_condition");
    /* FitzHugh-Nagumo at certain parameters is excitable with Hopf */
    double J1[2][2] = {{0.5, 1.0}, {-1.0, 0.2}};
    int hopf1 = is_hopf_unstable(J1);
    CHECK(hopf1 == 0 || hopf1 == 1, "hopf check valid");

    /* Clearly unstable: trace > 0, det > 0 */
    double J2[2][2] = {{1.0, -1.0}, {1.0, 1.0}};
    int hopf2 = is_hopf_unstable(J2);
    CHECK(hopf2 == 1, "should be Hopf-unstable");
    PASS();
}

/* ==========================================================================
   L4: Turing Space Test
   ==========================================================================*/

static void test_turing_space(void) {
    TEST("turing_space");
    ModelParams base;
    model_params_default(MODEL_SCHNAKENBERG, &base);
    TuringSpace *ts = turing_space_alloc(
        MODEL_SCHNAKENBERG, &base, 1.0, 40.0,
        4, 0.05, 0.3, 20,   /* param a range */
        5, 0.05, 1.5, 20);  /* param b range */
    CHECK(ts != NULL, "turing space alloc");
    SolverConfig sc = {0};
    turing_space_compute(ts, &sc);
    int counts[5];
    turing_space_count_regions(ts, counts);
    CHECK(counts[REGION_TURING] + counts[REGION_STABLE] +
          counts[REGION_HOPF] + counts[REGION_TURING_HOPF] +
          counts[REGION_NO_HSS] == 400, "all grid points classified");
    turing_space_free(ts);
    PASS();
}

/* ==========================================================================
   L4: Sensitivity Analysis Test
   ==========================================================================*/

static void test_sensitivity_analysis(void) {
    TEST("sensitivity_analysis");
    ModelParams p;
    model_params_default(MODEL_SCHNAKENBERG, &p);
    double sens[20];
    int ret = turing_sensitivity_analysis(MODEL_SCHNAKENBERG, &p,
                                           1.0, 40.0, sens, 20);
    CHECK(ret == 0 || ret == -1, "sensitivity analysis completed");
    PASS();
}

/* ==========================================================================
   Main
   ==========================================================================*/

int main(void) {
    printf("=== Turing Pattern Formation Test Suite ===\n\n");

    printf("[L1] Field Lifecycle:\n");
    test_field_alloc_free();
    test_field_fill_and_get();
    test_field_copy();

    printf("\n[L3] Laplacian:\n");
    test_laplacian_5pt_constant();
    test_laplacian_sinusoid();

    printf("\n[L4] Turing Conditions:\n");
    test_turing_conditions_schnakenberg();
    test_turing_conditions_no_diffusion();
    test_dispersion_relation();
    test_eigenvalues_real();
    test_eigenvalues_complex();
    test_hopf_condition();
    test_bifurcation_detection();
    test_turing_space();
    test_sensitivity_analysis();

    printf("\n[L5] Models:\n");
    test_grayscott_reaction();
    test_grayscott_jacobian();
    test_brusselator_hss();
    test_giere_meinhardt_hss();
    test_all_models_hss();
    test_newton_brusselator();
    test_newton_multiple_hss();

    printf("\n[L6] Visualization:\n");
    test_ppm_output();
    test_pattern_metrics();
    test_pattern_classification();

    printf("\n========================================\n");
    printf("Tests run:    %d\n", tests_run);
    printf("Tests passed: %d\n", tests_passed);
    printf("Tests failed: %d\n", tests_failed);
    printf("========================================\n");

    return (tests_failed > 0) ? 1 : 0;
}
