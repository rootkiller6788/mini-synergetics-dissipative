/**
 * @file test_core.c
 * @brief Tests for core entropy production and Onsager matrix operations.
 *
 * Tests cover:
 *   - L1: Structure initialization
 *   - L2: sigma computation from forces/fluxes
 *   - L3: Quadratic forms and Rayleigh function
 *   - L4: Prigogine theorem, Onsager reciprocity
 *   - L5: Matrix operations (LDL^T, solve, invert)
 */

#include "entropy_production.h"
#include "onsager_relations.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define TOL 1e-10

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) do { \
    tests_run++; \
    printf("  TEST %s... ", name); \
} while(0)

#define PASS() do { \
    tests_passed++; \
    printf("PASSED\n"); \
} while(0)

#define FAIL(msg) do { \
    printf("FAILED: %s\n", msg); \
} while(0)

/* ---------------------------------------------------------------------------
 * L1: Definitions -- Structure initialization
 * ------------------------------------------------------------------------- */

static void test_state_init(void)
{
    TEST("ThermodynamicState initialization");
    ThermodynamicState s;
    ep_state_init(&s, 300.0, 101325.0, 2);
    assert(fabs(s.temperature - 300.0) < TOL);
    assert(fabs(s.pressure - 101325.0) < TOL);
    assert(s.num_species == 2);
    PASS();
}

static void test_state_copy(void)
{
    TEST("ThermodynamicState copy");
    ThermodynamicState s1, s2;
    ep_state_init(&s1, 350.0, 200000.0, 3);
    s1.mole_numbers[0] = 1.0;
    s1.mole_numbers[1] = 2.0;
    ep_state_copy(&s2, &s1);
    assert(fabs(s2.temperature - 350.0) < TOL);
    assert(fabs(s2.mole_numbers[0] - 1.0) < TOL);
    PASS();
}

/* ---------------------------------------------------------------------------
 * L2: Core Concepts -- sigma computation
 * ------------------------------------------------------------------------- */

static void test_compute_sigma(void)
{
    TEST("Entropy production from forces and fluxes");
    double forces[] = {1.0, 2.0, 0.5};
    double fluxes[] = {3.0, 1.0, 4.0};
    double sigma = ep_compute_sigma(forces, fluxes, 3);
    /* sigma = 1*3 + 2*1 + 0.5*4 = 3 + 2 + 2 = 7 */
    assert(fabs(sigma - 7.0) < TOL);
    PASS();
}

static void test_compute_sigma_null(void)
{
    TEST("sigma with NULL pointers returns 0");
    double f[] = {1.0};
    double j[] = {1.0};
    double s = ep_compute_sigma(NULL, j, 1);
    assert(fabs(s) < TOL);
    s = ep_compute_sigma(f, NULL, 1);
    assert(fabs(s) < TOL);
    s = ep_compute_sigma(f, j, 0);
    assert(fabs(s) < TOL);
    PASS();
}

static void test_second_law_check(void)
{
    TEST("Second Law verification");
    assert(ep_check_second_law(5.0, 1e-12) == 1);
    assert(ep_check_second_law(0.0, 1e-12) == 1);
    assert(ep_check_second_law(-1e-13, 1e-12) == 1);   /* Within tolerance */
    assert(ep_check_second_law(-1e-10, 1e-12) == 0);   /* Violation */
    PASS();
}

/* ---------------------------------------------------------------------------
 * L3: Mathematical Structures -- Quadratic forms
 * ------------------------------------------------------------------------- */

static void test_quadratic_form(void)
{
    TEST("Quadratic form in forces (sigma = X^T L X)");
    /* Identity matrix: L = I */
    double L[] = {1, 0, 0, 1};
    double X[] = {3.0, 4.0};
    double sigma = ep_quadratic_form_forces(X, L, 2);
    /* sigma = 3^2 + 4^2 = 25 */
    assert(fabs(sigma - 25.0) < TOL);
    PASS();
}

static void test_rayleigh_dissipative(void)
{
    TEST("Rayleigh dissipative function");
    /* R = I, J = [1, 2] */
    double R[] = {1, 0, 0, 1};
    double J[] = {1.0, 2.0};
    double Phi = ep_rayleigh_dissipative(J, R, 2);
    /* Phi = 0.5 * (1*1 + 2*2) = 0.5 * 5 = 2.5 */
    assert(fabs(Phi - 2.5) < TOL);

    /* sigma from fluxes should equal 2*Phi */
    double sigma = ep_quadratic_form_fluxes(J, R, 2);
    assert(fabs(sigma - 5.0) < TOL);
    assert(fabs(sigma - 2.0 * Phi) < TOL);
    PASS();
}

static void test_ideal_gas_entropy(void)
{
    TEST("Ideal gas mixture entropy");
    ThermodynamicState s;
    ep_state_init(&s, 300.0, 1e5, 2);
    s.mole_numbers[0] = 1.0;
    s.mole_numbers[1] = 2.0;
    double std_S[] = {130.0, 200.0};  /* Standard molar entropies */
    double S = ep_ideal_gas_entropy(&s, std_S);
    /* Basic check: entropy should be positive */
    assert(S > 0.0);
    PASS();
}

/* ---------------------------------------------------------------------------
 * L4: Fundamental Laws -- Prigogine theorem
 * ------------------------------------------------------------------------- */

static void test_prigogine_dsigma_dt(void)
{
    TEST("Prigogine: d(sigma)/dt computation");
    /* L = I, X = [1, 2], dX/dt = [-1, -1] */
    double L[] = {1, 0, 0, 1};
    double X[] = {1.0, 2.0};
    double dXdt[] = {-1.0, -1.0};
    double dsigma = ep_prigogine_dsigma_dt(X, dXdt, L, 2);
    /* dsigma/dt = 2*(L11*X1*dX1 + L22*X2*dX2) = 2*(1*1*(-1) + 1*2*(-1)) = -6 */
    assert(fabs(dsigma - (-6.0)) < TOL);
    /* Negative => approaching steady state */
    assert(dsigma < 0.0);
    PASS();
}

static void test_mep_steady_state(void)
{
    TEST("MEP steady state check");
    /* L = I, X = [0, 0] is equilibrium (MEP with no constraints) */
    double L[] = {1, 0, 0, 1};
    double X[] = {0.0, 0.0};
    int constrained[] = {0, 0};
    int is_mep = ep_is_mep_steady_state(X, L, constrained, 2, 1e-10);
    assert(is_mep == 1);

    /* X = [1, 0] is NOT an MEP state with no constraints */
    double X2[] = {1.0, 0.0};
    is_mep = ep_is_mep_steady_state(X2, L, constrained, 2, 1e-10);
    assert(is_mep == 0);
    PASS();
}

static void test_find_mep_state(void)
{
    TEST("Find MEP steady state by optimization");
    /* L = I, X_0 constrained to 3.0, X_1 free.
     * Min sigma = X0^2 + X1^2 subject to X0=3 => X1=0, sigma=9 */
    double L[] = {1, 0, 0, 1};
    double X[] = {3.0, 10.0};  /* Initial guess */
    int constrained[] = {1, 0};
    double fixed_vals[] = {3.0, 0.0};
    int ret = ep_find_mep_steady_state(X, L, constrained, fixed_vals, 2, 1000, 1e-10);
    assert(ret == 0);
    assert(fabs(X[0] - 3.0) < TOL);
    assert(fabs(X[1]) < 1e-8);
    PASS();
}

/* ---------------------------------------------------------------------------
 * L4: Onsager reciprocity
 * ------------------------------------------------------------------------- */

static void test_onsager_alloc(void)
{
    TEST("Onsager matrix allocation");
    OnsagerMatrix *L = onsager_alloc(3);
    assert(L != NULL);
    assert(L->dim == 3);
    onsager_free(L);
    PASS();
}

static void test_onsager_set_get(void)
{
    TEST("Onsager matrix set/get");
    OnsagerMatrix *L = onsager_alloc(2);
    onsager_set(L, 0, 0, 10.0);
    onsager_set(L, 0, 1, 2.0);
    onsager_set(L, 1, 0, 2.0);
    onsager_set(L, 1, 1, 5.0);
    assert(fabs(onsager_get(L, 0, 0) - 10.0) < TOL);
    assert(fabs(onsager_get(L, 0, 1) - 2.0) < TOL);
    onsager_free(L);
    PASS();
}

static void test_onsager_reciprocity_verify(void)
{
    TEST("Onsager reciprocity verification");
    OnsagerMatrix *L = onsager_alloc(3);
    /* Symmetric matrix */
    onsager_set(L, 0, 0, 4.0);
    onsager_set(L, 0, 1, 1.0);
    onsager_set(L, 1, 0, 1.0);
    onsager_set(L, 0, 2, 2.0);
    onsager_set(L, 2, 0, 2.0);
    onsager_set(L, 1, 1, 3.0);
    onsager_set(L, 1, 2, 0.5);
    onsager_set(L, 2, 1, 0.5);
    onsager_set(L, 2, 2, 2.0);

    assert(onsager_verify_reciprocity(L, 1e-10) == 1);
    onsager_free(L);
    PASS();
}

static void test_onsager_diagonal(void)
{
    TEST("Diagonal Onsager matrix initialization");
    OnsagerMatrix *L = onsager_alloc(3);
    double diag[] = {2.0, 3.0, 1.0};
    onsager_init_diagonal(L, diag);
    assert(fabs(onsager_get(L, 0, 0) - 2.0) < TOL);
    assert(fabs(onsager_get(L, 0, 1)) < TOL);
    assert(fabs(onsager_get(L, 1, 1) - 3.0) < TOL);
    assert(fabs(onsager_get(L, 2, 2) - 1.0) < TOL);
    onsager_free(L);
    PASS();
}

static void test_onsager_two_process(void)
{
    TEST("Two-process coupled Onsager initialization");
    OnsagerMatrix *L = onsager_alloc(2);
    int ret = onsager_init_two_process(L, 4.0, 1.0, 1.5);
    assert(ret == 0);
    assert(fabs(onsager_get(L, 0, 0) - 4.0) < TOL);
    assert(fabs(onsager_get(L, 0, 1) - 1.5) < TOL);
    assert(fabs(onsager_get(L, 1, 0) - 1.5) < TOL);
    assert(fabs(onsager_get(L, 1, 1) - 1.0) < TOL);

    /* PSD check: 4*1 >= 1.5^2 = 2.25, yes */
    assert(onsager_verify_psd(L, 1e-10) == 1);
    onsager_free(L);
    PASS();
}

static void test_onsager_invalid_two_process(void)
{
    TEST("Invalid two-process (violates PSD)");
    OnsagerMatrix *L = onsager_alloc(2);
    /* L11*L22 = 1*1 = 1 < L12^2 = 4 -> not PSD */
    int ret = onsager_init_two_process(L, 1.0, 1.0, 2.0);
    assert(ret == -1);
    onsager_free(L);
    PASS();
}

/* ---------------------------------------------------------------------------
 * L5: Algorithms -- Linear algebra
 * ------------------------------------------------------------------------- */

static void test_onsager_ldlt(void)
{
    TEST("LDL^T decomposition");
    /* A = [4, 2; 2, 3] (positive definite) */
    double A[] = {4.0, 2.0, 2.0, 3.0};
    double D[4];
    int ret = onsager_ldlt_decompose(A, D, 2);
    assert(ret == 0);
    /* D should have positive diagonal */
    assert(D[0] > 0.0);
    assert(D[1] > 0.0);
    PASS();
}

static void test_onsager_solve(void)
{
    TEST("Solve L * x = b via LDL^T");
    double L[] = {4.0, 1.0, 1.0, 3.0};
    double b[] = {1.0, 2.0};
    double x[2];
    int ret = onsager_solve(L, b, x, 2);
    assert(ret == 0);
    /* L*x should equal b */
    double check0 = 4.0 * x[0] + 1.0 * x[1];
    double check1 = 1.0 * x[0] + 3.0 * x[1];
    assert(fabs(check0 - 1.0) < 1e-8);
    assert(fabs(check1 - 2.0) < 1e-8);
    PASS();
}

static void test_onsager_invert(void)
{
    TEST("Matrix inversion");
    double L[] = {2.0, 0.0, 0.0, 2.0};
    double R[4];
    int ret = onsager_invert(L, R, 2);
    assert(ret == 0);
    /* R should be [0.5, 0; 0, 0.5] */
    assert(fabs(R[0] - 0.5) < TOL);
    assert(fabs(R[1]) < TOL);
    assert(fabs(R[2]) < TOL);
    assert(fabs(R[3] - 0.5) < TOL);
    PASS();
}

static void test_onsager_fluxes(void)
{
    TEST("Compute fluxes from forces via Onsager matrix");
    OnsagerMatrix *L = onsager_alloc(2);
    onsager_set(L, 0, 0, 2.0);
    onsager_set(L, 0, 1, 0.5);
    onsager_set(L, 1, 0, 0.5);
    onsager_set(L, 1, 1, 1.0);

    double X[] = {3.0, 4.0};
    double J[2];
    onsager_compute_fluxes(L, X, J);
    /* J0 = 2*3 + 0.5*4 = 6 + 2 = 8
     * J1 = 0.5*3 + 1*4 = 1.5 + 4 = 5.5 */
    assert(fabs(J[0] - 8.0) < TOL);
    assert(fabs(J[1] - 5.5) < TOL);

    /* sigma = J dot X = 8*3 + 5.5*4 = 24 + 22 = 46 */
    double sigma = ep_quadratic_form_forces(X, L->data, 2);
    assert(fabs(sigma - 46.0) < TOL);

    onsager_free(L);
    PASS();
}

static void test_onsager_eigenvalues(void)
{
    TEST("Onsager matrix eigenvalues");
    OnsagerMatrix *L = onsager_alloc(2);
    onsager_set(L, 0, 0, 4.0);
    onsager_set(L, 0, 1, 1.0);
    onsager_set(L, 1, 0, 1.0);
    onsager_set(L, 1, 1, 3.0);

    double evals[2];
    onsager_eigenvalues(L, evals, 100);
    /* Eigenvalues of [4 1; 1 3]:
     * trace = 7, det = 4*3 - 1*1 = 11
     * lambda = (7 +/- sqrt(49 - 44))/2 = (7 +/- sqrt(5))/2 */
    double expected1 = (7.0 + sqrt(5.0)) / 2.0;   /* ~ 4.618 */
    double expected2 = (7.0 - sqrt(5.0)) / 2.0;   /* ~ 2.382 */
    assert(fabs(evals[0] - expected1) < 1e-6);
    assert(fabs(evals[1] - expected2) < 1e-6);
    onsager_free(L);
    PASS();
}

static void test_reaction_affinity(void)
{
    TEST("Chemical reaction affinity");
    ChemicalReaction rxn;
    memset(&rxn, 0, sizeof(rxn));
    /* A + B -> C: nu_A=-1, nu_B=-1, nu_C=+1 */
    rxn.stoichiometry[0] = -1.0;
    rxn.stoichiometry[1] = -1.0;
    rxn.stoichiometry[2] = +1.0;
    rxn.num_species_involved = 3;

    double mu[] = {-100.0, -50.0, -200.0};  /* Chemical potentials */
    double A = ep_reaction_affinity(&rxn, mu, 3);
    /* A = -(-1*(-100) + -1*(-50) + 1*(-200))
     *   = -(100 + 50 - 200) = -(-50) = 50 */
    assert(fabs(A - 50.0) < TOL);
    PASS();
}

/* ---------------------------------------------------------------------------
 * Main
 * ------------------------------------------------------------------------- */

int main(void)
{
    printf("=== Core Tests: mini-entropy-production-minimization ===\n\n");

    test_state_init();
    test_state_copy();
    test_compute_sigma();
    test_compute_sigma_null();
    test_second_law_check();
    test_quadratic_form();
    test_rayleigh_dissipative();
    test_ideal_gas_entropy();
    test_prigogine_dsigma_dt();
    test_mep_steady_state();
    test_find_mep_state();
    test_onsager_alloc();
    test_onsager_set_get();
    test_onsager_reciprocity_verify();
    test_onsager_diagonal();
    test_onsager_two_process();
    test_onsager_invalid_two_process();
    test_onsager_ldlt();
    test_onsager_solve();
    test_onsager_invert();
    test_onsager_fluxes();
    test_onsager_eigenvalues();
    test_reaction_affinity();

    printf("\n=== Results: %d/%d tests passed ===\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}