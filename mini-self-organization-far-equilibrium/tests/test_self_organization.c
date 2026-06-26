#include "sofe_core.h"
#include "sofe_thermodynamics.h"
#include "sofe_dynamics.h"
#include "sofe_patterns.h"
#include "sofe_bifurcation.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

static int tests_run = 0, tests_passed = 0, tests_failed = 0;

#define TEST(n) do { tests_run++; printf("  TEST %s... ", n); } while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)
#define FAIL(m) do { tests_failed++; printf("FAIL: %s\n", m); } while(0)
#define CHECK(c, m) do { if (c) PASS(); else FAIL(m); } while(0)

static void test_system_create(void) {
    TEST("system_create");
    SOFESystem* sys = sofe_system_create("test_sys");
    CHECK(sys != NULL, "sys null");
    CHECK(sys->temperature == 298.15, "default T");
    sofe_system_free(sys);
}

static void test_network_create(void) {
    TEST("network_create");
    SOFEReactionNetwork* net = sofe_network_create("test_net");
    int a = sofe_network_add_species(net, "A", 1.0, 0.01);
    int b = sofe_network_add_species(net, "B", 0.5, 0.02);
    CHECK(a == 0 && b == 1, "species idx");
    CHECK(net->n_species == 2, "species cnt");
    int r = sofe_network_add_reaction(net, "R1", 1.0, 0.0);
    int reactants[] = {0}; double stoich[] = {1.0};
    sofe_network_set_reactants(net, r, reactants, stoich, 1);
    sofe_network_compute_stoichiometric_matrix(net);
    CHECK(net->stoichiometric_matrix != NULL, "stoich");
    sofe_network_compute_jacobian(net);
    CHECK(net->jacobian != NULL, "jacobian");
    sofe_network_free(net);
}

static void test_reaction_rate(void) {
    TEST("reaction_rate");
    SOFEReactionNetwork* net = sofe_network_create("rt");
    int x = sofe_network_add_species(net, "X", 2.0, 0.01);
    int y = sofe_network_add_species(net, "Y", 1.0, 0.02);
    int r = sofe_network_add_reaction(net, "R", 1.5, 0.0);
    int reactants[] = {x, y}; double stoich[] = {1.0, 1.0};
    sofe_network_set_reactants(net, r, reactants, stoich, 2);
    double rate = sofe_network_reaction_rate(net, r);
    CHECK(fabs(rate - 3.0) < 0.01, "rate=3.0");
    sofe_network_free(net);
}

static void test_grid_operations(void) {
    TEST("grid_ops");
    SOFEGrid* g1d = sofe_grid_create_1d(10, 1.0, true);
    CHECK(g1d != NULL && g1d->n_points == 10, "1d grid");
    sofe_grid_free(g1d);
    SOFEGrid* g2d = sofe_grid_create_2d(5, 4, 1.0, 1.0, false, false);
    CHECK(g2d != NULL && g2d->n_points == 20, "2d grid");
    int nb[12]; int nf = 0;
    sofe_grid_get_neighbors(g2d, 6, nb, &nf);
    CHECK(nf >= 2, "neighbors");
    sofe_grid_free(g2d);
}

static void test_stability_classification(void) {
    TEST("stability");
    double stable[] = {-1.0, 0.0, -2.0, 0.0};
    CHECK(sofe_classify_stability(stable, 2) == SOFE_STABLE_NODE, "stable node");
    double saddle[] = {-1.0, 0.0, 1.0, 0.0};
    CHECK(sofe_classify_stability(saddle, 2) == SOFE_SADDLE, "saddle");
}

static void test_entropy_production(void) {
    TEST("entropy_prod");
    SOFEForceFluxSystem ffs;
    ffs.n_pairs = 2;
    ffs.pairs = (SOFEForceFluxPair*)malloc(2 * sizeof(SOFEForceFluxPair));
    ffs.pairs[0].force = 1.0; ffs.pairs[0].flux = 0.5;
    ffs.pairs[1].force = 2.0; ffs.pairs[1].flux = 0.3;
    double sigma = sofe_entropy_production_direct(&ffs);
    CHECK(sigma >= 0.0, "nonnegative");
    CHECK(fabs(sigma - 1.1) < 1e-10, "value");
    free(ffs.pairs);
}

static void test_onsager(void) {
    TEST("onsager");
    SOFEOnsagerMatrix ons;
    sofe_onsager_create(&ons, 3);
    sofe_onsager_set_coefficient(&ons, 0, 1, 2.0);
    sofe_onsager_set_coefficient(&ons, 1, 0, 2.0);
    CHECK(sofe_onsager_verify_reciprocity(&ons, 1e-10), "Lij=Lji");
    sofe_onsager_free(&ons);
}

static void test_fluctuation_theorems(void) {
    TEST("fluct_theorems");
    double kT = 2.5, dF = 1.0;
    double ratio = sofe_crooks_ratio(1.5, dF, kT);
    CHECK(ratio > 0.0, "Crooks");
    double works[] = {1.0, 1.2, 0.8, 1.5, 1.1};
    double dF_est = sofe_jarzynski_estimator(works, 5, kT);
    CHECK(dF_est > 0.0, "Jarzynski");
}

static void test_model_rates(void) {
    TEST("model_rates");
    double y[3], p[8], dydt[3];
    y[0]=1.0; y[1]=1.0; p[0]=1.0; p[1]=3.0;
    sofe_brusselator_rates(y, 0.0, p, 2, dydt);
    CHECK(fabs(dydt[0] - (-2.0)) < 1e-10, "Brusselator");
    sofe_lotka_volterra_rates(y, 0.0, p, 4, dydt);
    p[0]=1.0; p[1]=1.0; p[2]=1.0;
    sofe_ginzburg_landau_rates(y, 0.0, p, 3, dydt);
    CHECK(fabs(dydt[0] - (-2.0)) < 1e-10, "GL");
}

static void test_integration(void) {
    TEST("integration");
    double y[] = {1.0};
    double p[] = {0.0};
    sofe_integrate_step(y, 1, 0.0, 0.1, sofe_ginzburg_landau_rates, p, 3, SOFE_INT_RK4);
    CHECK(fabs(y[0]) < 10.0, "RK4 step bounded");
}

static void test_pattern_detection(void) {
    TEST("pattern");
    SOFEGrid* g = sofe_grid_create_1d(10, 1.0, false);
    SOFEField* f = sofe_field_create(g, 1);
    for (int i = 0; i < 10; i++) sofe_field_set(f, 0, i, sin(i * 0.5));
    SOFEPattern* pat = sofe_pattern_detect(f);
    CHECK(pat != NULL, "pattern obj");
    CHECK(pat->type == SOFE_PATTERN_STRIPES, "stripes");
    sofe_pattern_free(pat);
    sofe_field_free(f); sofe_grid_free(g);
}

int main(void) {
    printf("=== SOFE Test Suite ===\n");
    test_system_create();
    test_network_create();
    test_reaction_rate();
    test_grid_operations();
    test_stability_classification();
    test_entropy_production();
    test_onsager();
    test_fluctuation_theorems();
    test_model_rates();
    test_integration();
    test_pattern_detection();
    printf("\n=== %d/%d passed, %d failed ===\n", tests_passed, tests_run, tests_failed);
    return tests_failed > 0 ? 1 : 0;
}
