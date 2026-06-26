/**
 * @file test_steady_state.c
 * @brief Tests for steady-state analysis and MEP optimization.
 */

#include "steady_state.h"
#include "onsager_relations.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define TOL 1e-8

int main(void)
{
    printf("=== Steady State Tests ===\n\n");

    /* Test ss_minimum_sigma */
    OnsagerMatrix *L = onsager_alloc(2);
    onsager_set(L, 0, 0, 2.0);
    onsager_set(L, 0, 1, 0.5);
    onsager_set(L, 1, 0, 0.5);
    onsager_set(L, 1, 1, 3.0);

    double X[] = {1.0, 2.0};
    double s = ss_minimum_sigma(L, X, 2);
    /* sigma = 2*1^2 + 2*0.5*1*2 + 3*2^2 = 2 + 2 + 12 = 16 */
    assert(fabs(s - 16.0) < TOL);
    printf("  Minimum sigma: PASSED\n");

    /* Test steady fluxes */
    double J[2];
    ss_compute_steady_fluxes(L, X, J, 2);
    /* J0 = 2*1 + 0.5*2 = 3, J1 = 0.5*1 + 3*2 = 6.5 */
    assert(fabs(J[0] - 3.0) < TOL);
    assert(fabs(J[1] - 6.5) < TOL);
    printf("  Steady fluxes: PASSED\n");

    /* Test MEP with constraints */
    SteadyStateConstraints c;
    ss_init_constraints(&c, 2);
    ss_add_constraint(&c, 0, 3.0);  /* Fix X0 = 3.0 */

    double X_sol[] = {3.0, 0.0};
    int ret = ss_find_mep_state(L, &c, X_sol, 2, 1e-10);
    assert(ret == 0);
    assert(fabs(X_sol[0] - 3.0) < TOL);
    /* With X0=3 constrained, min sigma requires X1 = -L10/L11 * X0 = -0.5/3 * 3 = -0.5 */
    /* Actually: partial sigma / partial X1 = 2*(L10*X0 + L11*X1) = 2*(0.5*3 + 3*X1) = 0
     * X1 = -(0.5*3)/3 = -0.5 */
    assert(fabs(X_sol[1] - (-0.5)) < 1e-6);
    printf("  MEP with constraints: PASSED\n");

    /* Test linear stability */
    StabilityResult stab;
    double X_ss[] = {0.0, 0.0};
    ss_linear_stability(&stab, X_ss, L, 1.0, 2);
    /* Eigenvalues of L: (2+3 +/- sqrt((2-3)^2 + 4*0.5^2))/2
     * trace=5, det=6-0.25=5.75
     * lambda = (5 +/- sqrt(25-23))/2 = (5 +/- sqrt(2))/2 ~ 3.207, 1.793
     * Jacobian eigenvalues = -2*gamma*lambda = -6.414, -3.586 */
    assert(stab.is_asymptotically_stable == 1);
    assert(stab.max_real_part < 0.0);
    printf("  Linear stability: PASSED\n");

    /* Test excess entropy production */
    double X_dev[] = {1.0, 0.5};
    double X_ss2[] = {0.0, 0.0};
    double delta2 = ss_excess_entropy_production(X_dev, X_ss2, L, 2);
    /* delta^2 = X_dev^T L X_dev = 2*1^2 + 2*0.5*1*0.5 + 3*0.5^2 = 2 + 0.5 + 0.75 = 3.25 */
    assert(fabs(delta2 - 3.25) < TOL);
    assert(delta2 >= 0.0);
    printf("  Excess entropy production: PASSED\n");

    /* Test Glansdorff-Prigogine criterion */
    double dXdt[] = {-0.5, -0.3};
    double gp = ss_glansdorff_prigogine_criterion(X_dev, dXdt, L, 2);
    /* J = L*X_dev = [2*1+0.5*0.5=2.25, 0.5*1+3*0.5=2.0]
     * gp = J . dXdt = 2.25*(-0.5) + 2.0*(-0.3) = -1.125 - 0.6 = -1.725 */
    assert(gp < 0.0);  /* Should be negative (evolution criterion) */
    printf("  Glansdorff-Prigogine criterion: PASSED\n");

    /* Test projected gradient */
    double X_pg[] = {10.0, 5.0};
    int constrained[] = {1, 0};
    double fixed_vals[] = {3.0, 0.0};
    int n_iter = ss_projected_gradient(X_pg, L, constrained, fixed_vals, 2,
                                        0.01, 1000, 1e-8);
    assert(n_iter >= 0);
    assert(fabs(X_pg[0] - 3.0) < TOL);
    printf("  Projected gradient: PASSED\n");

    onsager_free(L);

    printf("\n=== All Steady State tests passed ===\n");
    return 0;
}