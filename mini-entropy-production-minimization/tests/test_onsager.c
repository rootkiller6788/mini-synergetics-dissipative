/**
 * @file test_onsager.c
 * @brief Additional tests for Onsager matrix operations and reciprocity.
 */

#include "onsager_relations.h"
#include "entropy_production.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define TOL 1e-10

int main(void)
{
    printf("=== Onsager Tests ===\n\n");

    /* Test Fourier initialization */
    OnsagerMatrix *L = onsager_alloc(1);
    assert(L != NULL);
    onsager_init_fourier(L, 0.6, 300.0);
    /* L = kappa * T^2 = 0.6 * 90000 = 54000 */
    assert(fabs(onsager_get(L, 0, 0) - 54000.0) < 1.0);
    printf("  Fourier init: PASSED\n");
    onsager_free(L);

    /* Test reciprocity verification with asymmetric matrix */
    L = onsager_alloc(2);
    onsager_set(L, 0, 0, 1.0);
    onsager_set(L, 0, 1, 2.0);
    onsager_set(L, 1, 0, 3.0);  /* Asymmetric! */
    onsager_set(L, 1, 1, 4.0);
    assert(onsager_verify_reciprocity(L, 1e-10) == 0);
    printf("  Asymmetry detection: PASSED\n");
    onsager_free(L);

    /* Test magnetic reversal */
    L = onsager_alloc(2);
    onsager_init_two_process(L, 4.0, 1.0, 1.0);
    int parity[] = {1, -1};  /* Force 0 even, force 1 odd under T-reversal */
    onsager_magnetic_reversal(L, parity, 2);
    /* L_{01} -> epsilon_0 * epsilon_1 * L_{01} = 1 * (-1) * 1.0 = -1.0 */
    assert(fabs(onsager_get(L, 0, 1) - (-1.0)) < TOL);
    printf("  Magnetic reversal: PASSED\n");
    onsager_free(L);

    /* Test condition number */
    L = onsager_alloc(2);
    onsager_init_diagonal(L, (double[]){100.0, 1.0});
    double cond = onsager_condition_number(L);
    assert(fabs(cond - 100.0) < 1e-6);
    printf("  Condition number: PASSED\n");
    onsager_free(L);

    /* Test force computation from fluxes */
    L = onsager_alloc(2);
    onsager_set(L, 0, 0, 2.0);
    onsager_set(L, 1, 1, 2.0);
    double J[] = {4.0, 6.0};
    double X[2];
    onsager_compute_forces(L, J, X);
    /* L = diag(2,2), R = diag(0.5, 0.5), X = R*J = [2, 3] */
    assert(fabs(X[0] - 2.0) < TOL);
    assert(fabs(X[1] - 3.0) < TOL);
    printf("  Force from fluxes: PASSED\n");
    onsager_free(L);

    /* Cross entropy production */
    double Lmat[] = {1, 0.5, 0.5, 1};
    double X1[] = {2.0, 0.0};
    double X2[] = {0.0, 3.0};
    double cross = ep_cross_entropy_production(X1, X2, Lmat, 2);
    /* X1^T L X2 = 2*1*0 + 2*0.5*3 + 0*0.5*0 + 0*1*3 = 3.0 */
    assert(fabs(cross - 3.0) < TOL);
    printf("  Cross entropy production: PASSED\n");

    printf("\n=== All Onsager tests passed ===\n");
    return 0;
}