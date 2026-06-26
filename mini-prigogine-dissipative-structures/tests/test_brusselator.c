#include "pds_brusselator.h"
#include "pds_stability.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

int main(void) {
    PDSBrusselatorParams params = {1.0, 3.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    double xs, ys;
    pds_brusselator_steady_state(&params, &xs, &ys);
    assert(fabs(xs - 1.0) < 1e-10);
    assert(fabs(ys - 3.0) < 1e-10);

    double bc = pds_brusselator_bifurcation_point(&params);
    assert(fabs(bc - 2.0) < 1e-10);

    PDSBrusselatorState state = {1.0, 3.0};
    double J[2][2];
    pds_brusselator_jacobian(&state, &params, J);
    assert(fabs(J[0][0] - 2.0) < 1e-10);
    assert(fabs(J[0][1] - 1.0) < 1e-10);
    assert(fabs(J[1][0] - (-3.0)) < 1e-10);
    assert(fabs(J[1][1] - (-1.0)) < 1e-10);

    double l1r, l1i, l2r, l2i;
    pds_eigenvalues_2x2(J[0][0], J[0][1], J[1][0], J[1][1],
                        &l1r, &l1i, &l2r, &l2i);
    assert(fabs(l1r - 0.5) < 1e-8);
    assert(fabs(l2r - 0.5) < 1e-8);

    PDSBrusselatorState initial = {0.5, 1.0};
    PDSBrusselatorParams pstable = {1.0, 1.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    PDSBrusselatorResult* result = pds_brusselator_simulate(
        &initial, &pstable, 5.0, 0.01, 10);
    assert(result != NULL);
    assert(result->n_points > 0);
    pds_brusselator_result_free(result);

    PDSBrusselatorResult* result2 = pds_brusselator_simulate(
        &initial, &params, 20.0, 0.01, 20);
    assert(result2 != NULL);
    pds_brusselator_result_free(result2);

    puts("ALL BRUSSELATOR TESTS PASSED");
    return 0;
}
