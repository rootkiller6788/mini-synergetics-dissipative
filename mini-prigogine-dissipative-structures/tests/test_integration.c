#include "pds_core.h"
#include "pds_thermo.h"
#include "pds_brusselator.h"
#include "pds_stability.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

int main(void) {
    PDSDissipativeSystem* sys = pds_system_create("int", PDS_BRUSSELATOR, 2, 1, 2);
    pds_system_set_parameter(sys, 0, 1.0);
    pds_system_set_parameter(sys, 1, 3.0);
    double x0[] = {0.5, 1.0};
    pds_system_set_state(sys, x0);
    double P = pds_entropy_production_rate(sys, 300.0);
    assert(P >= -1e-10);
    assert(pds_check_second_law(sys, 1e-6));
    PDSGlansdorffPrigogineStability* gps = pds_glansdorff_prigogine_stability(sys);
    assert(gps != NULL);
    PDSMatrix* Jmat = pds_matrix_create(2, 2);
    pds_matrix_set(Jmat, 0, 0, 2.0); pds_matrix_set(Jmat, 0, 1, 1.0);
    pds_matrix_set(Jmat, 1, 0, -3.0); pds_matrix_set(Jmat, 1, 1, -1.0);
    PDSStabilityResult* sr = pds_linear_stability_analysis(Jmat);
    assert(sr != NULL);
    pds_stability_result_free(sr);
    pds_matrix_free(Jmat);
    pds_gp_stability_free(gps);
    pds_system_free(sys);
    puts("ALL INTEGRATION TESTS PASSED");
    return 0;
}
