#include "pds_core.h"
#include "pds_thermo.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

int main(void) {
    PDSDissipativeSystem* sys = pds_system_create("t", PDS_BRUSSELATOR, 2, 1, 2);
    double x0[] = {1.0, 2.0}, x_eq[] = {0.5, 1.0};
    pds_system_set_state(sys, x0);
    pds_system_set_equilibrium(sys, x_eq);
    pds_entropy_balance_update(sys, 0.01, 300.0, NULL);
    assert(sys->entropy.d_i_S >= -1e-10);
    assert(pds_check_second_law(sys, 1e-8));
    PDSOnsagerMatrix* om = pds_onsager_matrix_create(2);
    pds_onsager_matrix_set(om, 0, 0, 2.0); pds_onsager_matrix_set(om, 0, 1, 0.5);
    pds_onsager_matrix_set(om, 1, 0, 0.5); pds_onsager_matrix_set(om, 1, 1, 2.0);
    pds_onsager_matrix_check_symmetry(om);
    assert(om->is_symmetric);
    assert(pds_onsager_matrix_is_positive_definite(om));
    double forces[] = {1.0, 2.0}, fluxes[2];
    pds_onsager_matrix_compute_fluxes(om, forces, fluxes);
    PDSGlansdorffPrigogineStability* gps = pds_glansdorff_prigogine_stability(sys);
    assert(gps != NULL);
    PDSExtendedThermo* eit = pds_extended_thermo_create(3);
    pds_extended_thermo_evolve(eit, 0.01);
    pds_gp_stability_free(gps);
    pds_onsager_matrix_free(om);
    pds_extended_thermo_free(eit);
    pds_system_free(sys);
    puts("ALL THERMO TESTS PASSED");
    return 0;
}
