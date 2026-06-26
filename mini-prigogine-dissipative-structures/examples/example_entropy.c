#include "pds_core.h"
#include "pds_thermo.h"
#include "pds_brusselator.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void) {
    puts("Entropy Production in Dissipative Structures");
    puts("============================================");
    PDSDissipativeSystem* sys = pds_system_create("entropy_demo", PDS_BRUSSELATOR, 2, 1, 2);
    pds_system_set_parameter(sys, 0, 1.0);
    pds_system_set_parameter(sys, 1, 3.0);
    double x0[] = {0.5, 1.0};
    double x_eq[] = {1.0, 3.0};
    pds_system_set_state(sys, x0);
    pds_system_set_equilibrium(sys, x_eq);
    puts("Far from equilibrium (b=3.0 > b_c=2.0):");
    for (int step = 0; step < 5; step++) {
        pds_entropy_balance_update(sys, 0.01, 300.0, NULL);
        double P = pds_entropy_production_rate(sys, 300.0);
        printf("  t=%.2f: P=%.6f, d_i_S=%.6f, d_e_S=%.6f, S=%.6f",
               step * 0.01, P, sys->entropy.d_i_S,
               sys->entropy.d_e_S, sys->entropy.S);
        putchar(0x0a);
    }
    putchar(0x0a);
    printf("Second Law check (d_i_S >= 0): %s",
           pds_check_second_law(sys, 1e-6) ? "PASS" : "FAIL");
    putchar(0x0a);
    PDSMinEntropyProduction* mep = pds_minimum_entropy_production_check(sys);
    printf("Minimum entropy production state: %s",
           mep->is_minimum ? "YES (near equilibrium)" : "NO (far from equilibrium)");
    putchar(0x0a);
    printf("Distance from minimum P: %.6f", mep->distance_from_minimum);
    putchar(0x0a);
    pds_min_entropy_production_free(mep);
    pds_system_free(sys);
    return 0;
}
