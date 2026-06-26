#include "sofe_core.h"
#include "sofe_thermodynamics.h"
#include "sofe_dynamics.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("============================================================\n");
    printf("  Self-Organization Far-from-Equilibrium -- Demo Overview\n");
    printf("============================================================\n\n");

    printf("This module implements the theory of self-organization in\n");
    printf("systems driven far from thermodynamic equilibrium.\n\n");

    SOFESystem* sys = sofe_system_create("demo_system");
    sofe_system_set_regime(sys, SOFE_REGIME_FAR_EQUILIBRIUM);
    sofe_system_set_temperature(sys, 300.0);

    SOFEReactionNetwork* net = sofe_network_create("brusselator");
    int x = sofe_network_add_species(net, "X", 1.0, 0.01);
    int y = sofe_network_add_species(net, "Y", 1.0, 0.02);
    sofe_species_set_autocatalytic(&net->species[x], true);

    int r = sofe_network_add_reaction(net, "autocat", 1.0, 0.0);
    int reactants[] = {x, x, y}; double rstoich[] = {2.0, 1.0};
    sofe_network_set_reactants(net, r, reactants, rstoich, 3);
    int products[] = {x}; double pstoich[] = {3.0};
    sofe_network_set_products(net, r, products, pstoich, 1);
    sofe_network_set_autocatalytic(net, r, true);
    sofe_network_set_kinetics(net, r, 0, 3, 1.0, 0.0);

    sys->network = net;
    sofe_network_compute_jacobian(net);
    sofe_network_compute_fluxes(net);

    sys->entropy_production_rate = sofe_entropy_production_reactions(net, sys->temperature);
    sys->has_macroscopic_order = (sys->entropy_production_rate > 0.1);

    sofe_system_add_order_parameter(sys, "osc_amp", SOFE_ORDER_SCALAR,
                                     SOFE_SYMMETRY_Z2, 1);
    sofe_system_set_order_amplitude(sys, 0, 0, 1.5);

    sofe_system_print_summary(sys);

    printf("\nBrusselator simulation (10 steps):\n");
    double state[] = {1.0, 1.0};
    double params[] = {1.0, 3.0};
    printf("t=0.0: X=%.4f Y=%.4f\n", state[0], state[1]);
    for (int i = 0; i < 10; i++) {
        sofe_integrate_step(state, 2, i*0.1, 0.1,
                            sofe_brusselator_rates, params, 2, SOFE_INT_RK4);
    }
    printf("t=1.0: X=%.4f Y=%.4f\n", state[0], state[1]);

    sofe_system_free(sys);
    printf("\nDemo complete.\n");
    return 0;
}
