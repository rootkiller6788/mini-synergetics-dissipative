#include "sofe_core.h"
#include "sofe_dynamics.h"
#include "sofe_patterns.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Turing Pattern Formation: Spatial Self-Organization ===\n");
    printf("Turing instability: diffusion-driven instability\n");
    printf("Condition: activator diffuses slower than inhibitor (D_a < D_i)\n\n");

    SOFEReactionNetwork* net = sofe_network_create("turing_net");
    int u = sofe_network_add_species(net, "U", 1.0, 0.01);
    int v = sofe_network_add_species(net, "V", 1.0, 0.1);
    int r = sofe_network_add_reaction(net, "R", 1.0, 0.0);
    int reactants[] = {u}; double stoich[] = {1.0};
    sofe_network_set_reactants(net, r, reactants, stoich, 1);
    int products[] = {v}; double pstoich[] = {1.0};
    sofe_network_set_products(net, r, products, pstoich, 1);
    sofe_network_compute_jacobian(net);

    SOFEGrid* g = sofe_grid_create_1d(64, 1.0, true);
    SOFEReactionDiffusionSystem* rd = sofe_rd_create(net, g);
    double D[] = {0.01, 0.1};
    sofe_rd_set_diffusion(rd, D);

    double thresh = sofe_turing_threshold(rd, u, v);
    printf("Turing threshold parameter: %.4f\n", thresh);
    printf("D_u / D_v = %.2f\n", D[0]/D[1]);

    if (thresh > 0.0)
        printf("Turing instability IS possible with these parameters.\n");
    else
        printf("No Turing instability with these parameters.\n");

    printf("\nInitializing with noise...\n");
    sofe_field_add_noise(rd->field, 0.1, 42);

    printf("Simulating 100 steps...\n");
    for (int step = 0; step < 100; step++) {
        sofe_rd_integrate_step(rd, 0.1, SOFE_INT_EULER_EXPLICIT);
    }

    SOFEPattern* pat = sofe_pattern_detect(rd->field);
    printf("Detected pattern: %s\n", sofe_pattern_name(pat->type));
    printf("Pattern amplitude: %.4f\n", pat->amplitude);

    sofe_pattern_free(pat);
    sofe_rd_free(rd);
    sofe_grid_free(g);
    sofe_network_free(net);

    printf("\nTuring patterns are a canonical example of\n");
    printf("spatial self-organization far from equilibrium.\n");
    return 0;
}
