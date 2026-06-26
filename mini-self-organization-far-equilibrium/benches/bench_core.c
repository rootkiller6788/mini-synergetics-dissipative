#include "sofe_core.h"
#include "sofe_thermodynamics.h"
#include "sofe_dynamics.h"
#include <stdio.h>
#include <time.h>
#include <math.h>

static double get_time_sec(void) {
    return (double)clock() / CLOCKS_PER_SEC;
}

int main(void) {
    printf("=== SOFE Benchmarks ===

");
    double t0, t1;

    /* Benchmark: Jacobian computation */
    printf("Benchmark 1: Jacobian computation (1000 runs)...
");
    SOFEReactionNetwork* net = sofe_network_create("bench");
    for (int i = 0; i < 10; i++) {
        char name[8]; sprintf(name, "S%d", i);
        sofe_network_add_species(net, name, 1.0, 0.01);
    }
    for (int r = 0; r < 5; r++) {
        char rname[8]; sprintf(rname, "R%d", r);
        int ridx = sofe_network_add_reaction(net, rname, 1.0, 0.0);
        int reactants[] = {r*2}; double stoich[] = {1.0};
        sofe_network_set_reactants(net, ridx, reactants, stoich, 1);
        int products[] = {r*2+1}; double pstoich[] = {1.0};
        sofe_network_set_products(net, ridx, products, pstoich, 1);
    }
    t0 = get_time_sec();
    for (int i = 0; i < 1000; i++) {
        sofe_network_compute_jacobian(net);
    }
    t1 = get_time_sec();
    printf("  Jacobian x1000: %.3f ms
", (t1-t0)*1000.0);
    sofe_network_free(net);

    /* Benchmark: Reaction-diffusion integration */
    printf("Benchmark 2: RD integration (100 steps on 64x64 grid)...
");
    SOFEReactionNetwork* net2 = sofe_network_create("bench2");
    sofe_network_add_species(net2, "U", 1.0, 0.01);
    sofe_network_add_species(net2, "V", 1.0, 0.05);
    int r2 = sofe_network_add_reaction(net2, "R", 1.0, 0.0);
    int reactants[] = {0}; double stoich[] = {1.0};
    sofe_network_set_reactants(net2, r2, reactants, stoich, 1);
    int products[] = {1}; double pstoich[] = {1.0};
    sofe_network_set_products(net2, r2, products, pstoich, 1);
    SOFEGrid* g = sofe_grid_create_2d(64, 64, 1.0, 1.0, true, true);
    SOFEReactionDiffusionSystem* rd = sofe_rd_create(net2, g);
    double D[] = {0.01, 0.05};
    sofe_rd_set_diffusion(rd, D);
    t0 = get_time_sec();
    for (int i = 0; i < 100; i++) {
        sofe_rd_integrate_step(rd, 0.01, SOFE_INT_EULER_EXPLICIT);
    }
    t1 = get_time_sec();
    printf("  RD 100 steps: %.3f ms
", (t1-t0)*1000.0);
    sofe_rd_free(rd);
    sofe_grid_free(g);
    sofe_network_free(net2);

    /* Benchmark: Entropy production */
    printf("Benchmark 3: Entropy production (10000 runs)...
");
    SOFEForceFluxSystem ffs;
    ffs.n_pairs = 5;
    ffs.pairs = (SOFEForceFluxPair*)malloc(5 * sizeof(SOFEForceFluxPair));
    for (int i = 0; i < 5; i++) {
        ffs.pairs[i].force = 1.0 + i * 0.5;
        ffs.pairs[i].flux = 0.2 + i * 0.1;
    }
    t0 = get_time_sec();
    for (int i = 0; i < 10000; i++) {
        sofe_entropy_production_direct(&ffs);
    }
    t1 = get_time_sec();
    printf("  Entropy x10000: %.3f ms
", (t1-t0)*1000.0);
    free(ffs.pairs);

    printf("
Benchmarks complete.
");
    return 0;
}
