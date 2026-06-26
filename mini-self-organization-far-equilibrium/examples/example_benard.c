#include "sofe_core.h"
#include "sofe_thermodynamics.h"
#include "sofe_dynamics.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Benard Convection: Rayleigh-Benard Instability ===\n");
    printf("A fluid layer heated from below spontaneously forms\n");
    printf("convection rolls when the temperature gradient exceeds\n");
    printf("a critical threshold (Rayleigh number > Ra_c).\n\n");

    double Ra_critical = 1707.8;
    double Ra_values[] = {500.0, 1500.0, 1707.8, 2000.0, 5000.0};

    printf("Rayleigh Number | Regime | Entropy Production\n");
    printf("----------------|--------|-------------------\n");

    for (int i = 0; i < 5; i++) {
        double Ra = Ra_values[i];
        const char* regime;
        double sigma_est = 0.0;

        if (Ra < Ra_critical) {
            regime = "Conduction (near-eq)";
            sigma_est = Ra / Ra_critical * 0.1;
        } else if (Ra < 10.0 * Ra_critical) {
            regime = "Convection rolls";
            sigma_est = 0.5 + (Ra - Ra_critical) / Ra_critical * 0.1;
        } else {
            regime = "Turbulent convection";
            sigma_est = 1.0 + log(Ra / Ra_critical);
        }

        printf("%14.1f | %-18s | %.4f\n", Ra, regime, sigma_est);
    }

    printf("\nThe Benard instability demonstrates:\n");
    printf("1. Threshold behavior (bifurcation at Ra_c)\n");
    printf("2. Symmetry breaking (roll direction)\n");
    printf("3. Increased entropy production in ordered state\n");
    printf("4. Hierarchical self-organization (rolls -> turbulence)\n");
    return 0;
}
