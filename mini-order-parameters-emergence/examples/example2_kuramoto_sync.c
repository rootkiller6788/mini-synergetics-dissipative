#include "oep_core.h"
#include "oep_collective_dynamics.h"
#include <stdio.h>
#include <math.h>

/* ==============================================================
 * Example 2: Kuramoto Synchronization Transition
 *
 * The Kuramoto model is a paradigmatic example of order parameter
 * emergence in a system of coupled oscillators.
 *
 * For coupling K < K_c: oscillators are incoherent (r ≈ 0)
 * For coupling K > K_c: spontaneous synchronization emerges (r > 0)
 *
 * This is a continuous (second-order) nonequilibrium phase transition.
 * The order parameter r measures the degree of phase coherence.
 * ============================================================== */

int main(void) {
    printf("=== Example 2: Kuramoto Synchronization ===\n\n");

    int N = 100;  /* Number of oscillators */
    int n_K = 21;
    double K_min = 0.0, K_max = 4.0;
    double freq_spread = 1.0;

    double* K_values = malloc((size_t)n_K * sizeof(double));
    double* r_values = malloc((size_t)n_K * sizeof(double));

    printf("Computing synchronization transition r(K)...\n");
    printf("N = %d oscillators, frequency spread σ = %.1f\n", N, freq_spread);
    printf("Critical coupling K_c (theory) = √(8/π)·σ ≈ %.3f\n",
           oep_kura_critical_coupling(freq_spread));

    oep_kura_sync_transition(K_min, K_max, n_K, freq_spread, N, K_values, r_values);

    printf("\nSynchronization curve:\n");
    printf("  K       r(K)      State\n");
    printf("  ------- --------- ----------\n");
    for (int i = 0; i < n_K; i++) {
        const char* state = (r_values[i] < 0.1) ? "INCOHERENT" :
                            (r_values[i] < 0.5) ? "TRANSITION" : "SYNCHRONIZED";
        printf("  %7.3f  %9.4f  %s\n", K_values[i], r_values[i], state);
    }

    /* Show phase distribution at different couplings */
    printf("\nPhase snapshots:\n");
    for (int k_idx = 0; k_idx < 3; k_idx++) {
        double K_test[] = {0.5, 2.0, 4.0};
        double K = K_test[k_idx];

        OEPKuramotoSystem* sys = oep_kura_create(N);
        oep_kura_set_frequencies_distribution(sys, 0.0, freq_spread);
        oep_kura_set_coupling(sys, K);

        /* Equilibrate */
        for (int t = 0; t < 500; t++) oep_kura_step(sys, 0.01);
        oep_kura_compute_order_parameter(sys);

        printf("\n  K = %.1f: r = %.4f", K, sys->mean_field_coherence);
        printf("  ψ = %.4f rad", sys->mean_field_phase);
        printf("  (first 10 phases: ");
        for (int i = 0; i < 10; i++)
            printf("%.2f ", sys->oscillators[i].phase);
        printf(")\n");

        oep_kura_free(sys);
    }

    /* Time evolution of order parameter */
    printf("\nOrder parameter time evolution (K = 2.5):\n");
    OEPKuramotoSystem* sys = oep_kura_create(N);
    oep_kura_set_frequencies_distribution(sys, 0.0, freq_spread);
    oep_kura_set_coupling(sys, 2.5);

    double* r_traj = malloc((size_t)100 * sizeof(double));
    for (int t = 0; t < 100; t++) {
        oep_kura_step(sys, 0.01);
        oep_kura_compute_order_parameter(sys);
        r_traj[t] = sys->mean_field_coherence;
        if (t % 10 == 0)
            printf("    t=%3d: r=%.4f\n", t, r_traj[t]);
    }

    free(r_traj);
    oep_kura_free(sys);
    free(K_values);
    free(r_values);

    printf("\n=== Example 2 Complete ===\n");
    return 0;
}
