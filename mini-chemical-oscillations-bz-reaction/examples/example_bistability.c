/**
 * @file example_bistability.c
 * @brief End-to-end example: BZ bistability and hysteresis
 *
 * Demonstrates:
 *   - Computing steady states in different parameter regimes
 *   - Detecting bistability
 *   - Computing hysteresis loop boundaries
 *   - Chemical switch (flip-flop) demonstration
 *
 * Knowledge: L6 Canonical Problem - BZ bistability
 *
 * Compile: gcc -I../include -o example_bistable \
 *          example_bistability.c ../src/bz_reaction.c \
 *          ../src/bz_oregonator.c ../src/bz_stability.c \
 *          ../src/bz_bistability.c -lm
 */

#include <stdio.h>
#include <math.h>
#include "../include/bz_reaction.h"
#include "../include/bz_oregonator.h"
#include "../include/bz_stability.h"
#include "../include/bz_bistability.h"

int main(void)
{
    bz_oregonator_params_t p;
    bz_oregonator_state_t ss[4];
    int n_ss, i;

    printf("=== BZ Bistability and Hysteresis Analysis ===\n\n");

    /* Setup parameters in a regime that can show bistability */
    bz_oregonator_params_init(&p, 1.5);
    p.q = 0.1;  /* larger q can enable bistability */

    printf("Parameters: f=%.2f, q=%.3f, eps=%.4f\n", p.f, p.q, p.epsilon);

    /* Check for bistability */
    printf("\n--- Bistability Check ---\n");
    if (bz_is_bistable(&p)) {
        printf("System IS bistable for these parameters.\n");
    } else {
        printf("System is NOT bistable for these parameters.\n");
        printf("(Bistability requires specific f,q ranges. Trying f sweep...)\n");
    }

    /* Sweep f to find bistable region */
    printf("\n--- f-Parameter Sweep ---\n");
    {
        int results[20];
        bz_oregonator_sweep_f(0.1, 3.0, 20, results, &p);
        printf("f-sweep oscillation map:\n  f: ");
        for (i = 0; i < 20; i++) {
            double f = 0.1 + 2.9 * i / 19.0;
            printf("%.2f ", f);
        }
        printf("\n  osc: ");
        for (i = 0; i < 20; i++) {
            printf(" %d  ", results[i]);
        }
        printf("\n");
    }

    /* Compute hysteresis for given q */
    printf("\n--- Hysteresis Loop ---\n");
    {
        bz_hysteresis_t hyst;
        if (bz_compute_hysteresis(&hyst, p.q)) {
            printf("Hysteresis found!\n");
            printf("  f_low  = %.6f (lower saddle-node)\n", hyst.f_low);
            printf("  f_high = %.6f (upper saddle-node)\n", hyst.f_high);
            printf("  width  = %.6f\n", hyst.width);
            printf("  Bistable region: %.6f < f < %.6f\n",
                   hyst.f_low, hyst.f_high);
        } else {
            printf("No hysteresis found for q=%.3f\n", p.q);
        }
    }

    /* Compute and classify steady states at f=1.5 */
    printf("\n--- Steady State Analysis (f=%.2f) ---\n", p.f);
    p.f = 1.5;
    n_ss = bz_oregonator_steady_states(ss, &p);
    printf("Number of steady states: %d\n", n_ss);

    {
        bz_stability_type_t stab[4];
        bz_classify_steady_states(stab, 4, ss, n_ss, &p);

        const char *stab_names[] = {
            "Stable Node", "Stable Focus", "Unstable Node",
            "Unstable Focus", "Saddle", "Center", "Degenerate"
        };

        for (i = 0; i < n_ss; i++) {
            printf("  SS%d: (x=%.6f, y=%.6f, z=%.6f) - %s\n",
                   i, ss[i].x, ss[i].y, ss[i].z,
                   (stab[i] >= 0 && stab[i] <= 6) ? stab_names[stab[i]] : "?");
        }
    }

    /* Chemical switch demonstration */
    printf("\n--- Chemical Switch Demo ---\n");
    if (n_ss >= 2) {
        bz_oregonator_state_t state;
        int result;

        /* Start near first steady state */
        state = ss[0];
        printf("Initial state: (x=%.6f, z=%.6f)\n", state.x, state.z);

        /* Try to switch to upper state */
        result = bz_chemical_switch(&state, &p, 0.5, 1);
        printf("Switch attempt (target=upper): %s\n",
               result ? "SUCCESS" : "FAILED");
        printf("State after switch: (x=%.6f, z=%.6f)\n", state.x, state.z);

        /* Try to switch back */
        result = bz_chemical_switch(&state, &p, 0.5, 0);
        printf("Switch attempt (target=lower): %s\n",
               result ? "SUCCESS" : "FAILED");
        printf("State after switch: (x=%.6f, z=%.6f)\n", state.x, state.z);
    } else {
        printf("Need 2+ steady states for switch demo.\n");
    }

    printf("\nDone.\n");
    return 0;
}
