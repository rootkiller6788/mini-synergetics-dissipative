#include <stdio.h>
#include <math.h>

/* Example 2: BZ Reaction Oregonator Simulation
 *
 * Simulates the Belousov-Zhabotinsky (BZ) reaction using the
 * Oregonator model (Field-Koros-Noyes 1974). The BZ reaction
 * is the canonical example of a chemical oscillator — a
 * non-equilibrium phase transition to sustained temporal
 * order via a Hopf bifurcation.
 *
 * Knowledge points: L6 canonical problem (BZ oscillation),
 *   L7 application (chemical BZ oscillator).
 */

/* Oregonator declarations (matching src/reaction_diffusion.c) */
typedef struct { double k1,k2,k3,k4,k5; double A,B,f; } OregonatorParams;
typedef struct { double X,Y,Z,t; } OregonatorState;
extern double oregonator_rk4_step(OregonatorState*,const OregonatorParams*,double);
extern void oregonator_set_defaults(OregonatorParams*);
extern void oregonator_init_state(OregonatorState*);

int main(void) {
    printf("============================================================\n");
    printf("  Example 2: BZ Reaction — Oregonator Simulation\n");
    printf("============================================================\n\n");

    OregonatorParams params;
    oregonator_set_defaults(&params);
    printf("Oregonator parameters:\n");
    printf("  k1=%.2f, k2=%.1e, k3=%.1f, k4=%.1e, k5=%.1f\n",
           params.k1, params.k2, params.k3, params.k4, params.k5);
    printf("  A=[BrO3-]=%.3f, B=[organic]=%.3f, f=%.2f\n\n",
           params.A, params.B, params.f);

    OregonatorState state;
    oregonator_init_state(&state);

    double dt = 0.01;
    int steps = 5000;
    int print_every = 500;

    printf("Simulating BZ reaction for %d steps (dt=%.3f)...\n", steps, dt);
    printf("%-10s %-12s %-12s %-12s %-12s\n",
           "Time", "X=[HBrO2]", "Y=[Br-]", "Z=[Ce4+]", "|State|");
    printf("--------------------------------------------------------------\n");

    double max_amplitude = 0.0, min_amplitude = 1e15;

    for (int step = 0; step < steps; step++) {
        double amp = oregonator_rk4_step(&state, &params, dt);

        if (step > 1000) { /* Skip transient */
            if (amp > max_amplitude) max_amplitude = amp;
            if (amp < min_amplitude) min_amplitude = amp;
        }

        if (step % print_every == 0) {
            printf("%-10.2f %-12.3e %-12.3e %-12.3e %-12.3e\n",
                   state.t, state.X, state.Y, state.Z, amp);
        }
    }

    printf("\nOscillation amplitude range: [%.3e, %.3e]\n",
           min_amplitude, max_amplitude);
    printf("Oscillation detected: %s\n",
           (max_amplitude - min_amplitude > 1e-8) ? "YES (limit cycle)" : "NO");

    printf("\nKey takeaways:\n");
    printf("  1. The BZ reaction undergoes a Hopf bifurcation to oscillations\n");
    printf("  2. This is a non-equilibrium phase transition to temporal order\n");
    printf("  3. The oscillation amplitude is the order parameter\n");
    printf("  4. Relaxation oscillations are characteristic of excitability\n");
    printf("  5. The BZ reaction is a dissipative structure (Prigogine)\n");

    return 0;
}
