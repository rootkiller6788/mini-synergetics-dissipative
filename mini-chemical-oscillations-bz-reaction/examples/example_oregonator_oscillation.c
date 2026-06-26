/**
 * @file example_oregonator_oscillation.c
 * @brief End-to-end example: Oregonator oscillation simulation
 *
 * Demonstrates:
 *   - Setting up Oregonator parameters
 *   - Integrating the ODE system with RK4
 *   - Detecting oscillations and measuring period/amplitude
 *
 * Knowledge: L6 Canonical Problem - BZ oscillation
 *
 * Compile: gcc -I../include -o example_oregonator \
 *          example_oregonator_oscillation.c ../src/bz_reaction.c \
 *          ../src/bz_oregonator.c -lm
 */

#include <stdio.h>
#include <math.h>
#include "../include/bz_reaction.h"
#include "../include/bz_oregonator.h"

int main(void)
{
    bz_oregonator_params_t p;
    bz_oregonator_state_t s;
    int n_peaks;
    double period, amplitude;
    double T_analytic;

    printf("=== BZ Oregonator Oscillation Simulation ===\n\n");

    /* Setup oscillatory regime parameters */
    bz_oregonator_params_init(&p, 0.8);
    printf("Parameters: f=%.2f, q=%.2e, eps=%.2e, eps'=%.2e\n",
           p.f, p.q, p.epsilon, p.epsilon_p);

    /* Initialize near the unstable steady state */
    bz_oregonator_init_oscillatory(&s, &p);
    printf("Initial state: x=%.6f, y=%.6f, z=%.6f\n", s.x, s.y, s.z);

    /* Analytical period estimate */
    T_analytic = bz_oregonator_period_approx(&p);
    printf("Analytical period estimate: T ~ %.3f (dimensionless)\n", T_analytic);

    /* Check oscillatory regime */
    printf("Oscillatory regime: %s\n",
           bz_oregonator_is_oscillatory(&p) ? "YES" : "NO");

    /* Integrate and analyze */
    printf("\nIntegrating for dimensionless time T=50 with 5000 steps...\n");
    bz_oregonator_integrate_and_analyze(&s, &p, 50.0, 5000,
                                         &n_peaks, &period, &amplitude);

    printf("\n=== Results ===\n");
    printf("Number of oscillation peaks: %d\n", n_peaks);
    printf("Average period: %.6f (dimensionless)\n", period);
    printf("Amplitude (x peak-to-peak): %.6f\n", amplitude);
    printf("Final state: x=%.6f, y=%.6f, z=%.6f\n", s.x, s.y, s.z);

    /* Compare with analytical estimate */
    if (n_peaks >= 2) {
        printf("\nPeriod comparison:\n");
        printf("  Numerical:  T_num = %.6f\n", period);
        printf("  Analytical: T_ana = %.6f\n", T_analytic);
        if (T_analytic > 0) {
            printf("  Relative error: %.2f%%\n",
                   100.0 * fabs(period - T_analytic) / T_analytic);
        }
    }

    /* Classification */
    {
        /* Quick classification based on relaxation ratio */
        double rel_ratio = bz_oregonator_relaxation_ratio(&p);
        printf("\nOscillation type: ");
        if (rel_ratio > 100) {
            printf("RELAXATION (sharp spikes, slow recovery)\n");
        } else if (rel_ratio > 10) {
            printf("MIXED (some relaxation character)\n");
        } else {
            printf("SINUSOIDAL (near Hopf onset)\n");
        }
        printf("Relaxation ratio (eps/eps'): %.1f\n", rel_ratio);
    }

    printf("\nDone.\n");
    return 0;
}
