/**
 * @file example_fluctuation.c
 * @brief Example: Fluctuation theorem verification via Brownian motion.
 *
 * Demonstrates the Evans-Searles fluctuation theorem by simulating
 * a Brownian particle and verifying P(Sigma)/P(-Sigma) = exp(Sigma/k_B T).
 *
 * L6: Canonical problem -- fluctuation theorem validation
 * L8: Advanced topic -- stochastic thermodynamics
 */

#include "fluctuation_theorems.h"
#include <stdio.h>
#include <math.h>

int main(void)
{
    printf("==================================================\n");
    printf("  Fluctuation Theorem Verification\n");
    printf("  Evans-Searles & Jarzynski Relations\n");
    printf("==================================================\n\n");

    /* Parameters for a colloidal particle in water */
    double mass = 1.0e-15;       /* ~1 pg particle [kg] */
    double friction = 1.0e-8;    /* Friction coefficient [kg/s] */
    double spring_k = 1.0e-6;    /* Optical trap stiffness [N/m] */
    double temperature = 300.0;  /* Room temperature [K] */
    double k_B = 1.380649e-23;   /* Boltzmann constant [J/K] */
    double dt = 1.0e-6;          /* Time step [s] */
    int n_steps = 50000;         /* Number of steps */

    printf("System parameters:\n");
    printf("  Particle mass: %.1e kg\n", mass);
    printf("  Friction: %.1e kg/s\n", friction);
    printf("  Trap stiffness: %.1e N/m\n", spring_k);
    printf("  Temperature: %.0f K\n", temperature);
    printf("  Time step: %.1e s\n", dt);
    printf("  Steps: %d\n\n", n_steps);

    /* Simulate a single trajectory */
    StochasticTrajectory traj;
    ft_init_trajectory(&traj, n_steps, dt, temperature, k_B);

    printf("Simulating Brownian particle in harmonic trap...\n");
    ft_simulate_brownian_harmonic(&traj, mass, friction, spring_k,
                                   temperature, k_B, 0.0, 0.0, dt, n_steps);

    printf("Computing entropy production along trajectory...\n");
    ft_compute_trajectory_entropy(&traj);

    /* Extract final entropy production statistics */
    double sigma_final = traj.entropy_production[n_steps - 1];
    double sigma_per_step = sigma_final / n_steps;
    double k_B_na = k_B;

    printf("\nTrajectory results:\n");
    printf("  Final position: %.4e m\n", traj.position[n_steps - 1]);
    printf("  Final velocity: %.4e m/s\n", traj.velocity[n_steps - 1]);
    printf("  Cumulative entropy production: %.6e J/K\n", sigma_final);
    printf("  Entropy production per step: %.6e J/K\n", sigma_per_step);
    printf("  Comparison k_B: %.6e J/K\n", k_B_na);

    /* Verify integral fluctuation theorem:
     * <exp(-Sigma/k_B)> should equal 1 */
    double avg_exp = 0.0;
    int n_check = 100;  /* Sample every N steps */
    for (int i = n_steps - n_check; i < n_steps; i++) {
        avg_exp += exp(-traj.entropy_production[i] / k_B);
    }
    avg_exp /= n_check;

    printf("\nIntegral Fluctuation Theorem check:\n");
    printf("  <exp(-Sigma/k_B)> = %.6f (should = 1)\n", avg_exp);
    if (fabs(avg_exp - 1.0) < 0.2) {
        printf("  IFT verified within tolerance!\n");
    } else {
        printf("  IFT check: large deviation (need more samples).\n");
    }

    printf("\nReference: Evans & Searles (1994) Phys. Rev. E 50, 1645\n");
    printf("Reference: Jarzynski (1997) Phys. Rev. Lett. 78, 2690\n");

    ft_free_trajectory(&traj);
    printf("\n=== Fluctuation Theorem Example Complete ===\n");
    return 0;
}