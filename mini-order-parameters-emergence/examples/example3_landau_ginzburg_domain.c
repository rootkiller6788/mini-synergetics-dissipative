#include "oep_core.h"
#include "oep_landau_ginzburg.h"
#include <stdio.h>
#include <math.h>

/* ==============================================================
 * Example 3: Landau-Ginzburg Domain Wall Dynamics
 *
 * Demonstrates the spatiotemporal evolution of an order parameter
 * field governed by the time-dependent Ginzburg-Landau equation:
 *
 *   ∂ξ/∂t = -Γ [aξ + bξ³ - c∇²ξ]
 *
 * Starting from a random initial condition, the system evolves to
 * form domains of positive and negative order parameter, separated
 * by domain walls. This illustrates:
 *
 *   1. Phase ordering kinetics (coarsening)
 *   2. Domain wall structure (kink solution: ξ ∝ tanh)
 *   3. Free energy minimization
 *   4. Correlation length and pattern formation
 *
 * Reference: Bray (1994) "Theory of phase-ordering kinetics"
 * ============================================================== */

int main(void) {
    printf("=== Example 3: Landau-Ginzburg Domain Walls ===\n\n");

    /* System parameters */
    int nx = 64, ny = 64;
    double Lx = 10.0, Ly = 10.0;
    double a = -1.0;   /* Below critical point → ordered phase */
    double b = 1.0;
    double c = 0.5;    /* Gradient stiffness */
    double Gamma = 1.0;
    double dt = 0.001;

    /* Create the Landau-Ginzburg system */
    OEPLandauGinzburg* lg = oep_lg_create(nx, ny, Lx, Ly, a, b, c, Gamma);
    OEPSpatialField* field = oep_field_create(nx, ny, Lx, Ly);
    if (!lg || !field) {
        printf("ERROR: Failed to create system\n");
        return 1;
    }

    /* Initialize with small random fluctuations around zero */
    oep_field_set_random(field, 0.1);

    printf("Initial state:\n");
    printf("  Spatial mean:     %.6f\n", oep_field_spatial_mean(field));
    printf("  Spatial variance: %.6f\n", oep_field_spatial_variance(field));

    /* Compute theoretical values */
    double xi_eq = sqrt(fabs(a) / b);
    double xi_corr = oep_lg_correlation_length(a, c);
    double wall_width = oep_lg_domain_wall_width(a, c);
    double wall_energy = oep_lg_domain_wall_energy(a, b, c);

    printf("\nTheoretical predictions:\n");
    printf("  Equilibrium amplitude ξ₀:  %.4f\n", xi_eq);
    printf("  Correlation length ξ:      %.4f\n", xi_corr);
    printf("  Domain wall width δ:       %.4f\n", wall_width);
    printf("  Domain wall energy σ:      %.4f\n", wall_energy);

    /* Evolve the system */
    int n_steps = 1000;
    printf("\nEvolving for %d steps (dt = %.4f)...\n", n_steps, dt);

    double* mean_ts = malloc((size_t)n_steps * sizeof(double));
    double* var_ts = malloc((size_t)n_steps * sizeof(double));

    /* Create 1D domain wall for comparison */
    OEPSpatialField* wall = oep_field_create(nx, 1, Lx, 1.0);
    oep_lg_domain_wall_1d(wall, Lx * 0.5, a, b, c, nx);

    printf("\nDomain wall profile (1D, analytic):\n");
    printf("  x        ξ(x)\n");
    for (int i = 0; i < nx; i += 4) {
        printf("  %5.2f   %8.4f\n",
               (double)i * wall->dx, oep_field_get(wall, i, 0));
    }

    /* Evolve the 2D field */
    oep_lg_evolve(lg, field, dt, n_steps, mean_ts, var_ts);

    printf("\nEvolution summary (every 200 steps):\n");
    printf("  Step    <ξ>          Var(ξ)       Free Energy\n");
    printf("  ------  ------------ ------------ ------------\n");
    for (int step = 0; step < n_steps; step += 200) {
        double F = oep_field_free_energy(field, a, b, c, field->dx, field->dy);
        printf("  %6d  %12.6f %12.6f %12.6f\n",
               step, mean_ts[step], var_ts[step], F);
    }

    printf("\nFinal state:\n");
    printf("  Spatial mean:     %.6f (theory: ±%.4f)\n",
           oep_field_spatial_mean(field), xi_eq);
    printf("  Spatial variance: %.6f\n", oep_field_spatial_variance(field));
    printf("  Free energy:      %.4f\n",
           oep_field_free_energy(field, a, b, c, field->dx, field->dy));

    /* Superconductivity parameters */
    printf("\nSuperconductivity application:\n");
    double Tc = 9.2;    /* Niobium T_c */
    double T = 4.2;     /* Liquid helium */
    double xi0 = 38.0;  /* nm */
    double lambda0 = 39.0;  /* nm */

    double gap = oep_lg_superconductor_gap(T, Tc, 1.0, 1.0);
    double xi_T = oep_lg_coherence_length(T, Tc, xi0);
    double lambda_T = oep_lg_penetration_depth(T, Tc, lambda0);
    double kappa = oep_lg_ginzburg_parameter(lambda_T, xi_T);

    printf("  T_c = %.1f K, T = %.1f K\n", Tc, T);
    printf("  Gap |Ψ| = %.4f\n", gap);
    printf("  Coherence length ξ(T) = %.2f nm\n", xi_T);
    printf("  Penetration depth λ(T) = %.2f nm\n", lambda_T);
    printf("  Ginzburg parameter κ = %.4f\n", kappa);
    printf("  Type: %s (κ %s 1/√2)\n",
           oep_lg_is_type_II(kappa) ? "TYPE II" : "TYPE I",
           oep_lg_is_type_II(kappa) ? ">" : "<");

    /* Cleanup */
    free(mean_ts);
    free(var_ts);
    oep_field_free(field);
    oep_field_free(wall);
    oep_lg_free(lg);

    printf("\n=== Example 3 Complete ===\n");
    return 0;
}
