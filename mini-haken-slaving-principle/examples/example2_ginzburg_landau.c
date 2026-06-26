/* ==============================================================
 * example2_ginzburg_landau.c — Ginzburg-Landau Equation Analysis
 *
 * Demonstrates the Ginzburg-Landau theory of second-order phase
 * transitions in the context of the slaving principle.
 *
 * The GL free energy:  F(ξ) = (a/2)ξ² + (b/4)ξ⁴  (b > 0)
 * Equation of motion:  dξ/dt = -∂F/∂ξ = -a·ξ - b·ξ³
 *
 * Below threshold (a > 0): single stable fixed point at ξ = 0
 * Above threshold (a < 0): three fixed points, ξ = 0 unstable,
 *                           ξ = ±√(-a/b) stable (symmetry broken)
 *
 * This example shows:
 *   1. Free energy landscape and potential wells
 *   2. Fixed point structure vs control parameter
 *   3. Critical exponents and scaling
 *   4. Correlation length divergence
 *
 * Reference: Landau & Lifshitz (1980) Statistical Physics, Part 1, §146
 * ============================================================== */

#include "haken_order_param.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

int main(void) {
    printf("=== Ginzburg-Landau Phase Transition Analysis ===\n\n");

    /* Parameter sweep: a from +1 (below threshold) to -1 (above) */
    printf("GL equation: dξ/dt = -a·ξ - b·ξ³\n");
    printf("Parameter sweep: a ∈ [-1, 1], b = 1\n\n");

    printf("%10s %10s %14s %10s %10s\n",
           "a", "n_fixed", "ξ_eq (stable)", "suscept", "corr_len");

    for (int i = 0; i <= 20; i++) {
        double a = 1.0 - 0.1 * i;  /* from 1.0 to -1.0 */
        Haken_GLParameters params = {a, 1.0, 0.0, 1.0, 1.0, 0.5, 2.0, 0.5};
        double fixed[3];
        int n_fixed = haken_gl_fixed_points(&params, fixed);

        char eq_str[32] = "none";
        double stable_eq = 0.0;
        if (n_fixed == 1) {
            snprintf(eq_str, sizeof(eq_str), "0");
            stable_eq = 0.0;
        } else if (n_fixed == 3) {
            snprintf(eq_str, sizeof(eq_str), "±%.4f", fixed[1]);
            stable_eq = fixed[1];
        }

        /* Susceptibility χ = 1/|a| (mean-field) */
        double chi = (fabs(a) > 1e-10) ? (1.0 / fabs(a)) : INFINITY;
        /* Correlation length ξ_corr ~ 1/√|a| */
        double xi = (fabs(a) > 1e-10) ? (1.0 / sqrt(fabs(a))) : INFINITY;

        printf("%10.4f %10d %14s %10.4f %10.4f\n",
               a, n_fixed, eq_str,
               (chi < 100) ? chi : INFINITY,
               (xi < 100) ? xi : INFINITY);
    }

    /* Free energy landscape */
    printf("\n\nFree energy landscape F(ξ) = (a/2)ξ² + (b/4)ξ⁴:\n\n");

    double a_vals[3] = {0.5, 0.0, -0.5};
    for (int ai = 0; ai < 3; ai++) {
        double a = a_vals[ai];
        printf("a = %+.1f:\n", a);
        printf("  %8s %12s\n", "ξ", "F(ξ)");
        for (int j = 0; j <= 10; j++) {
            double xi = -1.5 + 0.3 * j;
            Haken_GLParameters params = {a, 1.0, 0.0, 1.0, 1.0, 0.5, 2.0, 0.5};
            double F = haken_gl_free_energy(xi, &params);
            printf("  %8.3f %12.6f\n", xi, F);
        }
        printf("\n");
    }

    /* Critical exponents */
    printf("Critical exponents (mean-field theory):\n");
    printf("  β = 1/2  → ξ_eq ~ |a|^(1/2) above threshold\n");
    printf("  γ = 1    → χ ~ |a|^(-1)\n");
    printf("  ν = 1/2  → ξ_corr ~ |a|^(-1/2)\n");

    printf("\n=== Done ===\n");
    return 0;
}
