#include "oep_core.h"
#include "oep_order_parameter.h"
#include "oep_phase_transition.h"
#include <stdio.h>
#include <math.h>

/* ==============================================================
 * Example 1: Supercritical Pitchfork Bifurcation
 *
 * Demonstrates the canonical order parameter equation:
 *   dξ/dt = αξ - βξ³
 *
 * As the control parameter α crosses zero from negative to positive:
 *   α < 0: single stable equilibrium at ξ = 0 (disordered phase)
 *   α = 0: critical point, onset of ordering
 *   α > 0: two stable equilibria at ξ = ±√(α/β) (ordered phase)
 *
 * This is the simplest example of spontaneous symmetry breaking
 * and order parameter emergence in synergetic systems.
 * ============================================================== */

int main(void) {
    printf("=== Example 1: Pitchfork Bifurcation ===\n\n");

    double beta = 1.0;

    /* Sweep control parameter α from -2 to 2 */
    printf("Control parameter α sweep:\n");
    printf("  α        ξ_eq(stable)    λ(eigenvalue)    Phase\n");
    printf("  -------- ---------------- ---------------- --------\n");

    for (int i = -20; i <= 20; i++) {
        double alpha = (double)i * 0.1;
        double eq0, eq1, eq2;
        int n_eq;
        oep_op_scalar_equilibria(alpha, beta, &eq0, &eq1, &eq2, &n_eq);

        /* Find stable equilibria */
        if (n_eq == 1) {
            /* α ≤ 0: only ξ = 0 is stable */
            double lambda = oep_op_scalar_stability(0.0, alpha, beta);
            printf("  %8.2f  %16.4f  %16.4f  DISORDER\n", alpha, 0.0, lambda);
        } else {
            /* α > 0: ξ = 0 is unstable, ±√(α/β) are stable */
            double xi_stable = eq1;  /* Positive branch */
            double lambda = oep_op_scalar_stability(xi_stable, alpha, beta);
            printf("  %8.2f  %16.4f  %16.4f  ORDER\n", alpha, xi_stable, lambda);
        }
    }

    /* Evolve order parameter at different α values */
    printf("\nOrder parameter evolution (ξ₀ = 0.1):\n\n");

    double alpha_values[] = {-0.5, 0.0, 0.5, 1.0};
    int n_steps = 100;
    double dt = 0.05;

    for (int a = 0; a < 4; a++) {
        double alpha = alpha_values[a];
        OEPScalarOP* op = oep_op_scalar_create(alpha, beta, 0.0);
        op->value = 0.1;

        double* traj = malloc((size_t)n_steps * sizeof(double));
        oep_op_scalar_evolve(op, dt, n_steps, traj);

        printf("  α = %5.1f: ξ(t) → %8.4f", alpha, traj[n_steps - 1]);
        if (alpha > 0.0) {
            double xi_eq = sqrt(alpha / beta);
            printf("  (expected: ±%8.4f)", xi_eq);
        }
        printf("\n");

        free(traj);
        oep_op_scalar_free(op);
    }

    /* Show potential landscape */
    printf("\nLandau potential V(ξ) = -αξ²/2 + βξ⁴/4:\n");
    printf("  For α = 1.0, β = 1.0:\n");
    for (int i = -30; i <= 30; i++) {
        double xi = (double)i * 0.1;
        double V = oep_op_scalar_potential(xi, 1.0, 1.0);
        printf("    V(%5.1f) = %8.4f\n", xi, V);
    }

    /* Bifurcation diagram generation */
    printf("\nGenerating bifurcation diagram...\n");
    OEPBifDiag* bd = oep_bifdiag_pitchfork(-2.0, 2.0, 41, 0.0, 1.0, beta);
    printf("  Total branches: %d\n", bd->n_branches);
    printf("  Branch 0 (ξ=0): first point p=%.2f, ξ=%.4f\n",
           bd->branches[0][0].param, bd->branches[0][0].equilibrium);
    printf("  Branch 0 (ξ=0): last  point p=%.2f, ξ=%.4f\n",
           bd->branches[0][40].param, bd->branches[0][40].equilibrium);
    oep_bifdiag_free(bd);

    printf("\n=== Example 1 Complete ===\n");
    return 0;
}
