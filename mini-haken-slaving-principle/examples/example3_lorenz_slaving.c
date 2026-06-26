/* ==============================================================
 * example3_lorenz_slaving.c — Lorenz Model Near First Instability
 *
 * Edward Lorenz (1963) discovered deterministic chaos in a
 * simplified convection model. Near the first instability
 * (r → 1+), the Lorenz system reduces to a single order
 * parameter obeying a Ginzburg-Landau equation via Haken's
 * slaving principle.
 *
 * Full Lorenz equations:
 *   ẋ = σ(y - x)
 *   ẏ = r·x - y - xz
 *   ż = xy - b·z
 *
 * Parameters: σ = 10 (Prandtl), b = 8/3 (geometry), r (Rayleigh)
 *
 * Haken (1975) showed that near r=1, x and y are order parameters
 * (O(ε)) while z is slaved (O(ε²)). The reduced equation is a GL
 * form: dξ/dt = (r-1)ξ - ξ³.
 *
 * This example:
 *   1. Creates the Lorenz system at r = 1.1 (just above onset)
 *   2. Analyzes eigenmodes and identifies critical mode(s)
 *   3. Performs adiabatic elimination
 *   4. Compares GL prediction with full dynamics
 *
 * Reference: Haken (1975) Phys. Lett. 53A, 77-78
 * ============================================================== */

#include "haken_core.h"
#include "haken_dynamics.h"
#include "haken_adiabatic.h"
#include "haken_order_param.h"
#include "haken_stability.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void) {
    printf("=== Lorenz Model — Slaving Near First Instability ===\n\n");

    double sigma = 10.0;
    double b_geo = 8.0 / 3.0;
    double r = 1.1;  /* Just above onset (r_c = 1) */

    printf("Lorenz parameters: σ=%.1f, b=%.4f, r=%.2f\n\n", sigma, b_geo, r);

    Haken_System* sys = haken_create_lorenz_instability(r, sigma, b_geo);
    if (!sys) { printf("Failed to create system\n"); return 1; }

    /* Compute eigenmodes at the origin (fixed point) */
    haken_update_linear_matrix(sys);
    haken_compute_eigenmodes(sys);

    printf("Eigenvalues at origin (r=%.2f):\n", r);
    for (int k = 0; k < sys->n_dim; k++) {
        printf("  λ%d = %+.6f", k, sys->modes[k].lambda_re);
        if (fabs(sys->modes[k].lambda_im) > 1e-6)
            printf(" ± %.6fi", fabs(sys->modes[k].lambda_im));
        printf("  τ = %.4f\n", sys->modes[k].relaxation_time);
    }

    /* Classify modes */
    haken_classify_modes(sys, 0.5, 1.0);
    printf("\nMode classification:\n");
    printf("  Critical (order parameter) modes: %d\n", sys->n_critical);
    printf("  Stable (slaved) modes: %d\n", sys->n_stable);
    printf("  Unstable modes: %d\n", sys->n_unstable);

    /* Spectral gap */
    double gap = haken_spectral_gap(sys);
    printf("\nSpectral gap γ = %.4f\n", gap);
    printf("Near r_c=1, the most unstable eigenvalue is λ₁ ≈ r-1 = %.4f\n", r - 1.0);
    printf("The remaining eigenvalues are O(1): λ₂ ≈ -%.1f, λ₃ ≈ -%.1f\n",
           sigma + 1.0, b_geo);

    /* GL prediction for the order parameter */
    printf("\nGinzburg-Landau reduced equation:\n");
    printf("  dξ/dt = %.4f·ξ - (1/σ)·ξ³ = %.4f·ξ - %.4f·ξ³\n",
           r - 1.0, r - 1.0, 1.0 / sigma);
    printf("  Stationary amplitude: ξ_eq = ±√(σ(r-1)) = ±%.4f\n",
           sqrt(sigma * (r - 1.0)));

    /* Numerical integration */
    printf("\nTime evolution from small initial perturbation:\n");
    double q0[3] = {0.01, 0.01, 0.0};
    haken_set_state(sys, q0);
    printf("  Initial: (%.6f, %.6f, %.6f)\n", q0[0], q0[1], q0[2]);
    printf("  %8s  %12s  %12s  %12s\n", "t", "x", "y", "z");

    for (int step = 0; step <= 100; step++) {
        if (step % 20 == 0) {
            printf("  %8.3f  %12.6f  %12.6f  %12.6f\n",
                   0.01 * step, sys->state[0], sys->state[1], sys->state[2]);
        }
        if (step < 100) haken_step_rk4(sys, 0.01);
    }

    /* Adiabatic elimination at the fixed point */
    Haken_AdiabaticConfig config = {HAKEN_ADIABATIC_LINEAR, 1e-8, 50, true, 0.5};
    Haken_SlavingResult result;
    int ret = haken_adiabatic_eliminate(sys, &config, &result);
    printf("\nAdiabatic elimination: %s\n", (ret == 0) ? "OK" : "FAILED");
    if (ret == 0) {
        printf("  Slaving valid: %s\n", result.slaving_holds ? "YES" : "NO");
        printf("  Adiabatic error: %.6f\n", result.adiabatic_error);
        printf("  Bifurcation type: %d\n", result.bif_type);
    }

    /* Check mean-field order parameter */
    double mf = haken_mean_field_op(sys->state, 3);
    printf("\nMean-field order parameter: %.6f\n", mf);

    haken_system_free(sys);
    printf("\n=== Done ===\n");
    return 0;
}
