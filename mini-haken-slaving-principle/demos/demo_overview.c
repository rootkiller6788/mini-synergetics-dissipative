/* ==============================================================
 * demo_overview.c — Interactive Demo of Haken's Slaving Principle
 *
 * Demonstrates all major features of the library:
 *   1. System creation and configuration
 *   2. Eigenmode computation and classification
 *   3. Spectral gap analysis
 *   4. Adiabatic elimination
 *   5. Ginzburg-Landau parameter fitting
 *   6. Visual output of free energy landscape
 *
 * Run with: ./build/demo_demo_overview
 * ============================================================== */

#include "haken_core.h"
#include "haken_mode_decomp.h"
#include "haken_adiabatic.h"
#include "haken_order_param.h"
#include "haken_stability.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static void print_separator(const char* title) {
    printf("\n========================================\n");
    printf("  %s\n", title);
    printf("========================================\n");
}

int main(void) {
    printf("╔══════════════════════════════════════════════╗\n");
    printf("║  Haken Slaving Principle — Interactive Demo  ║\n");
    printf("║  Hermann Haken (1977) Synergetics           ║\n");
    printf("╚══════════════════════════════════════════════╝\n");

    /* ---- Part 1: System Definition ---- */
    print_separator("Part 1: System Definition");

    Haken_System* sys = haken_system_create(3, NULL, NULL);
    printf("Created 3-dimensional synergetic system.\n");
    printf("  State dimension: %d\n", sys->n_dim);
    printf("  Control parameter: α = %.3f\n", sys->control_param);

    /* ---- Part 2: Linear Stability ---- */
    print_separator("Part 2: Linear Stability Analysis");

    /* Set up a system with one critical mode and two stable modes */
    double L[9] = {
        -0.02, 0.0,  0.0,    /* λ₁ ≈ -0.02 (critical/slow) */
         0.0, -2.0,  0.0,    /* λ₂ = -2.0 (stable/fast) */
         0.0,  0.0, -5.0     /* λ₃ = -5.0 (stable/very fast) */
    };
    memcpy(sys->linear_matrix, L, 9 * sizeof(double));
    haken_compute_eigenmodes(sys);

    printf("Eigenvalue spectrum:\n");
    for (int k = 0; k < 3; k++) {
        printf("  λ%d = %+.4f  ", k, sys->modes[k].lambda_re);
        printf("τ%d = ", k);
        if (sys->modes[k].relaxation_time < INFINITY)
            printf("%.4f", sys->modes[k].relaxation_time);
        else
            printf("∞");
        printf("\n");
    }

    haken_classify_modes(sys, 0.1, 0.5);
    printf("\nMode classification (ε_crit=0.1, ε_stable=0.5):\n");
    printf("  Critical modes: %d\n", sys->n_critical);
    printf("  Stable modes:   %d\n", sys->n_stable);
    printf("  Unstable modes: %d\n", sys->n_unstable);

    /* ---- Part 3: Spectral Gap ---- */
    print_separator("Part 3: Spectral Gap Analysis");

    double gap = haken_spectral_gap(sys);
    printf("Spectral gap γ = %.4f\n", gap);
    printf("Timescale ratio = %.2f\n", gap / fabs(sys->modes[1].lambda_re));
    if (gap > 0.5) {
        printf("✓ Gap sufficient for slaving principle\n");
    } else {
        printf("✗ Gap too small — slaving may not hold\n");
    }

    /* ---- Part 4: Adiabatic Elimination ---- */
    print_separator("Part 4: Adiabatic Elimination");

    Haken_AdiabaticConfig config = {
        .order = HAKEN_ADIABATIC_LINEAR,
        .tolerance = 1e-8,
        .max_iterations = 50,
        .check_validity = true,
        .max_allowed_error = 0.5
    };

    Haken_SlavingResult result;
    int ret = haken_adiabatic_eliminate(sys, &config, &result);

    if (ret == 0) {
        printf("Adiabatic elimination successful.\n");
        printf("  Slaving holds: %s\n", result.slaving_holds ? "YES ✓" : "NO ✗");
        printf("  Adiabatic error: %.6f\n", result.adiabatic_error);
        printf("  Fast modes slaved: %d\n", result.slaving.n_fast);
        printf("  Slow (order param) modes: %d\n", result.slaving.n_slow);
    } else {
        printf("Adiabatic elimination failed (code %d)\n", ret);
    }

    /* ---- Part 5: Order Parameter Identification ---- */
    print_separator("Part 5: Order Parameter Identification");

    Haken_OrderParam op;
    haken_identify_spectral(sys, 2, 0.1, &op);
    printf("Identified %d order parameter(s):\n", op.n_op);
    for (int a = 0; a < op.n_op; a++) {
        printf("  OP%d: λ_eff = %.6f, amplitude = %.6f\n",
               a, op.effective_lambda[a], op.amplitude[a]);
    }

    double var_expl = haken_variance_explained(sys, &op);
    printf("Variance explained by OPs: %.1f%%\n", var_expl * 100.0);

    /* ---- Part 6: Ginzburg-Landau Free Energy ---- */
    print_separator("Part 6: GL Free Energy Landscape");

    Haken_GLParameters gl_params = {
        .a = 0.02,
        .b = 1.0,
        .critical_alpha = 0.0,
        .xi0 = 1.0,
        .tau0 = 1.0,
        .nu = 0.5,
        .z = 2.0,
        .beta = 0.5
    };

    printf("GL parameters: a=%.2f, b=%.1f\n", gl_params.a, gl_params.b);

    double xi_grid[21], F_grid[21];
    haken_gl_potential_landscape(&gl_params, -1.0, 1.0, 21, xi_grid, F_grid);

    printf("Free energy F(ξ):\n");
    for (int i = 0; i < 21; i++) {
        printf("  ξ = %+5.2f  →  F(ξ) = %+8.6f", xi_grid[i], F_grid[i]);
        /* Visual bar */
        int bar_len = (int)(fabs(F_grid[i]) * 40);
        if (bar_len > 40) bar_len = 40;
        printf("  ");
        for (int b = 0; b < bar_len; b++) printf("█");
        printf("\n");
    }

    /* Fixed points */
    double fixed[3];
    int n_fixed = haken_gl_fixed_points(&gl_params, fixed);
    printf("\nGL fixed points (a=%.2f): %d\n", gl_params.a, n_fixed);
    for (int k = 0; k < n_fixed; k++) {
        double stab_val;
        haken_gl_stability(&gl_params, &fixed[k], 1, &stab_val);
        printf("  ξ* = %+.4f  (λ = %+.4f, %s)\n",
               fixed[k], stab_val,
               stab_val < 0 ? "stable" : "unstable");
    }

    /* ---- Part 7: Summary ---- */
    print_separator("Summary");

    printf("The slaving principle in action:\n");
    printf("  1. Near bifurcation (α → α_c): one mode slows down (λ → 0)\n");
    printf("  2. Spectral gap γ = %.2f separates slow from fast modes\n", gap);
    printf("  3. Fast modes are adiabatically eliminated: q_f = h(q_s)\n");
    printf("  4. Reduced GL equation: dξ/dt = a·ξ - b·ξ³\n");
    printf("  5. Order parameter ξ governs macroscopic behavior\n");

    haken_system_free(sys);
    printf("\n╔══════════════════════════════════════════════╗\n");
    printf("║              Demo Complete ✓                  ║\n");
    printf("╚══════════════════════════════════════════════╝\n");

    return 0;
}
