/* ==============================================================
 * example1_haken_prototype.c — The Simplest Slaving Example
 *
 * Haken (1977) §5.1 prototype system:
 *   dq1/dt = α·q1 - a·q1·q2     (order parameter, slow)
 *   dq2/dt = -β·q2 + b·q1²       (slaved mode, fast)
 *
 * Near α ≈ 0, q1 is the order parameter and q2 adiabatically
 * follows: q2 ≈ (b/β)·q1².
 *
 * Substituting: dq1/dt = α·q1 - (ab/β)·q1³  (Ginzburg-Landau)
 *
 * This example demonstrates:
 *   1. Timescale separation (α ≈ 0 vs β ≫ 0)
 *   2. Adiabatic elimination
 *   3. Supercritical pitchfork bifurcation at α = 0
 * ============================================================== */

#include "haken_core.h"
#include "haken_dynamics.h"
#include "haken_adiabatic.h"
#include "haken_order_param.h"
#include "haken_stability.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Prototype model with explicit RHS */
static void proto_rhs(const double* q, int n, double* dq) {
    /* Parameters */
    double alpha = 0.1;  /* Near critical */
    double beta  = 2.0;
    double a     = 1.0;
    double b     = 1.0;

    dq[0] = alpha * q[0] - a * q[0] * q[1];
    dq[1] = -beta * q[1] + b * q[0] * q[0];
    (void)n;
}

int main(void) {
    printf("=== Haken Slaving Principle — Prototype Model ===\n\n");
    printf("System:  dq1/dt = α·q1 - a·q1·q2\n");
    printf("         dq2/dt = -β·q2 + b·q1²\n\n");

    double alpha = 0.1;
    double beta  = 2.0;
    double a     = 1.0;
    double b_coeff = 1.0;

    /* Create system */
    Haken_System* sys = haken_create_prototype_model(alpha, beta, a, b_coeff);
    if (!sys) { printf("Failed to create system\n"); return 1; }

    /* Set initial conditions */
    double q0[2] = {0.2, 0.0};
    haken_set_state(sys, q0);

    printf("Initial state: q1=%.4f, q2=%.4f\n", q0[0], q0[1]);
    printf("Control param: α=%.4f (critical at α_c=0.0)\n", alpha);
    printf("Fast mode β=%.4f (relaxation time τ=%.4f)\n\n", beta, 1.0/beta);

    /* Compute eigenmodes */
    haken_update_linear_matrix(sys);
    haken_compute_eigenmodes(sys);
    haken_classify_modes(sys, 0.1, 0.5);

    printf("Mode analysis:\n");
    for (int k = 0; k < sys->n_dim; k++) {
        printf("  Mode %d: λ=%.4f, τ=%.4f, type=%d\n",
               k, sys->modes[k].lambda_re,
               sys->modes[k].relaxation_time,
               sys->modes[k].type);
    }

    double gap = haken_spectral_gap(sys);
    printf("\nSpectral gap: γ = %.4f\n", gap);
    printf("Timescale ratio: τ_slow/τ_fast = %.4f\n\n", gap / fabs(sys->modes[1].lambda_re));

    /* Adiabatic elimination */
    Haken_SlavingResult result;
    Haken_AdiabaticConfig config = {HAKEN_ADIABATIC_LINEAR, 1e-8, 50, true, 0.1};
    int ret = haken_adiabatic_eliminate(sys, &config, &result);
    if (ret == 0) {
        printf("Adiabatic elimination:\n");
        printf("  Slaving holds: %s\n", result.slaving_holds ? "YES" : "NO");
        printf("  Adiabatic error: %.6f\n", result.adiabatic_error);
        printf("  Spectral gap: %.4f\n", result.spectral_gap);

        /* GL prediction */
        double gl_cubic = a * b_coeff / beta;
        printf("\nGinzburg-Landau reduced equation:\n");
        printf("  dξ/dt = %.4f·ξ - %.4f·ξ³\n", alpha, gl_cubic);
        printf("  Stationary states: ξ = 0, ±√(α·β/(a·b)) = ±%.4f\n",
               sqrt(alpha * beta / (a * b_coeff)));
    }

    /* Integrate the system */
    printf("\nTime evolution (RK4, dt=0.01):\n");
    printf("  %8s  %12s  %12s\n", "t", "q1", "q2");
    for (int step = 0; step <= 50; step++) {
        if (step % 10 == 0) {
            printf("  %8.3f  %12.6f  %12.6f\n",
                   0.01 * step, sys->state[0], sys->state[1]);
        }
        if (step < 50) haken_step_rk4(sys, 0.01);
    }

    /* Order parameter identification */
    Haken_OrderParam op;
    haken_identify_spectral(sys, 1, 0.1, &op);
    printf("\nIdentified order parameter: λ_eff=%.6f\n", op.effective_lambda[0]);

    double variance = haken_variance_explained(sys, &op);
    printf("Variance explained by OP: %.2f%%\n", variance * 100);

    haken_system_free(sys);
    printf("\n=== Done ===\n");
    return 0;
}
