/**
 * example_giere_meinhardt.c — Gierer-Meinhardt Activator-Inhibitor Patterns
 *
 * Demonstrates the classical activator-inhibitor model for biological
 * morphogenesis. Patterns form via "short-range activation, long-range
 * inhibition" when the inhibitor diffusion rate greatly exceeds the
 * activator diffusion rate.
 *
 * Reference: Gierer, A. & Meinhardt, H. (1972)
 *            "A Theory of Biological Pattern Formation"
 *            Kybernetik 12(1):30-39
 *
 * Knowledge Coverage: L6 — Canonical problem: Gierer-Meinhardt morphogenesis
 */

#include "turing_common.h"
#include "turing_models.h"
#include "turing_analysis.h"
#include "turing_visualize.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main(void) {
    printf("=== Gierer-Meinhardt Activator-Inhibitor Pattern Formation ===\n\n");

    int N = 128;
    double L = 2.0;

    ModelParams params;
    model_params_default(MODEL_GIERER_MEINHARDT, &params);
    params.a     = 0.1;
    params.b     = 1.0;
    params.c     = 0.001;
    params.rho_a = 1.0;
    params.rho_h = 1.0;
    params.mu_a  = 0.01;
    params.mu_h  = 0.01;

    /* Key: Dv >> Du — inhibitor diffuses much faster */
    double Du = 0.002;
    double Dv = 0.8;
    double diffusion_ratio = Dv / Du;

    printf("Model: Gierer-Meinhardt (Activator-Inhibitor)\n");
    printf("Parameters:\n");
    printf("  a=%.4f, mu_a=%.4f, rho_a=%.4f\n", params.a, params.mu_a, params.rho_a);
    printf("  b=%.4f, mu_h=%.4f, rho_h=%.4f\n\n", params.b, params.mu_h, params.rho_h);
    printf("Diffusion: Du=%.4f, Dv=%.4f (Dv/Du=%.1f)\n",
           Du, Dv, diffusion_ratio);
    printf("  ➜ Short-range activation (low Du) + Long-range inhibition (high Dv)\n");
    printf("Grid: %dx%d, Domain: %.2fx%.2f\n", N, N, L, L);

    /* Stability analysis */
    printf("\n--- Turing Instability Analysis ---\n");
    HomogeneousSteadyState hss;
    double u_guess = 1.0, v_guess = 1.0;

    hss.converged = 0;
    hss.u_star = u_guess;
    hss.v_star = v_guess;

    /* Quick Newton solve from initial guess */
    for (int iter = 0; iter < 100; iter++) {
        double f, g, J[2][2];
        giere_meinhardt_reaction(hss.u_star, hss.v_star, &params, &f, &g);
        giere_meinhardt_jacobian(hss.u_star, hss.v_star, &params, J);

        double det = J[0][0]*J[1][1] - J[0][1]*J[1][0];
        if (fabs(det) < 1e-15) break;

        double du = -( J[1][1]*f - J[0][1]*g) / det;
        double dv = -(-J[1][0]*f + J[0][0]*g) / det;
        hss.u_star += du;
        hss.v_star += dv;

        if (fabs(du) < 1e-12 && fabs(dv) < 1e-12) {
            hss.converged = 1;
            break;
        }
    }

    if (hss.converged) {
        printf("HSS: u*=%.6f, v*=%.6f\n", hss.u_star, hss.v_star);

        double J[2][2];
        giere_meinhardt_jacobian(hss.u_star, hss.v_star, &params, J);
        printf("Jacobian: [[%.4f, %.4f], [%.4f, %.4f]]\n",
               J[0][0], J[0][1], J[1][0], J[1][1]);
        printf("  f_u=%.4f > 0 → autocatalytic activation\n", J[0][0]);
        printf("  f_v=%.4f < 0 → inhibitor suppresses activator\n", J[0][1]);
        printf("  g_v=%.4f → inhibitor decay\n\n", J[1][1]);

        TuringConditions tc;
        turing_conditions_compute(J, Du, Dv, &tc);
        turing_conditions_print(&tc);

        if (tc.turing_unstable) {
            DispersionRelation *dr = turing_dispersion_compute(J, Du, Dv,
                200, 0.1, 20.0);
            if (dr) {
                printf("\nDispersion Relation:\n");
                printf("  λ_max = %.6f (at k_c = %.4f)\n",
                       dr->lambda_max, dr->k_critical);
                printf("  Predicted wavelength: %.4f\n",
                       dr->predicted_wavelength);
                printf("  Unstable wavenumber band: [%.4f, %.4f]\n",
                       dr->k_min, dr->k_max);
                dispersion_free(dr);
            }
        } else {
            printf("\n⚠  NOT Turing-unstable with these parameters.\n");
            printf("   Consider increasing Dv/Du ratio or adjusting a, mu_a.\n");
        }
    } else {
        printf("Warning: HSS not converged. Using u=1.0, v=1.0 as initial.\n");
        hss.u_star = 1.0;
        hss.v_star = 1.0;
    }

    printf("\n--- Simulation ---\n");

    ScalarField2D *u = field_alloc(N, N, L, L, 1);
    ScalarField2D *v = field_alloc(N, N, L, L, 1);
    ScalarField2D *lap_u = field_alloc(N, N, L, L, 1);
    ScalarField2D *lap_v = field_alloc(N, N, L, L, 1);

    if (!u || !v || !lap_u || !lap_v) {
        printf("ERROR: Allocation failed\n");
        return 1;
    }

    /* Random perturbation around HSS */
    field_init_random(u, hss.u_star, 0.05 * hss.u_star, 42);
    field_init_random(v, hss.v_star, 0.05 * hss.v_star, 142);
    field_apply_bc(u, BC_NEUMANN, BC_NEUMANN);
    field_apply_bc(v, BC_NEUMANN, BC_NEUMANN);

    double dt = 0.01;
    int max_steps = 2000;
    int save_interval = 200;

    printf("dt=%.4f, max_steps=%d\n", dt, max_steps);

    clock_t t_start = clock();

    for (int step = 0; step <= max_steps; step++) {
        field_apply_bc(u, BC_NEUMANN, BC_NEUMANN);
        field_apply_bc(v, BC_NEUMANN, BC_NEUMANN);
        laplacian_compute(u, lap_u, LAPLACIAN_5PT);
        laplacian_compute(v, lap_v, LAPLACIAN_5PT);

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                double u_val = field_get(u, i, j);
                double v_val = field_get(v, i, j);
                double f, g;
                giere_meinhardt_reaction(u_val, v_val, &params, &f, &g);

                field_set(u, i, j,
                    u_val + dt * (f + Du * field_get(lap_u, i, j)));
                field_set(v, i, j,
                    v_val + dt * (g + Dv * field_get(lap_v, i, j)));
            }
        }

        if (step % save_interval == 0) {
            PatternMetrics pm_u = field_compute_metrics(u);
            double t = step * dt;
            printf("  Step %5d, t=%.3f: mean(u)=%.4f, amp(u)=%.4f, type=%s\n",
                   step, t, pm_u.mean, pm_u.pattern_amplitude,
                   pattern_type_name(pm_u.classified_type));

            char fname[256];
            snprintf(fname, sizeof(fname),
                     "build/gm_%05d.ppm", step / save_interval);
            field_write_ppm_pair(u, v, fname, COLORMAP_VIRIDIS, 0.0, 2.0, 0.0, 5.0);
        }
    }

    clock_t t_end = clock();
    printf("\nSimulation complete: %.2f CPU seconds\n",
           (double)(t_end - t_start) / CLOCKS_PER_SEC);

    /* Final analysis */
    PatternMetrics pm_u = field_compute_metrics(u);
    PatternMetrics pm_v = field_compute_metrics(v);
    printf("\n--- Final Pattern ---\n");
    printf("Activator (u):\n"); pattern_metrics_print(&pm_u);
    printf("\nInhibitor (v):\n"); pattern_metrics_print(&pm_v);

    field_write_ppm_pair(u, v, "build/gm_final.ppm",
                         COLORMAP_VIRIDIS, 0.0, 2.0, 0.0, 5.0);
    printf("\nFinal state: build/gm_final.ppm\n");

    field_free(u); field_free(v);
    field_free(lap_u); field_free(lap_v);

    printf("Done.\n");
    return 0;
}
