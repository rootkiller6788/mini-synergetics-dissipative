/**
 * example_grayscott.c — Gray-Scott Pattern Formation Simulation
 *
 * Demonstrates the Gray-Scott reaction-diffusion model on a 2D grid.
 * Produces spot or stripe patterns depending on (F, k) parameters.
 *
 * Reference: Pearson, J.E. (1993) "Complex Patterns in a Simple System"
 *            Science 261:189-192
 *
 * Knowledge Coverage: L6 — Canonical problem: Gray-Scott pattern formation
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
    printf("=== Gray-Scott Turing Pattern Simulation ===\n\n");

    /* Domain: 128×128 grid, physical size 2.5×2.5 */
    int N = 128;
    double L = 2.5;

    /* Parameters: spot-forming regime (Pearson α-type)
     * F=0.04, k=0.06 produces stable spots
     * Du=2e-5, Dv=1e-5 with Dv/Du=0.5 (ratio tuned for spots) */
    ModelParams params;
    model_params_default(MODEL_GRAY_SCOTT, &params);
    params.F = 0.04;
    params.k = 0.06;

    double Du = 2.0e-5;
    double Dv = 1.0e-5;

    printf("Model: Gray-Scott\n");
    printf("Parameters: F=%.4f, k=%.4f\n", params.F, params.k);
    printf("Diffusion: Du=%.2e, Dv=%.2e (ratio=%.2f)\n", Du, Dv, Dv/Du);
    printf("Grid: %dx%d, Domain: %.2fx%.2f\n", N, N, L, L);

    /* Step 1: Linear stability analysis */
    printf("\n--- Linear Stability Analysis ---\n");
    HomogeneousSteadyState hss;
    if (grayscott_find_hss(&params, &hss) == 0 && hss.converged) {
        printf("HSS: u*=%.6f, v*=%.6f\n", hss.u_star, hss.v_star);

        double J[2][2];
        grayscott_jacobian(hss.u_star, hss.v_star, &params, J);
        printf("Jacobian at HSS: [[%.4f, %.4f], [%.4f, %.4f]]\n",
               J[0][0], J[0][1], J[1][0], J[1][1]);

        TuringConditions tc;
        turing_conditions_compute(J, Du, Dv, &tc);
        turing_conditions_print(&tc);

        if (tc.turing_unstable) {
            DispersionRelation *dr = turing_dispersion_compute(J, Du, Dv,
                200, 0.01, 10.0);
            if (dr) {
                turing_dispersion_print(dr);
                printf("  Predicted pattern wavelength: %.4f\n",
                       dr->predicted_wavelength);
                dispersion_free(dr);
            }
        }
    }

    /* Step 2: Simulate */
    printf("\n--- Simulation ---\n");

    /* Allocate fields */
    ScalarField2D *u = field_alloc(N, N, L, L, 1);
    ScalarField2D *v = field_alloc(N, N, L, L, 1);
    ScalarField2D *lap_u = field_alloc(N, N, L, L, 1);
    ScalarField2D *lap_v = field_alloc(N, N, L, L, 1);

    if (!u || !v || !lap_u || !lap_v) {
        printf("ERROR: Failed to allocate fields\n");
        return 1;
    }

    /* Initialize: HSS + small random perturbation */
    if (hss.converged) {
        field_init_random(u, hss.u_star, 0.01 * hss.u_star, 12345);
        field_init_random(v, hss.v_star, 0.01 * hss.v_star, 67890);
    } else {
        field_fill_interior(u, 1.0);
        field_fill_interior(v, 0.0);
        field_init_random(u, 1.0, 0.01, 12345);
        field_init_random(v, 0.0, 0.01, 67890);
    }

    field_apply_bc(u, BC_NEUMANN, BC_NEUMANN);
    field_apply_bc(v, BC_NEUMANN, BC_NEUMANN);

    /* Time stepping */
    double dt = 0.5;
    int max_steps = 5000;
    int save_interval = 500;
    double dx = L / (N - 1);

    printf("Time step: dt=%.4f, max steps: %d\n", dt, max_steps);
    printf("CFL stability limit: %.4f\n", 0.5 * dx * dx / (2.0 * Dv));

    clock_t t_start = clock();

    for (int step = 0; step <= max_steps; step++) {
        /* Compute Laplacians */
        field_apply_bc(u, BC_NEUMANN, BC_NEUMANN);
        field_apply_bc(v, BC_NEUMANN, BC_NEUMANN);
        laplacian_compute(u, lap_u, LAPLACIAN_5PT);
        laplacian_compute(v, lap_v, LAPLACIAN_5PT);

        /* Forward Euler update */
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                double u_val = field_get(u, i, j);
                double v_val = field_get(v, i, j);
                double lap_u_val = field_get(lap_u, i, j);
                double lap_v_val = field_get(lap_v, i, j);

                double f, g;
                grayscott_reaction(u_val, v_val, &params, &f, &g);

                field_set(u, i, j, u_val + dt * (f + Du * lap_u_val));
                field_set(v, i, j, v_val + dt * (g + Dv * lap_v_val));
            }
        }

        /* Output and diagnostics */
        if (step % save_interval == 0) {
            /* Compute pattern metrics */
            PatternMetrics pm = field_compute_metrics(u);

            double t = step * dt;
            printf("  Step %5d, t=%.4f: mean(u)=%.4f, amp=%.4f, type=%s\n",
                   step, t, pm.mean, pm.pattern_amplitude,
                   pattern_type_name(pm.classified_type));

            /* Save PPM image */
            char filename[256];
            snprintf(filename, sizeof(filename),
                     "build/grayscott_%05d.ppm", step / save_interval);
            field_write_ppm(u, filename, COLORMAP_VIRIDIS, 0.0, 1.0);
        }
    }

    clock_t t_end = clock();
    double cpu_time = (double)(t_end - t_start) / CLOCKS_PER_SEC;
    printf("\nSimulation complete: %.2f CPU seconds\n", cpu_time);

    /* Final pattern analysis */
    printf("\n--- Final Pattern Analysis ---\n");
    PatternMetrics pm_final = field_compute_metrics(u);
    pattern_metrics_print(&pm_final);

    /* Save final state */
    field_write_ppm(u, "build/grayscott_final.ppm", COLORMAP_VIRIDIS, 0.0, 1.0);
    printf("\nFinal state saved to build/grayscott_final.ppm\n");

    field_free(u);
    field_free(v);
    field_free(lap_u);
    field_free(lap_v);

    printf("\nDone.\n");
    return 0;
}
