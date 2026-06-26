/**
 * @file example_spiral_wave.c
 * @brief End-to-end example: BZ spiral wave simulation in 2D
 *
 * Demonstrates:
 *   - Setting up a 2D reaction-diffusion domain
 *   - Seeding a spiral wave
 *   - Simulating wave propagation with forward Euler
 *   - Detecting the spiral tip
 *
 * Knowledge: L6 Canonical Problem - BZ spiral wave
 *
 * Compile: gcc -I../include -o example_spiral \
 *          example_spiral_wave.c ../src/bz_reaction_diffusion.c \
 *          ../src/bz_reaction.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/bz_reaction_diffusion.h"

int main(void)
{
    bz_rd_2d_t domain;
    bz_barkley_params_t bp;
    double *u, *v, *u_new, *v_new;
    int nx = 80, ny = 80;
    int n_steps = 200;
    int step;
    double dt = 0.05;
    double tip_x, tip_y;
    int n_total;

    printf("=== BZ Spiral Wave Simulation (2D) ===\n\n");

    /* Setup domain */
    domain.nx = nx;
    domain.ny = ny;
    domain.Lx = 10.0;
    domain.Ly = 10.0;
    domain.dx = domain.Lx / (nx - 1);
    domain.dy = domain.Ly / (ny - 1);
    domain.Du = 1.0;    /* activator diffusion */
    domain.Dv = 0.0;    /* inhibitor does not diffuse in Barkley */
    domain.bc = BZ_BC_NEUMANN;

    n_total = nx * ny;
    u = (double*)malloc(n_total * sizeof(double));
    v = (double*)malloc(n_total * sizeof(double));
    u_new = (double*)malloc(n_total * sizeof(double));
    v_new = (double*)malloc(n_total * sizeof(double));

    if (!u || !v || !u_new || !v_new) {
        printf("Memory allocation failed!\n");
        return 1;
    }

    /* Barkley parameters for excitable media */
    bp.a = 0.75;
    bp.b = 0.02;
    bp.epsilon = 0.02;

    printf("Domain: %d x %d grid, dx=%.4f, dt=%.4f\n", nx, ny, domain.dx, dt);
    printf("Barkley params: a=%.2f, b=%.3f, eps=%.3f\n", bp.a, bp.b, bp.epsilon);

    /* Initialize with resting state */
    {
        int i;
        for (i = 0; i < n_total; i++) {
            u[i] = 0.0;
            v[i] = 0.0;
        }
    }

    /* Seed spiral wave */
    bz_spiral_seed_2d(u, v, nx, ny, nx/2, ny/3, 10);
    printf("Spiral seed placed at (%d, %d)\n", nx/2, ny/3);

    /* Simulate */
    printf("\nSimulating %d steps...\n", n_steps);
    for (step = 0; step < n_steps; step++) {
        bz_rd_2d_step(u_new, v_new, u, v, &domain, &bp, dt);

        /* Swap buffers */
        {
            double *tmp_u, *tmp_v;
            tmp_u = u; u = u_new; u_new = tmp_u;
            tmp_v = v; v = v_new; v_new = tmp_v;
        }

        /* Report every 50 steps */
        if (step % 50 == 0) {
            /* Find spiral tip */
            int found = bz_find_spiral_tip(&tip_x, &tip_y, u, v, nx, ny);
            if (found) {
                printf("  Step %3d: spiral tip at (%.1f, %.1f)\n",
                       step, tip_x, tip_y);
            }

            /* Count active (excited) cells */
            {
                int i, active = 0;
                for (i = 0; i < n_total; i++) {
                    if (u[i] > 0.5) active++;
                }
                printf("           active cells: %d / %d (%.1f%%)\n",
                       active, n_total, 100.0*active/n_total);
            }
        }
    }

    /* Final analysis */
    {
        int found = bz_find_spiral_tip(&tip_x, &tip_y, u, v, nx, ny);
        printf("\n=== Final State ===\n");
        if (found) {
            printf("Spiral tip location: (%.2f, %.2f)\n", tip_x, tip_y);
        } else {
            printf("No spiral tip detected.\n");
        }

        /* Count spiral arms */
        if (found) {
            int arms = bz_spiral_arm_count(tip_x, tip_y, u, nx, ny);
            printf("Spiral arms: %d\n", arms);
        }

        /* u field statistics */
        {
            int i;
            double u_max = u[0], u_min = u[0], u_avg = 0.0;
            for (i = 0; i < n_total; i++) {
                u_avg += u[i];
                if (u[i] > u_max) u_max = u[i];
                if (u[i] < u_min) u_min = u[i];
            }
            u_avg /= n_total;
            printf("u field: min=%.4f, max=%.4f, avg=%.4f\n",
                   u_min, u_max, u_avg);
        }
    }

    free(u); free(v); free(u_new); free(v_new);
    printf("\nDone.\n");
    return 0;
}
