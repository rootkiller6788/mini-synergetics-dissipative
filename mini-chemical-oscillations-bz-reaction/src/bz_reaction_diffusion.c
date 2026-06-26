
/**
 * @file bz_reaction_diffusion.c
 * @brief Reaction-diffusion PDE solvers: 1D/2D Laplacian, Turing instability,
 *        RD time-stepping, spiral wave seeding.
 *
 * Implements:
 *   L3: Finite-difference Laplacian (1D/2D)
 *   L4: Turing instability condition checker
 *   L5: Forward Euler RD time-stepping (1D/2D)
 *   L6: Spiral wave seed, random initialization
 *
 * References:
 *   - Turing (1952) Phil. Trans. R. Soc. B 237, 37-72
 *   - Murray (2002) "Mathematical Biology" Vol II
 *   - Barkley (1991) Physica D 49, 61-70
 */

#include "bz_reaction_diffusion.h"
#include "bz_reaction.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* =========================================================================
 * L3: Finite-difference Laplacian operators
 * ========================================================================= */

void bz_laplacian_1d(double *lap, const double *u, int nx, double dx,
                      bz_rd_boundary_t bc)
{
    int i;
    double inv_dx2;

    if (!lap || !u || nx < 3 || dx <= 0.0) return;

    inv_dx2 = 1.0 / (dx * dx);

    /* Interior points: central difference */
    for (i = 1; i < nx - 1; i++) {
        lap[i] = (u[i+1] - 2.0 * u[i] + u[i-1]) * inv_dx2;
    }

    /* Boundary conditions */
    switch (bc) {
    case BZ_BC_NEUMANN:
        lap[0]    = (u[1] - 2.0 * u[0] + u[1]) * inv_dx2;
        lap[nx-1] = (u[nx-2] - 2.0 * u[nx-1] + u[nx-2]) * inv_dx2;
        break;
    case BZ_BC_DIRICHLET:
        lap[0]    = (u[1] - 2.0 * u[0] + 0.0) * inv_dx2;
        lap[nx-1] = (0.0 - 2.0 * u[nx-1] + u[nx-2]) * inv_dx2;
        break;
    case BZ_BC_PERIODIC:
        lap[0]    = (u[1] - 2.0 * u[0] + u[nx-2]) * inv_dx2;
        lap[nx-1] = (u[0] - 2.0 * u[nx-1] + u[nx-2]) * inv_dx2;
        break;
    default:
        lap[0]    = (u[1] - 2.0 * u[0] + u[1]) * inv_dx2;
        lap[nx-1] = (u[nx-2] - 2.0 * u[nx-1] + u[nx-2]) * inv_dx2;
        break;
    }
}

void bz_laplacian_2d(double *lap, const double *u, int nx, int ny,
                      double dx, bz_rd_boundary_t bc)
{
    int i, j, idx;
    double inv_dx2;
    double u_im1, u_ip1, u_jm1, u_jp1;

    if (!lap || !u || nx < 3 || ny < 3 || dx <= 0.0) return;
    inv_dx2 = 1.0 / (dx * dx);

    for (j = 0; j < ny; j++) {
        for (i = 0; i < nx; i++) {
            idx = j * nx + i;

            if (i == 0) {
                u_im1 = (bc == BZ_BC_PERIODIC) ? u[j*nx + (nx-2)] :
                        (bc == BZ_BC_NEUMANN)   ? u[j*nx + 1] : 0.0;
            } else {
                u_im1 = u[j*nx + (i-1)];
            }

            if (i == nx - 1) {
                u_ip1 = (bc == BZ_BC_PERIODIC) ? u[j*nx + 1] :
                        (bc == BZ_BC_NEUMANN)   ? u[j*nx + (nx-2)] : 0.0;
            } else {
                u_ip1 = u[j*nx + (i+1)];
            }

            if (j == 0) {
                u_jm1 = (bc == BZ_BC_PERIODIC) ? u[(ny-2)*nx + i] :
                        (bc == BZ_BC_NEUMANN)   ? u[1*nx + i] : 0.0;
            } else {
                u_jm1 = u[(j-1)*nx + i];
            }

            if (j == ny - 1) {
                u_jp1 = (bc == BZ_BC_PERIODIC) ? u[1*nx + i] :
                        (bc == BZ_BC_NEUMANN)   ? u[(ny-2)*nx + i] : 0.0;
            } else {
                u_jp1 = u[(j+1)*nx + i];
            }

            lap[idx] = (u_ip1 + u_im1 + u_jp1 + u_jm1 - 4.0 * u[idx]) * inv_dx2;
        }
    }
}

/* =========================================================================
 * L4: Turing instability condition
 * ========================================================================= */

int bz_turing_unstable(double fu, double fv, double gu, double gv,
                        double Du, double Dv)
{
    double trace = fu + gv;
    double determinant = fu * gv - fv * gu;

    if (trace >= 0.0) return 0;
    if (determinant <= 0.0) return 0;
    if (Dv * fu + Du * gv <= 0.0) return 0;

    {
        double cross = Dv * fu + Du * gv;
        if (cross * cross <= 4.0 * Du * Dv * determinant) return 0;
    }
    return 1;
}

double bz_turing_critical_wavenumber(double fu, double fv,
                                      double gu, double gv,
                                      double Du, double Dv)
{
    double det = fu * gv - fv * gu;
    if (det <= 0.0 || Du <= 0.0 || Dv <= 0.0) return 0.0;
    return sqrt(sqrt(det / (Du * Dv)));
}

/* =========================================================================
 * L5: Reaction-diffusion time stepping
 * ========================================================================= */

void bz_rd_1d_step(double *u_new, double *v_new,
                    const double *u, const double *v,
                    const bz_rd_1d_t *domain, double dt,
                    const bz_oregonator_params_t *params)
{
    int i;
    double *lap_u, *lap_v;
    double dudt, dvdt;

    if (!u_new || !v_new || !u || !v || !domain) return;

    lap_u = (double*)malloc(domain->nx * sizeof(double));
    lap_v = (double*)malloc(domain->nx * sizeof(double));

    if (!lap_u || !lap_v) {
        free(lap_u);
        free(lap_v);
        return;
    }

    bz_laplacian_1d(lap_u, u, domain->nx, domain->dx, domain->bc);
    bz_laplacian_1d(lap_v, v, domain->nx, domain->dx, domain->bc);

    for (i = 0; i < domain->nx; i++) {
        bz_tyson_2var_rhs(&dudt, &dvdt, u[i], v[i], params);
        u_new[i] = u[i] + dt * (dudt + domain->Du * lap_u[i]);
        v_new[i] = v[i] + dt * (dvdt + domain->Dv * lap_v[i]);

        if (u_new[i] < 0.0) u_new[i] = 0.0;
        if (v_new[i] < 0.0) v_new[i] = 0.0;
        if (u_new[i] > 1.0) u_new[i] = 1.0;
    }

    free(lap_u);
    free(lap_v);
}

void bz_rd_2d_step(double *u_new, double *v_new,
                    const double *u, const double *v,
                    const bz_rd_2d_t *domain,
                    const bz_barkley_params_t *bp, double dt)
{
    int n, i;
    double *lap_u, *lap_v;
    double u_val, v_val, reaction_u, reaction_v;

    if (!u_new || !v_new || !u || !v || !domain || !bp) return;

    n = domain->nx * domain->ny;
    lap_u = (double*)malloc(n * sizeof(double));
    lap_v = (double*)malloc(n * sizeof(double));

    if (!lap_u || !lap_v) {
        free(lap_u);
        free(lap_v);
        return;
    }

    bz_laplacian_2d(lap_u, u, domain->nx, domain->ny, domain->dx, domain->bc);
    bz_laplacian_2d(lap_v, v, domain->nx, domain->ny, domain->dx, domain->bc);

    for (i = 0; i < n; i++) {
        u_val = u[i];
        v_val = v[i];

        reaction_u = (1.0 / bp->epsilon) * u_val * (1.0 - u_val)
                     * (u_val - (v_val + bp->b) / bp->a);
        reaction_v = u_val - v_val;

        u_new[i] = u_val + dt * (reaction_u + domain->Du * lap_u[i]);
        v_new[i] = v_val + dt * (reaction_v + domain->Dv * lap_v[i]);

        if (u_new[i] < 0.0) u_new[i] = 0.0;
        if (u_new[i] > 1.0) u_new[i] = 1.0;
        if (v_new[i] < 0.0) v_new[i] = 0.0;
    }

    free(lap_u);
    free(lap_v);
}

void bz_spiral_seed_2d(double *u, double *v, int nx, int ny,
                        int cx, int cy, int length)
{
    int i, j, idx;

    if (!u || !v) return;

    for (j = 0; j < ny; j++) {
        for (i = 0; i < nx; i++) {
            idx = j * nx + i;
            if (j < cy) {
                if (i < cx) {
                    u[idx] = 1.0;
                    v[idx] = 0.0;
                } else {
                    u[idx] = 0.0;
                    v[idx] = 0.0;
                }
            } else {
                u[idx] = 0.0;
                v[idx] = 0.0;
            }
        }
    }

    for (i = 0; i < length && (cx + i) < nx; i++) {
        for (j = cy; j < cy + 3 && j < ny; j++) {
            idx = j * nx + (cx + i);
            v[idx] = 1.0;
        }
    }
}

void bz_rd_init_random(double *u, double *v, int n, double low, double high)
{
    int i;
    double range;

    if (!u || !v || n <= 0) return;
    range = high - low;
    for (i = 0; i < n; i++) {
        u[i] = low + range * ((double)rand() / (double)RAND_MAX);
        v[i] = low + range * ((double)rand() / (double)RAND_MAX);
    }
}

/* =========================================================================
 * L6: Spiral tip detection via phase singularity
 * ========================================================================= */

int bz_find_spiral_tip(double *tip_x, double *tip_y,
                        const double *u, const double *v,
                        int nx, int ny)
{
    double max_grad = 0.0;
    int max_i = 0, max_j = 0;
    int i, j, idx;
    double gx, gy, gmag;

    if (!tip_x || !tip_y || !u || !v) return 0;

    for (j = 1; j < ny - 1; j++) {
        for (i = 1; i < nx - 1; i++) {
            idx = j * nx + i;
            gx = u[idx + 1] - u[idx - 1];
            gy = u[idx + nx] - u[idx - nx];
            gmag = sqrt(gx*gx + gy*gy);
            if (gmag > max_grad) {
                max_grad = gmag;
                max_i = i;
                max_j = j;
            }
        }
    }

    *tip_x = (double)max_i;
    *tip_y = (double)max_j;
    return (max_grad > 0.01) ? 1 : 0;
}
