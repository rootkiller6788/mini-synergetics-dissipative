/**
 * @file bz_reaction_diffusion.h
 * @brief Reaction-diffusion coupling for BZ waves: Turing instability,
 *        spiral waves, target patterns.
 *
 * Knowledge coverage:
 *   L2: Reaction-diffusion coupling, Turing mechanism, wave propagation
 *   L3: PDE discretization (finite difference), Laplacian operator
 *   L4: Turing instability condition (D_u != D_v, activator-inhibitor)
 *   L5: Finite-difference PDE solver, spiral wave simulation
 *   L6: BZ spiral wave formation on a 2D domain
 *
 * Reaction-diffusion equation:
 *   du/dt = f(u,v) + D_u * nabla^2 u
 *   dv/dt = g(u,v) + D_v * nabla^2 v
 *
 * References:
 *   - Turing (1952) Phil. Trans. R. Soc. B 237, 37-72
 *   - Winfree (1972) Science 175, 634-636
 *   - Murray (2002) "Mathematical Biology" Vol II, Ch 2
 */

#ifndef BZ_REACTION_DIFFUSION_H
#define BZ_REACTION_DIFFUSION_H

#include "bz_reaction.h"

typedef enum {
    BZ_RD_DIM_1D = 1,
    BZ_RD_DIM_2D = 2
} bz_rd_dimension_t;

typedef enum {
    BZ_BC_NEUMANN    = 0,
    BZ_BC_DIRICHLET  = 1,
    BZ_BC_PERIODIC   = 2
} bz_rd_boundary_t;

typedef struct {
    int     nx;
    double  L;
    double  dx;
    double *u;
    double *v;
    double  Du;
    double  Dv;
    bz_rd_boundary_t bc;
    double  u_bc;
    double  v_bc;
} bz_rd_1d_t;

typedef struct {
    int     nx;
    int     ny;
    double  Lx;
    double  Ly;
    double  dx;
    double  dy;
    double *u;
    double *v;
    double  Du;
    double  Dv;
    bz_rd_boundary_t bc;
} bz_rd_2d_t;

typedef enum {
    BZ_RD_KINETICS_OREGONATOR = 0,
    BZ_RD_KINETICS_FITZHUGH   = 1,
    BZ_RD_KINETICS_BARKLEY    = 2
} bz_rd_kinetics_t;

typedef struct {
    double a;
    double b;
    double epsilon;
} bz_barkley_params_t;

void bz_laplacian_1d(double *lap, const double *u, int nx, double dx,
                      bz_rd_boundary_t bc);

void bz_laplacian_2d(double *lap, const double *u, int nx, int ny,
                      double dx, bz_rd_boundary_t bc);

int bz_turing_unstable(double fu, double fv, double gu, double gv,
                        double Du, double Dv);

double bz_turing_critical_wavenumber(double fu, double fv,
                                      double gu, double gv,
                                      double Du, double Dv);

void bz_rd_1d_step(double *u_new, double *v_new,
                    const double *u, const double *v,
                    const bz_rd_1d_t *domain, double dt,
                    const bz_oregonator_params_t *params);

void bz_rd_2d_step(double *u_new, double *v_new,
                    const double *u, const double *v,
                    const bz_rd_2d_t *domain,
                    const bz_barkley_params_t *bp, double dt);

void bz_spiral_seed_2d(double *u, double *v, int nx, int ny,
                        int cx, int cy, int length);

void bz_rd_init_random(double *u, double *v, int n, double low, double high);

int bz_find_spiral_tip(double *tip_x, double *tip_y,
                        const double *u, const double *v,
                        int nx, int ny);

#endif /* BZ_REACTION_DIFFUSION_H */
