/**
 * turing_solver.c — Numerical Solvers for Reaction-Diffusion PDEs
 *
 * Implements multiple time-integration schemes for the coupled
 * reaction-diffusion system:
 *   ∂u/∂t = f(u,v) + D_u ∇²u
 *   ∂v/∂t = g(u,v) + D_v ∇²v
 *
 * Methods:
 *   - Explicit Euler (O(dt), simple but conditionally stable)
 *   - RK2 Midpoint (O(dt²), good balance of accuracy/cost)
 *   - Classical RK4 (O(dt⁴), workhorse for non-stiff systems)
 *   - Semi-Implicit Euler (implicit diffusion, explicit reaction)
 *   - IMEX Euler (Implicit-Explicit, handles stiff diffusion)
 *   - ETD-RK4 (Exponential Time Differencing, exact linear part)
 *
 * Reference: LeVeque, R.J. (2007) Finite Difference Methods for ODEs and PDEs
 *            Press et al. (2007) Numerical Recipes, 3rd ed.
 *            Cox & Matthews (2002) J. Comput. Phys. 176(2):430-455
 *            Hundsdorfer & Verwer (2003) Num. Sol. of Time-Dependent ADR Eqs.
 *
 * Knowledge Coverage:
 *   L5: Time-stepping algorithms, adaptive step-size control, CFL condition
 *   L3: Finite difference discretization, operator splitting
 *   L4: Von Neumann stability analysis, convergence criteria
 */

#include "turing_solver.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==========================================================================
 * Solver Lifecycle
 * ==========================================================================*/

SolverState *solver_alloc(ModelType model, int Nx, int Ny,
                          double Lx, double Ly,
                          double Du, double Dv,
                          const ModelParams *params,
                          const SolverConfig *config,
                          IntegrationMethod method) {
    if (!params || !config) return NULL;
    if (Nx < 4 || Ny < 4) return NULL;
    if (Du < 0.0 || Dv < 0.0) return NULL;

    SolverState *s = (SolverState *)calloc(1, sizeof(SolverState));
    if (!s) return NULL;

    s->model = model;
    s->params = *params;
    s->config = *config;
    s->method = method;

    /* Allocate main state */
    s->rd = rd_alloc(Nx, Ny, Lx, Ly, Du, Dv, config->dt);
    if (!s->rd) { free(s); return NULL; }

    s->rd_prev = rd_alloc(Nx, Ny, Lx, Ly, Du, Dv, config->dt);
    if (!s->rd_prev) { rd_free(s->rd); free(s); return NULL; }

    /* Set reaction/Jacobian function pointers */
    /* (resolved per-call via dispatchers for flexibility) */

    /* Allocate workspace fields */
    int ghost = 1;
    s->ku1 = field_alloc(Nx, Ny, Lx, Ly, ghost);
    s->kv1 = field_alloc(Nx, Ny, Lx, Ly, ghost);
    s->ku2 = field_alloc(Nx, Ny, Lx, Ly, ghost);
    s->kv2 = field_alloc(Nx, Ny, Lx, Ly, ghost);
    s->ku3 = field_alloc(Nx, Ny, Lx, Ly, ghost);
    s->kv3 = field_alloc(Nx, Ny, Lx, Ly, ghost);
    s->ku4 = field_alloc(Nx, Ny, Lx, Ly, ghost);
    s->kv4 = field_alloc(Nx, Ny, Lx, Ly, ghost);
    s->lap_u = field_alloc(Nx, Ny, Lx, Ly, ghost);
    s->lap_v = field_alloc(Nx, Ny, Lx, Ly, ghost);
    s->u_temp = field_alloc(Nx, Ny, Lx, Ly, ghost);
    s->v_temp = field_alloc(Nx, Ny, Lx, Ly, ghost);

    /* Check allocations */
    if (!s->ku1 || !s->kv1 || !s->ku2 || !s->kv2 || !s->ku3 || !s->kv3 ||
        !s->ku4 || !s->kv4 || !s->lap_u || !s->lap_v || !s->u_temp || !s->v_temp) {
        solver_free(s);
        return NULL;
    }

    s->dt_current = config->dt;
    s->total_steps = 0;
    s->converged = 0;

    return s;
}

void solver_free(SolverState *s) {
    if (!s) return;
    rd_free(s->rd);
    rd_free(s->rd_prev);
    field_free(s->ku1); field_free(s->kv1);
    field_free(s->ku2); field_free(s->kv2);
    field_free(s->ku3); field_free(s->kv3);
    field_free(s->ku4); field_free(s->kv4);
    field_free(s->lap_u); field_free(s->lap_v);
    field_free(s->u_temp); field_free(s->v_temp);
    free(s);
}

/* ==========================================================================
 * Initial Conditions
 * ==========================================================================*/

int solver_set_initial(SolverState *s, const ScalarField2D *u_init,
                       const ScalarField2D *v_init) {
    if (!s || !u_init || !v_init) return -1;
    if (u_init->Nx != s->rd->u.Nx || u_init->Ny != s->rd->u.Ny) return -1;

    field_copy(u_init, &s->rd->u);
    field_copy(v_init, &s->rd->v);
    field_apply_bc(&s->rd->u, s->config.bc_x, s->config.bc_y);
    field_apply_bc(&s->rd->v, s->config.bc_x, s->config.bc_y);

    s->rd->t = 0.0;
    s->rd->step = 0;
    s->converged = 0;

    return 0;
}

int solver_set_random_initial(SolverState *s, const HomogeneousSteadyState *hss,
                              double noise_amplitude, uint64_t seed) {
    if (!s || !hss) return -1;

    field_init_random(&s->rd->u, hss->u_star, noise_amplitude, seed);
    field_init_random(&s->rd->v, hss->v_star, noise_amplitude, seed + 1);

    field_apply_bc(&s->rd->u, s->config.bc_x, s->config.bc_y);
    field_apply_bc(&s->rd->v, s->config.bc_x, s->config.bc_y);

    s->rd->t = 0.0;
    s->rd->step = 0;
    s->converged = 0;

    return 0;
}

/* ==========================================================================
 * RHS Evaluation
 * ==========================================================================*/

int solver_eval_rhs(SolverState *s, ScalarField2D *dudt, ScalarField2D *dvdt) {
    if (!s || !dudt || !dvdt) return -1;

    int Nx = s->rd->u.Nx, Ny = s->rd->u.Ny;

    /* Compute Laplacians */
    laplacian_compute(&s->rd->u, s->lap_u, s->config.stencil);
    laplacian_compute(&s->rd->v, s->lap_v, s->config.stencil);

    /* RHS = reaction + diffusion */
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            double u = field_get(&s->rd->u, i, j);
            double v = field_get(&s->rd->v, i, j);
            double lap_u = field_get(s->lap_u, i, j);
            double lap_v = field_get(s->lap_v, i, j);

            double f, g;
            model_reaction(s->model, u, v, &s->params, &f, &g);

            field_set(dudt, i, j, f + s->rd->Du * lap_u);
            field_set(dvdt, i, j, g + s->rd->Dv * lap_v);
        }
    }
    return 0;
}

int solver_eval_reaction(SolverState *s, ScalarField2D *r_u, ScalarField2D *r_v) {
    if (!s || !r_u || !r_v) return -1;
    int Nx = s->rd->u.Nx, Ny = s->rd->u.Ny;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            double u = field_get(&s->rd->u, i, j);
            double v = field_get(&s->rd->v, i, j);
            double f, g;
            model_reaction(s->model, u, v, &s->params, &f, &g);
            field_set(r_u, i, j, f);
            field_set(r_v, i, j, g);
        }
    }
    return 0;
}

int solver_eval_diffusion(SolverState *s, ScalarField2D *d_u, ScalarField2D *d_v) {
    if (!s || !d_u || !d_v) return -1;
    laplacian_compute(&s->rd->u, s->lap_u, s->config.stencil);
    laplacian_compute(&s->rd->v, s->lap_v, s->config.stencil);
    int Nx = s->rd->u.Nx, Ny = s->rd->u.Ny;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            field_set(d_u, i, j, s->rd->Du * field_get(s->lap_u, i, j));
            field_set(d_v, i, j, s->rd->Dv * field_get(s->lap_v, i, j));
        }
    }
    return 0;
}

/* ==========================================================================
 * Time-Stepping Kernels
 * ==========================================================================*/

int step_euler_explicit(SolverState *s) {
    if (!s) return -1;

    /* Save current state */
    field_copy(&s->rd->u, &s->rd_prev->u);
    field_copy(&s->rd->v, &s->rd_prev->v);

    /* Compute RHS */
    solver_eval_rhs(s, s->ku1, s->kv1);

    /* Forward Euler: u^{n+1} = u^n + dt * RHS */
    double dt = s->dt_current;
    int Nx = s->rd->u.Nx, Ny = s->rd->u.Ny;

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            double u_new = field_get(&s->rd->u, i, j) + dt * field_get(s->ku1, i, j);
            double v_new = field_get(&s->rd->v, i, j) + dt * field_get(s->kv1, i, j);
            field_set(&s->rd->u, i, j, u_new);
            field_set(&s->rd->v, i, j, v_new);
        }
    }

    field_apply_bc(&s->rd->u, s->config.bc_x, s->config.bc_y);
    field_apply_bc(&s->rd->v, s->config.bc_x, s->config.bc_y);

    return 0;
}

int step_rk2_midpoint(SolverState *s) {
    if (!s) return -1;

    double dt = s->dt_current;
    int Nx = s->rd->u.Nx, Ny = s->rd->u.Ny;

    /* Save current state */
    field_copy(&s->rd->u, &s->rd_prev->u);
    field_copy(&s->rd->v, &s->rd_prev->v);

    /* k1 = RHS(u^n) */
    solver_eval_rhs(s, s->ku1, s->kv1);

    /* u_half = u^n + dt/2 * k1 */
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            double uh = field_get(&s->rd_prev->u, i, j)
                      + 0.5 * dt * field_get(s->ku1, i, j);
            double vh = field_get(&s->rd_prev->v, i, j)
                      + 0.5 * dt * field_get(s->kv1, i, j);
            field_set(&s->rd->u, i, j, uh);
            field_set(&s->rd->v, i, j, vh);
        }
    }
    field_apply_bc(&s->rd->u, s->config.bc_x, s->config.bc_y);
    field_apply_bc(&s->rd->v, s->config.bc_x, s->config.bc_y);

    /* k2 = RHS(u_half) */
    solver_eval_rhs(s, s->ku2, s->kv2);

    /* u^{n+1} = u^n + dt * k2 */
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            double un = field_get(&s->rd_prev->u, i, j)
                      + dt * field_get(s->ku2, i, j);
            double vn = field_get(&s->rd_prev->v, i, j)
                      + dt * field_get(s->kv2, i, j);
            field_set(&s->rd->u, i, j, un);
            field_set(&s->rd->v, i, j, vn);
        }
    }
    field_apply_bc(&s->rd->u, s->config.bc_x, s->config.bc_y);
    field_apply_bc(&s->rd->v, s->config.bc_x, s->config.bc_y);

    return 0;
}

int step_rk4(SolverState *s) {
    if (!s) return -1;

    double dt = s->dt_current;
    int Nx = s->rd->u.Nx, Ny = s->rd->u.Ny;

    /* Save original state */
    field_copy(&s->rd->u, &s->rd_prev->u);
    field_copy(&s->rd->v, &s->rd_prev->v);

    /* --- k1 = RHS(u^n) --- */
    solver_eval_rhs(s, s->ku1, s->kv1);

    /* --- k2: u = u^n + dt/2 * k1 --- */
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            field_set(&s->rd->u, i, j,
                field_get(&s->rd_prev->u, i, j) + 0.5*dt*field_get(s->ku1,i,j));
            field_set(&s->rd->v, i, j,
                field_get(&s->rd_prev->v, i, j) + 0.5*dt*field_get(s->kv1,i,j));
        }
    }
    field_apply_bc(&s->rd->u, s->config.bc_x, s->config.bc_y);
    field_apply_bc(&s->rd->v, s->config.bc_x, s->config.bc_y);
    solver_eval_rhs(s, s->ku2, s->kv2);

    /* --- k3: u = u^n + dt/2 * k2 --- */
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            field_set(&s->rd->u, i, j,
                field_get(&s->rd_prev->u, i, j) + 0.5*dt*field_get(s->ku2,i,j));
            field_set(&s->rd->v, i, j,
                field_get(&s->rd_prev->v, i, j) + 0.5*dt*field_get(s->kv2,i,j));
        }
    }
    field_apply_bc(&s->rd->u, s->config.bc_x, s->config.bc_y);
    field_apply_bc(&s->rd->v, s->config.bc_x, s->config.bc_y);
    solver_eval_rhs(s, s->ku3, s->kv3);

    /* --- k4: u = u^n + dt * k3 --- */
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            field_set(&s->rd->u, i, j,
                field_get(&s->rd_prev->u, i, j) + dt*field_get(s->ku3,i,j));
            field_set(&s->rd->v, i, j,
                field_get(&s->rd_prev->v, i, j) + dt*field_get(s->kv3,i,j));
        }
    }
    field_apply_bc(&s->rd->u, s->config.bc_x, s->config.bc_y);
    field_apply_bc(&s->rd->v, s->config.bc_x, s->config.bc_y);
    solver_eval_rhs(s, s->ku4, s->kv4);

    /* --- u^{n+1} = u^n + dt/6 * (k1 + 2k2 + 2k3 + k4) --- */
    double dt6 = dt / 6.0;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            double un = field_get(&s->rd_prev->u, i, j)
                + dt6 * (field_get(s->ku1,i,j) + 2.0*field_get(s->ku2,i,j)
                       + 2.0*field_get(s->ku3,i,j) + field_get(s->ku4,i,j));
            double vn = field_get(&s->rd_prev->v, i, j)
                + dt6 * (field_get(s->kv1,i,j) + 2.0*field_get(s->kv2,i,j)
                       + 2.0*field_get(s->kv3,i,j) + field_get(s->kv4,i,j));
            field_set(&s->rd->u, i, j, un);
            field_set(&s->rd->v, i, j, vn);
        }
    }
    field_apply_bc(&s->rd->u, s->config.bc_x, s->config.bc_y);
    field_apply_bc(&s->rd->v, s->config.bc_x, s->config.bc_y);

    return 0;
}

int step_semi_implicit(SolverState *s) {
    if (!s) return -1;

    double dt = s->dt_current;
    double Du = s->rd->Du, Dv = s->rd->Dv;
    double dx = s->rd->u.dx, dy = s->rd->u.dy;
    int Nx = s->rd->u.Nx, Ny = s->rd->u.Ny;

    /* Save state */
    field_copy(&s->rd->u, &s->rd_prev->u);
    field_copy(&s->rd->v, &s->rd_prev->v);

    /* Step 1: Explicit reaction half-step to temporary */
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            double u = field_get(&s->rd_prev->u, i, j);
            double v = field_get(&s->rd_prev->v, i, j);
            double f, g;
            model_reaction(s->model, u, v, &s->params, &f, &g);
            field_set(s->u_temp, i, j, u + 0.5 * dt * f);
            field_set(s->v_temp, i, j, v + 0.5 * dt * g);
        }
    }
    field_apply_bc(s->u_temp, s->config.bc_x, s->config.bc_y);
    field_apply_bc(s->v_temp, s->config.bc_x, s->config.bc_y);

    /* Step 2: Implicit diffusion via Jacobi iteration
     * Solve (I - D*dt*∇²) u_{new} = u_temp */
    double alpha_u = Du * dt / (dx * dx);
    double alpha_v = Dv * dt / (dx * dx);

    /* Jacobi for u */
    for (int iter = 0; iter < 20; iter++) {
        field_copy(s->u_temp, &s->rd->u);
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                double rhs = field_get(s->u_temp, i, j);
                double un = (rhs + alpha_u *
                    (field_get(&s->rd->u, i+1, j) + field_get(&s->rd->u, i-1, j)
                   + field_get(&s->rd->u, i, j+1) + field_get(&s->rd->u, i, j-1)))
                    / (1.0 + 4.0 * alpha_u);
                field_set(s->u_temp, i, j, un);
            }
        }
        field_apply_bc(s->u_temp, s->config.bc_x, s->config.bc_y);
    }
    field_copy(s->u_temp, &s->rd->u);

    /* Jacobi for v */
    for (int iter = 0; iter < 20; iter++) {
        field_copy(s->v_temp, &s->rd->v);
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                double rhs = field_get(s->v_temp, i, j);
                double vn = (rhs + alpha_v *
                    (field_get(&s->rd->v, i+1, j) + field_get(&s->rd->v, i-1, j)
                   + field_get(&s->rd->v, i, j+1) + field_get(&s->rd->v, i, j-1)))
                    / (1.0 + 4.0 * alpha_v);
                field_set(s->v_temp, i, j, vn);
            }
        }
        field_apply_bc(s->v_temp, s->config.bc_x, s->config.bc_y);
    }
    field_copy(s->v_temp, &s->rd->v);

    /* Step 3: Explicit reaction half-step to finalize */
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            double u = field_get(&s->rd->u, i, j);
            double v = field_get(&s->rd->v, i, j);
            double f, g;
            model_reaction(s->model, u, v, &s->params, &f, &g);
            field_set(&s->rd->u, i, j, u + 0.5 * dt * f);
            field_set(&s->rd->v, i, j, v + 0.5 * dt * g);
        }
    }
    field_apply_bc(&s->rd->u, s->config.bc_x, s->config.bc_y);
    field_apply_bc(&s->rd->v, s->config.bc_x, s->config.bc_y);

    return 0;
}

int step_imex_euler(SolverState *s) {
    if (!s) return -1;

    double dt = s->dt_current;
    double Du = s->rd->Du, Dv = s->rd->Dv;
    double dx = s->rd->u.dx;
    int Nx = s->rd->u.Nx, Ny = s->rd->u.Ny;

    field_copy(&s->rd->u, &s->rd_prev->u);
    field_copy(&s->rd->v, &s->rd_prev->v);

    /* Compute reaction term explicitly */
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            double u = field_get(&s->rd_prev->u, i, j);
            double v = field_get(&s->rd_prev->v, i, j);
            double f, g;
            model_reaction(s->model, u, v, &s->params, &f, &g);
            field_set(s->u_temp, i, j, u + dt * f);
            field_set(s->v_temp, i, j, v + dt * g);
        }
    }

    /* Solve implicit diffusion: (I - D*dt*∇²) u_{n+1} = u_temp */
    jacobi_solve_helmholtz(&s->rd->u, s->u_temp,
                           Du * dt / (dx * dx), dx, dx,
                           s->config.bc_x, s->config.bc_y, 30, 1e-8);
    jacobi_solve_helmholtz(&s->rd->v, s->v_temp,
                           Dv * dt / (dx * dx), dx, dx,
                           s->config.bc_x, s->config.bc_y, 30, 1e-8);

    return 0;
}

/* ==========================================================================
 * ETD-RK4 (Exponential Time Differencing)
 *
 * For systems of the form: du/dt = Lu + N(u)
 * where L is the stiff linear diffusion operator and N is the nonlinear
 * reaction term. ETD methods treat L exactly via the matrix exponential.
 *
 * ETD-RK4 coefficients (Cox & Matthews 2002):
 *   a_n = u_n + α(dt*L)(u_n + dt*N(u_n))
 *   ... (4 stages using φ functions)
 *
 * For the spectral (Fourier) case, L is diagonal: L_k = -D*k²
 * The exponential e^{L*dt} is just e^{-D*k²*dt} in Fourier space.
 *
 * Reference: Cox & Matthews (2002) J. Comput. Phys. 176(2):430-455
 *            Kassam & Trefethen (2005) SIAM J. Sci. Comput. 26(4):1214-1233
 * ==========================================================================*/

/**
 * φ₀(z) = e^z
 * φ₁(z) = (e^z - 1)/z
 * φ₂(z) = (e^z - 1 - z)/z²
 * φ₃(z) = (e^z - 1 - z - z²/2)/z³
 *
 * For small |z|, use Taylor expansions to avoid catastrophic cancellation.
 */
static double phi0(double z) {
    if (z < -700.0) return 0.0;
    if (z > 700.0) return exp(700.0);
    return exp(z);
}

static double phi1(double z) {
    if (fabs(z) < 0.01) {
        /* Taylor: φ₁(z) = 1 + z/2 + z²/6 + z³/24 + ... */
        return 1.0 + z/2.0 + z*z/6.0 + z*z*z/24.0;
    }
    return (exp(z) - 1.0) / z;
}

static double phi2(double z) {
    if (fabs(z) < 0.01) {
        return 0.5 + z/6.0 + z*z/24.0 + z*z*z/120.0;
    }
    return (exp(z) - 1.0 - z) / (z * z);
}

static double phi3(double z) {
    if (fabs(z) < 0.01) {
        return 1.0/6.0 + z/24.0 + z*z/120.0 + z*z*z/720.0;
    }
    double ez = exp(z);
    return (ez - 1.0 - z - 0.5*z*z) / (z * z * z);
}

int step_etd_rk4(SolverState *s) {
    if (!s) return -1;

    double dt = s->dt_current;
    double Du = s->rd->Du, Dv = s->rd->Dv;
    double dx = s->rd->u.dx, dy = s->rd->u.dy;
    int Nx = s->rd->u.Nx, Ny = s->rd->u.Ny;

    /* For simplicity, use the "slab" approximation:
     * Treat L as the 5-point Laplacian, diagonalized in Fourier space.
     * k²_max ≈ π²(1/dx² + 1/dy²) for the highest resolved mode. */
    double kx_max = M_PI / dx;
    double ky_max = M_PI / dy;
    double k2_u = -Du * (kx_max*kx_max + ky_max*ky_max);
    double k2_v = -Dv * (kx_max*kx_max + ky_max*ky_max);

    /* Save state */
    field_copy(&s->rd->u, &s->rd_prev->u);
    field_copy(&s->rd->v, &s->rd_prev->v);

    /* ETD-RK4 using integrating factor approach:
     *   N(u) = f(u,v) (reaction part, handled explicitly)
     *   Each stage: compute N, then apply e^{L*dt} or φ functions.
     *
     * Simplified: spectral ETD-RK4 using φ functions applied per point.
     */

    /* For non-spectral (finite difference) implementation, use the fact that
     * the linear part is Du*∇²u, and approximate e^{L*dt} by
     * the implicit Euler operator (I - dt*Du*∇²)^{-1} */

    /* Stage 1: compute N(u^n) */
    solver_eval_reaction(s, s->ku1, s->kv1);

    /* a = e^{L*dt/2} * u^n + φ₁(L*dt/2) * dt/2 * N(u^n) */
    /* For simplicity, use implicit Euler to approximate e^{L*dt/2} */
    double hdt = 0.5 * dt;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            double u_n = field_get(&s->rd_prev->u, i, j);
            double v_n = field_get(&s->rd_prev->v, i, j);
            double u_a = u_n + hdt * field_get(s->ku1, i, j);
            double v_a = v_n + hdt * field_get(s->kv1, i, j);
            field_set(s->u_temp, i, j, u_a);
            field_set(s->v_temp, i, j, v_a);
        }
    }
    /* Apply implicit diffusion half-step */
    jacobi_solve_helmholtz(s->u_temp, s->u_temp,
                           Du * hdt / (dx*dx), dx, dy,
                           s->config.bc_x, s->config.bc_y, 10, 1e-6);
    jacobi_solve_helmholtz(s->v_temp, s->v_temp,
                           Dv * hdt / (dx*dx), dx, dy,
                           s->config.bc_x, s->config.bc_y, 10, 1e-6);

    /* Stage 2: N(a) at half-step */
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            double u = field_get(s->u_temp, i, j);
            double v = field_get(s->v_temp, i, j);
            double f, g;
            model_reaction(s->model, u, v, &s->params, &f, &g);
            field_set(s->ku2, i, j, f);
            field_set(s->kv2, i, j, g);
        }
    }

    /* b = e^{L*dt/2} * u^n + φ₁(L*dt/2) * dt/2 * N(a) */
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            double u_n = field_get(&s->rd_prev->u, i, j);
            double v_n = field_get(&s->rd_prev->v, i, j);
            field_set(s->u_temp, i, j, u_n + hdt * field_get(s->ku2, i, j));
            field_set(s->v_temp, i, j, v_n + hdt * field_get(s->kv2, i, j));
        }
    }
    jacobi_solve_helmholtz(s->u_temp, s->u_temp,
                           Du * hdt / (dx*dx), dx, dy,
                           s->config.bc_x, s->config.bc_y, 10, 1e-6);
    jacobi_solve_helmholtz(s->v_temp, s->v_temp,
                           Dv * hdt / (dx*dx), dx, dy,
                           s->config.bc_x, s->config.bc_y, 10, 1e-6);

    /* Stage 3: N(b) */
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            double u = field_get(s->u_temp, i, j);
            double v = field_get(s->v_temp, i, j);
            double f, g;
            model_reaction(s->model, u, v, &s->params, &f, &g);
            field_set(s->ku3, i, j, f);
            field_set(s->kv3, i, j, g);
        }
    }

    /* c = e^{L*dt} * u^n + φ₁(L*dt) * dt * (2*N(b) - N(u^n)) */
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            double u_n = field_get(&s->rd_prev->u, i, j);
            double v_n = field_get(&s->rd_prev->v, i, j);
            double n_u = 2.0*field_get(s->ku3,i,j) - field_get(s->ku1,i,j);
            double n_v = 2.0*field_get(s->kv3,i,j) - field_get(s->kv1,i,j);
            field_set(s->u_temp, i, j, u_n + dt * n_u);
            field_set(s->v_temp, i, j, v_n + dt * n_v);
        }
    }
    jacobi_solve_helmholtz(s->u_temp, s->u_temp,
                           Du * dt / (dx*dx), dx, dy,
                           s->config.bc_x, s->config.bc_y, 10, 1e-6);
    jacobi_solve_helmholtz(s->v_temp, s->v_temp,
                           Dv * dt / (dx*dx), dx, dy,
                           s->config.bc_x, s->config.bc_y, 10, 1e-6);

    /* Stage 4: N(c) */
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            double u = field_get(s->u_temp, i, j);
            double v = field_get(s->v_temp, i, j);
            double f, g;
            model_reaction(s->model, u, v, &s->params, &f, &g);
            field_set(s->ku4, i, j, f);
            field_set(s->kv4, i, j, g);
        }
    }

    /* Final assembly (simplified) */
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            double u_n = field_get(&s->rd_prev->u, i, j);
            double v_n = field_get(&s->rd_prev->v, i, j);
            double n_avg = (field_get(s->ku1,i,j) + 2.0*field_get(s->ku2,i,j)
                          + 2.0*field_get(s->ku3,i,j) + field_get(s->ku4,i,j))/6.0;
            double m_avg = (field_get(s->kv1,i,j) + 2.0*field_get(s->kv2,i,j)
                          + 2.0*field_get(s->kv3,i,j) + field_get(s->kv4,i,j))/6.0;
            field_set(&s->rd->u, i, j, u_n + dt * n_avg);
            field_set(&s->rd->v, i, j, v_n + dt * m_avg);
        }
    }
    /* Apply implicit diffusion to final state */
    jacobi_solve_helmholtz(&s->rd->u, &s->rd->u,
                           Du * dt / (dx*dx), dx, dy,
                           s->config.bc_x, s->config.bc_y, 10, 1e-6);
    jacobi_solve_helmholtz(&s->rd->v, &s->rd->v,
                           Dv * dt / (dx*dx), dx, dy,
                           s->config.bc_x, s->config.bc_y, 10, 1e-6);

    return 0;
}

/* ==========================================================================
 * Adaptive Time-Stepping
 * ==========================================================================*/

double solver_estimate_error(SolverState *s) {
    if (!s) return -1.0;

    /* Richard extrapolation: compare full step (dt) vs two half-steps (dt/2)
     * Error ≈ ||u_full - u_two_half|| / (2^p - 1) where p = order */

    double dt_save = s->dt_current;
    double dt_half = 0.5 * dt_save;

    /* Save state */
    field_copy(&s->rd->u, s->u_temp);
    field_copy(&s->rd->v, s->v_temp);

    /* Two half-steps */
    s->dt_current = dt_half;
    solver_step(s);
    solver_step(s);

    /* Save result */
    field_copy(&s->rd->u, s->ku1);
    field_copy(&s->rd->v, s->kv1);

    /* Restore and do one full step */
    field_copy(s->u_temp, &s->rd->u);
    field_copy(s->v_temp, &s->rd->v);
    s->dt_current = dt_save;
    solver_step(s);

    /* Compute error */
    double error = field_diff_norm(s->ku1, &s->rd->u);

    /* p = 4 for RK4, p = 2 for midpoint, p = 1 for Euler */
    int order = (s->method == METHOD_RK4) ? 4 :
                (s->method == METHOD_RK2_MIDPOINT) ? 2 : 1;

    error /= (pow(2.0, order) - 1.0);
    s->dt_current = dt_save;
    s->error_estimate = error;

    return error;
}

double solver_adjust_dt(SolverState *s, double error, double dt_current) {
    if (!s) return dt_current;
    if (error <= 0.0) return dt_current;

    double safety = s->config.safety_factor;
    double tol = s->config.tolerance;

    /* Standard PID controller for step size (Press et al. 2007, §17.2) */
    double factor = safety * pow(tol / (error + 1e-15), 0.25);  /* 1/p for p=4 */

    /* Limit changes */
    if (factor > 1.5) factor = 1.5;
    if (factor < 0.5) factor = 0.5;

    double new_dt = dt_current * factor;

    /* Clamp to allowed range */
    if (new_dt > s->config.dt_max) new_dt = s->config.dt_max;
    if (new_dt < s->config.dt_min) new_dt = s->config.dt_min;

    /* CFL check */
    double cfl = solver_cfl_limit(s);
    if (new_dt > cfl) new_dt = cfl;

    return new_dt;
}

double solver_cfl_limit(const SolverState *s) {
    if (!s) return 1e10;

    double dx = s->rd->u.dx;
    double dy = s->rd->u.dy;
    double D_max = (s->rd->Du > s->rd->Dv) ? s->rd->Du : s->rd->Dv;

    /* CFL condition for explicit diffusion: D*dt/(dx²) ≤ 0.5 */
    double dx2 = dx * dx;
    double dy2 = dy * dy;
    double h2_inv = 1.0 / dx2 + 1.0 / dy2;

    return 0.5 / (D_max * h2_inv);
}

/* ==========================================================================
 * Unified Stepping + Run Control
 * ==========================================================================*/

int solver_step(SolverState *s) {
    if (!s) return -1;

    int ret;
    switch (s->method) {
        case METHOD_EULER_EXPLICIT: ret = step_euler_explicit(s); break;
        case METHOD_RK2_MIDPOINT:   ret = step_rk2_midpoint(s);   break;
        case METHOD_RK4:            ret = step_rk4(s);             break;
        case METHOD_SEMI_IMPLICIT:  ret = step_semi_implicit(s);  break;
        case METHOD_IMEX_EULER:     ret = step_imex_euler(s);     break;
        case METHOD_ETD_RK4:        ret = step_etd_rk4(s);        break;
        default:                    return -1;
    }

    if (ret == 0) {
        s->rd->t += s->dt_current;
        s->rd->step++;
        s->total_steps++;

        /* Check for numerical blow-up */
        if (solver_detect_blowup(s, 1e6)) {
            return -1; /* numerical instability */
        }

        /* Adaptive time-stepping */
        if (s->config.adaptive_dt) {
            double error = solver_estimate_error(s);
            s->dt_current = solver_adjust_dt(s, error, s->dt_current);
        }

        /* Steady state check */
        if (solver_is_steady(s)) {
            s->converged = 1;
            return 1;
        }
    }

    return ret;
}

int solver_run(SolverState *s, int max_steps) {
    if (!s) return -1;

    int steps = 0;
    for (int i = 0; i < max_steps; i++) {
        int ret = solver_step(s);
        steps++;
        if (ret != 0) break; /* error or steady state */
    }
    return steps;
}

int solver_run_until(SolverState *s, double t_target) {
    if (!s) return -1;

    int steps = 0;
    int max_steps = s->config.max_steps;

    for (int i = 0; i < max_steps; i++) {
        if (s->rd->t + s->dt_current > t_target) {
            s->dt_current = t_target - s->rd->t;
        }
        int ret = solver_step(s);
        steps++;
        if (ret != 0) break;
        if (s->rd->t >= t_target - 1e-12) break;
    }
    return steps;
}

int solver_run_to_steady(SolverState *s) {
    if (!s) return -1;

    int steps = 0;
    int max_steps = s->config.max_steps;

    for (int i = 0; i < max_steps; i++) {
        int ret = solver_step(s);
        steps++;
        if (ret != 0) break;
    }
    return steps;
}

/* ==========================================================================
 * Diagnostics
 * ==========================================================================*/

double solver_max_derivative(const SolverState *s) {
    if (!s) return -1.0;

    solver_eval_rhs((SolverState *)s, s->ku1, s->kv1);
    double max_val = 0.0;
    int Nx = s->rd->u.Nx, Ny = s->rd->u.Ny;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            double du = fabs(field_get(s->ku1, i, j));
            double dv = fabs(field_get(s->kv1, i, j));
            if (du > max_val) max_val = du;
            if (dv > max_val) max_val = dv;
        }
    }
    return max_val;
}

double solver_rms_derivative(const SolverState *s) {
    if (!s) return -1.0;

    solver_eval_rhs((SolverState *)s, s->ku1, s->kv1);
    double sum_sq = 0.0;
    int Nx = s->rd->u.Nx, Ny = s->rd->u.Ny;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            double du = field_get(s->ku1, i, j);
            double dv = field_get(s->kv1, i, j);
            sum_sq += du*du + dv*dv;
        }
    }
    return sqrt(sum_sq / (2.0 * Nx * Ny));
}

int solver_is_steady(const SolverState *s) {
    if (!s) return 0;
    double rms = solver_rms_derivative(s);
    return (rms < s->config.tolerance) ? 1 : 0;
}

int solver_detect_blowup(const SolverState *s, double threshold) {
    if (!s) return 1;
    int Nx = s->rd->u.Nx, Ny = s->rd->u.Ny;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            double u = field_get(&s->rd->u, i, j);
            double v = field_get(&s->rd->v, i, j);
            if (!is_finite_val(u) || !is_finite_val(v)) return 1;
            if (fabs(u) > threshold || fabs(v) > threshold) return 1;
        }
    }
    return 0;
}

const char *solver_status_string(const SolverState *s) {
    if (!s) return "NULL";
    if (solver_detect_blowup(s, 1e6)) return "UNSTABLE";
    if (s->converged) return "STEADY";
    return "RUNNING";
}

/* ==========================================================================
 * Jacobi Solver for Helmholtz Equation
 *
 * Solves (I - α∇²) u = rhs using Jacobi iteration.
 *
 * Iteration: u_new[i][j] = (rhs[i][j] + α/Δx² * (u[i+1][j] + u[i-1][j]
 *                          + u[i][j+1] + u[i][j-1])) / (1 + 4α/Δx²)
 *
 * Convergence guaranteed for α > 0 (strictly diagonally dominant).
 * Spectral radius: ρ = 4α/(Δx²*(1+4α/Δx²)) * cos²(π/(2N)) < 1.
 * ==========================================================================*/

int jacobi_solve_helmholtz(ScalarField2D *u, const ScalarField2D *rhs,
                           double alpha, double dx, double dy,
                           BoundaryCondition bc_x, BoundaryCondition bc_y,
                           int max_iter, double tol) {
    if (!u || !rhs) return -1;
    if (u->Nx != rhs->Nx || u->Ny != rhs->Ny) return -1;

    int Nx = u->Nx, Ny = u->Ny;
    double idx2 = 1.0 / (dx * dx);
    double idy2 = 1.0 / (dy * dy);
    double coeff = 1.0 / (1.0 + 2.0 * alpha * (idx2 + idy2));

    /* Workspace */
    ScalarField2D *u_old = field_alloc(Nx, Ny, u->Lx, u->Ly, u->ghost);
    if (!u_old) return -1;

    for (int iter = 0; iter < max_iter; iter++) {
        field_copy(u, u_old);
        field_apply_bc(u_old, bc_x, bc_y);

        double max_change = 0.0;

        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                double rhs_val = field_get(rhs, i, j);
                double u_xx = (field_get(u_old, i+1, j) + field_get(u_old, i-1, j)) * idx2;
                double u_yy = (field_get(u_old, i, j+1) + field_get(u_old, i, j-1)) * idy2;
                double u_new = coeff * (rhs_val + alpha * (u_xx + u_yy));
                double change = fabs(u_new - field_get(u_old, i, j));
                if (change > max_change) max_change = change;
                field_set(u, i, j, u_new);
            }
        }
        field_apply_bc(u, bc_x, bc_y);

        if (max_change < tol) {
            field_free(u_old);
            return iter + 1;
        }
    }
    field_free(u_old);
    return max_iter;  /* not fully converged */
}
