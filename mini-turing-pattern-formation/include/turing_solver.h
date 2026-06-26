/**
 * turing_solver.h — Numerical Solvers for Reaction-Diffusion PDEs
 *
 * Reference: LeVeque, R.J. (2007) Finite Difference Methods for ODEs and PDEs
 *            Strikwerda, J.C. (2004) Finite Difference Schemes and PDEs, 2nd ed.
 *            Hundsdorfer & Verwer (2003) Numerical Solution of Time-Dependent
 *            Advection-Diffusion-Reaction Equations
 *            Press, W.H. et al. (2007) Numerical Recipes, 3rd ed.
 *
 * Provides multiple numerical integration schemes for the reaction-diffusion
 * system ∂u/∂t = f(u,v) + D_u ∇²u, ∂v/∂t = g(u,v) + D_v ∇²v.
 *
 * Knowledge Coverage:
 *   L5: Explicit Euler, RK4, Semi-Implicit, Adaptive time-stepping
 *   L3: FDM discretization, CFL condition, von Neumann stability
 *   L4: Convergence properties, stability criteria
 */

#ifndef TURING_SOLVER_H
#define TURING_SOLVER_H

#include "turing_common.h"
#include "turing_models.h"

/* ---------------------------------------------------------------------------
 * L5: Integration Method Enum
 * -------------------------------------------------------------------------*/
typedef enum {
    METHOD_EULER_EXPLICIT = 0,
    METHOD_RK4            = 1,
    METHOD_RK2_MIDPOINT   = 2,
    METHOD_SEMI_IMPLICIT  = 3,
    METHOD_IMEX_EULER     = 4,   /* IMplicit-EXplicit: implicit diffusion,
                                    explicit reaction */
    METHOD_ETD_RK4        = 5    /* Exponential Time Differencing RK4 */
} IntegrationMethod;

/* ---------------------------------------------------------------------------
 * L5: Solver State — mutable simulation context
 * -------------------------------------------------------------------------*/
typedef struct {
    ModelType           model;
    ModelParams         params;
    SolverConfig        config;
    IntegrationMethod   method;
    ReactionDiffusion2D *rd;         /* current state */
    ReactionDiffusion2D *rd_prev;    /* previous state (for adaptive methods) */
    ReactionFunc        reaction_f;
    JacobianFunc        jacobian_f;

    /* Workspace buffers for multi-stage methods */
    ScalarField2D  *ku1, *kv1;       /* RK stage 1 */
    ScalarField2D  *ku2, *kv2;       /* RK stage 2 */
    ScalarField2D  *ku3, *kv3;       /* RK stage 3 */
    ScalarField2D  *ku4, *kv4;       /* RK stage 4 */
    ScalarField2D  *lap_u, *lap_v;   /* Laplacian workspace */
    ScalarField2D  *u_temp, *v_temp; /* temporary fields */

    /* Adaptive time-stepping */
    double  dt_current;              /* current time step */
    double  error_estimate;          /* local truncation error estimate */
    int     dt_rejections;           /* count of rejected steps */

    /* Statistics */
    int     total_steps;             /* total steps taken */
    int     rejected_steps;          /* total rejected steps */
    double  cpu_time;                /* cumulative wall/CPU time */
    int     converged;               /* 1 if steady state reached */
    double  final_delta;             /* ||u_new - u_old|| / dt at final step */
} SolverState;

/* ---------------------------------------------------------------------------
 * Solver Lifecycle
 * -------------------------------------------------------------------------*/

/**
 * Allocate and initialize a solver for the given model, domain, and config.
 * Returns NULL on invalid parameters.
 */
SolverState *solver_alloc(ModelType model, int Nx, int Ny,
                          double Lx, double Ly,
                          double Du, double Dv,
                          const ModelParams *params,
                          const SolverConfig *config,
                          IntegrationMethod method);

/** Free solver and all associated memory. */
void solver_free(SolverState *s);

/* ---------------------------------------------------------------------------
 * Initial Conditions
 * -------------------------------------------------------------------------*/

/** Set initial conditions for the solver. Returns 0 on success. */
int solver_set_initial(SolverState *s, const ScalarField2D *u_init,
                       const ScalarField2D *v_init);

/** Set random initial conditions around HSS with specified noise amplitude. */
int solver_set_random_initial(SolverState *s, const HomogeneousSteadyState *hss,
                              double noise_amplitude, uint64_t seed);

/* ---------------------------------------------------------------------------
 * Time Stepping
 * -------------------------------------------------------------------------*/

/**
 * Advance the solver by one time step.
 * Returns:
 *   0 — step completed successfully
 *   1 — steady state reached (||∂u/∂t|| < tolerance)
 *  -1 — numerical instability detected
 *  -2 — NaN/Inf detected
 */
int solver_step(SolverState *s);

/**
 * Run solver for n steps or until t_max or steady state.
 * Returns number of steps actually taken.
 */
int solver_run(SolverState *s, int max_steps);

/**
 * Run solver until a specified simulation time.
 * Returns number of steps taken.
 */
int solver_run_until(SolverState *s, double t_target);

/**
 * Run solver until steady state is detected.
 * Returns number of steps taken.
 */
int solver_run_to_steady(SolverState *s);

/* ---------------------------------------------------------------------------
 * Single-Step Kernels — exposed for analysis and testing
 *
 * Each kernel computes one full step of the corresponding method.
 * Derivatives dudt, dvdt are computed and stored in ku1, kv1.
 * -------------------------------------------------------------------------*/

/** Explicit Euler step: u^{n+1} = u^n + dt * (f + D*∇²u). O(dt) accuracy. */
int step_euler_explicit(SolverState *s);

/**
 * Classical 4th-order Runge-Kutta step. O(dt⁴) accuracy.
 *
 * k1 = f(u^n)
 * k2 = f(u^n + dt/2 * k1)
 * k3 = f(u^n + dt/2 * k2)
 * k4 = f(u^n + dt * k3)
 * u^{n+1} = u^n + dt/6 * (k1 + 2k2 + 2k3 + k4)
 */
int step_rk4(SolverState *s);

/** 2nd-order midpoint method: k1=f(u^n), k2=f(u^n+dt/2*k1), u^{n+1}=u^n+dt*k2. */
int step_rk2_midpoint(SolverState *s);

/**
 * Semi-implicit (IMEX) step: diffusion treated implicitly, reaction explicit.
 * (I - dt*D*∇²) u^{n+1} = u^n + dt * f(u^n)
 * Solves the resulting linear system via Jacobi iteration.
 */
int step_semi_implicit(SolverState *s);

/** IMEX Euler: implicit diffusion, explicit reaction. No splitting error. */
int step_imex_euler(SolverState *s);

/**
 * Exponential Time Differencing RK4 (ETD-RK4).
 * Handles stiff linear diffusion term exactly via matrix exponential.
 * Suitable for large diffusion coefficients.
 *
 * Reference: Cox, S.M. & Matthews, P.C. (2002) J. Comput. Phys. 176(2):430-455
 */
int step_etd_rk4(SolverState *s);

/* ---------------------------------------------------------------------------
 * Right-Hand Side (RHS) Evaluation
 *
 * Computes the full PDE right-hand side:
 *   dudt = f(u,v) + Du * ∇²u
 *   dvdt = g(u,v) + Dv * ∇²v
 * -------------------------------------------------------------------------*/

/**
 * Evaluate the reaction-diffusion RHS at the current state.
 * Stores results in dudt_field, dvdt_field.
 * Complexity: O(Nx·Ny) grid operations
 */
int solver_eval_rhs(SolverState *s, ScalarField2D *dudt, ScalarField2D *dvdt);

/**
 * Evaluate reaction term only (without diffusion) for operator splitting.
 */
int solver_eval_reaction(SolverState *s, ScalarField2D *r_u, ScalarField2D *r_v);

/**
 * Evaluate diffusion term only (without reaction).
 */
int solver_eval_diffusion(SolverState *s, ScalarField2D *d_u, ScalarField2D *d_v);

/* ---------------------------------------------------------------------------
 * Adaptive Time-Stepping
 *
 * Adjusts dt based on local truncation error estimate.
 * Uses step doubling (Richardson extrapolation) or embedded RK.
 * -------------------------------------------------------------------------*/

/** Estimate local error by comparing full and half steps. */
double solver_estimate_error(SolverState *s);

/** Adjust dt based on error estimate and tolerance. Returns new dt. */
double solver_adjust_dt(SolverState *s, double error, double dt_current);

/** Compute the CFL stability limit for explicit diffusion. */
double solver_cfl_limit(const SolverState *s);

/* ---------------------------------------------------------------------------
 * Diagnostics and Convergence
 * -------------------------------------------------------------------------*/

/** Compute the maximum absolute time derivative (∞-norm of ∂u/∂t, ∂v/∂t). */
double solver_max_derivative(const SolverState *s);

/** Compute the L2 norm of the time derivative over the domain. */
double solver_rms_derivative(const SolverState *s);

/** Check if steady state has been reached. */
int solver_is_steady(const SolverState *s);

/** Check for numerical blow-up (any value exceeding threshold). */
int solver_detect_blowup(const SolverState *s, double threshold);

/** Get solver status string for logging. */
const char *solver_status_string(const SolverState *s);

/* ---------------------------------------------------------------------------
 * Jacobi Iteration for Implicit Methods
 *
 * Solves (I - α∇²) u_new = u_rhs using Jacobi iteration.
 * Used internally by semi-implicit and IMEX methods.
 * -------------------------------------------------------------------------*/

/**
 * Jacobi iteration for the Helmholtz-like system (I - α∇²)u = rhs.
 * Parameters:
 *   alpha: diffusion coefficient * dt (= D*dt/dx² for implicit diffusion)
 *   max_iter: maximum iterations
 *   tol: convergence tolerance
 * Returns number of iterations used.
 */
int jacobi_solve_helmholtz(ScalarField2D *u, const ScalarField2D *rhs,
                           double alpha, double dx, double dy,
                           BoundaryCondition bc_x, BoundaryCondition bc_y,
                           int max_iter, double tol);

#endif /* TURING_SOLVER_H */
