/**
 * turing_common.h — Core Types and Utilities for Turing Pattern Formation
 *
 * Reference: Turing, A.M. (1952) "The Chemical Basis of Morphogenesis"
 *            Phil. Trans. R. Soc. Lond. B 237(641):37-72
 *
 * This header defines the fundamental data structures used throughout
 * the Turing pattern formation numerical framework:
 *   - 2D scalar fields on uniform grids
 *   - Boundary condition types
 *   - Reaction-diffusion PDE parameter sets
 *   - Pattern characterization metrics
 *   - Initial condition generators
 *
 * Knowledge Coverage: L1 (Definitions), L2 (Core Concepts), L3 (Math Structures)
 */

#ifndef TURING_COMMON_H
#define TURING_COMMON_H

#include <stddef.h>
#include <stdint.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---------------------------------------------------------------------------
 * L1: Boundary Condition Types
 * Dirichlet: fixed value at boundary (e.g., u=0)
 * Neumann: zero flux across boundary (∂u/∂n=0)
 * Periodic: wrap-around (torus topology)
 * Robin: mixed condition αu + β∂u/∂n = γ
 * -------------------------------------------------------------------------*/
typedef enum {
    BC_DIRICHLET = 0,
    BC_NEUMANN   = 1,
    BC_PERIODIC  = 2,
    BC_ROBIN     = 3
} BoundaryCondition;

/* ---------------------------------------------------------------------------
 * L1: Laplacian Discretization Stencil Types
 * 5-point: standard 2nd-order accurate (∂²/∂x² + ∂²/∂y²)
 * 7-point: includes rotated axis points for isotropy
 * 9-point: full isotropic with 4th-order correction (Mehrstellen)
 * -------------------------------------------------------------------------*/
typedef enum {
    LAPLACIAN_5PT  = 0,
    LAPLACIAN_7PT  = 1,
    LAPLACIAN_9PT  = 2
} LaplacianStencil;

/* ---------------------------------------------------------------------------
 * L1: Pattern Type Classification
 * Homogeneous: no spatial structure, uniform concentration
 * Spots: localized high-concentration regions on low background
 * Stripes/Labyrinths: elongated connected structures
 * Hexagons: regular hexagonal lattice (H⁻ or H⁺)
 * Solitons: localized traveling structures
 * Chaos: spatiotemporal irregularity
 * -------------------------------------------------------------------------*/
typedef enum {
    PATTERN_HOMOGENEOUS = 0,
    PATTERN_SPOTS       = 1,
    PATTERN_STRIPES     = 2,
    PATTERN_LABYRINTHS  = 3,
    PATTERN_HEXAGONS    = 4,
    PATTERN_SOLITONS    = 5,
    PATTERN_CHAOS       = 6,
    PATTERN_UNKNOWN     = 7
} PatternType;

/* ---------------------------------------------------------------------------
 * L3: 2D Scalar Field — fundamental data container
 *
 * Mathematical structure: u ∈ ℝ^(Nx × Ny)
 *
 * Memory layout: row-major, u[i][j] = data[i * Ny + j]
 * Ghost cells: optional padding for boundary condition handling
 * -------------------------------------------------------------------------*/
typedef struct {
    double  *data;       /* field values, size Nx*Ny */
    int      Nx;         /* grid points in x-direction */
    int      Ny;         /* grid points in y-direction */
    double   Lx;         /* physical domain length in x */
    double   Ly;         /* physical domain length in y */
    double   dx;         /* grid spacing in x */
    double   dy;         /* grid spacing in y */
    int      ghost;      /* number of ghost cell layers */
} ScalarField2D;

/* ---------------------------------------------------------------------------
 * L3: Reaction-Diffusion System — two coupled PDEs
 *
 * System:
 *   ∂u/∂t = f(u, v; p) + D_u ∇²u
 *   ∂v/∂t = g(u, v; p) + D_v ∇²v
 *
 * where p = (p1, p2, ..., pn) are kinetic parameters.
 * -------------------------------------------------------------------------*/
typedef struct {
    ScalarField2D u;     /* activator / short-range species */
    ScalarField2D v;     /* inhibitor / long-range species */
    double        Du;    /* diffusion coefficient of u */
    double        Dv;    /* diffusion coefficient of v */
    double        dt;    /* time step */
    double        t;     /* current simulation time */
    int           step;  /* current step number */
} ReactionDiffusion2D;

/* ---------------------------------------------------------------------------
 * L3: Kinetic Model Parameters — polymorphic parameter set
 *
 * Each model uses a subset of these parameters.
 * Gray-Scott:     uses F, k
 * Gierer-Meinhardt: uses a,b,c, rho_a, rho_h, mu_a, mu_h
 * Brusselator:    uses A, B
 * Schnakenberg:   uses a, b
 * FitzHugh-Nagumo: uses ε, β, γ, I_ext
 * Lengyel-Epstein: uses a, b, c
 * -------------------------------------------------------------------------*/
typedef struct {
    /* universal parameters */
    double  F;          /* feed rate (Gray-Scott) */
    double  k;          /* kill rate (Gray-Scott) */

    /* Gierer-Meinhardt parameters */
    double  a;          /* activator production rate */
    double  b;          /* inhibitor production rate */
    double  c;          /* basal production */
    double  rho_a;      /* activator saturation */
    double  rho_h;      /* inhibitor saturation */
    double  mu_a;       /* activator decay */
    double  mu_h;       /* inhibitor decay */

    /* Brusselator parameters */
    double  A;          /* reservoir concentration A (Brusselator) */
    double  B;          /* reservoir concentration B (Brusselator) */

    /* FitzHugh-Nagumo parameters */
    double  epsilon;    /* timescale separation */
    double  beta;       /* threshold parameter */
    double  gamma;      /* recovery rate */
    double  I_ext;      /* external stimulus */

    /* Lengyel-Epstein parameters */
    /* a, b reused; c is the starch indicator concentration */

    /* Thomas model parameters */
    double  rho;        /* reaction rate (Thomas) */
    double  K;          /* inhibition constant (Thomas) */
    double  alpha;      /* ratio of diffusion/time scales (Thomas) */

    /* Threshold/switch parameters */
    double  threshold;  /* activation threshold */
    double  saturation; /* saturation level */
} ModelParams;

/* ---------------------------------------------------------------------------
 * L1: Model Type Enum
 * -------------------------------------------------------------------------*/
typedef enum {
    MODEL_GRAY_SCOTT       = 0,
    MODEL_GIERER_MEINHARDT = 1,
    MODEL_BRUSSELATOR      = 2,
    MODEL_SCHNAKENBERG     = 3,
    MODEL_FITZHUGH_NAGUMO  = 4,
    MODEL_LENGYEL_EPSTEIN  = 5,
    MODEL_THOMAS           = 6,
    MODEL_USER_DEFINED     = 7
} ModelType;

/* ---------------------------------------------------------------------------
 * L3: Reaction Function Signature — universal kinetics interface
 *
 * Each model provides its own f(u,v) and g(u,v).
 * Returns 0 on success, -1 on error (e.g., NaN detected).
 * -------------------------------------------------------------------------*/
typedef int (*ReactionFunc)(double u, double v, const ModelParams *p,
                            double *f_out, double *g_out);

/* ---------------------------------------------------------------------------
 * L3: Jacobian Function Signature — for stability analysis
 *
 * Computes the 2×2 Jacobian matrix of the reaction terms:
 *   J = [[∂f/∂u, ∂f/∂v],
 *        [∂g/∂u, ∂g/∂v]]
 *
 * Evaluated at a given (u,v).
 * -------------------------------------------------------------------------*/
typedef int (*JacobianFunc)(double u, double v, const ModelParams *p,
                            double J[2][2]);

/* ---------------------------------------------------------------------------
 * L2: Homogeneous Steady State (HSS)
 *
 * A spatially uniform solution (u*, v*) such that:
 *   f(u*, v*) = 0  and  g(u*, v*) = 0
 *
 * The HSS is stable to homogeneous perturbations but unstable
 * to certain spatial modes — this is the essence of Turing instability.
 * -------------------------------------------------------------------------*/
typedef struct {
    double  u_star;     /* activator steady-state concentration */
    double  v_star;     /* inhibitor steady-state concentration */
    int     converged;  /* 1 if solver converged, 0 otherwise */
    int     iterations; /* number of Newton iterations used */
    double  residual;   /* final residual ||f(u*,v*)|| */
} HomogeneousSteadyState;

/* ---------------------------------------------------------------------------
 * L2: Dispersion Relation Result
 *
 * λ(k) = max(eigenvalues of J - k²D)
 *
 * where J is the Jacobian at HSS, D = diag(Du, Dv).
 * If Re(λ(k)) > 0 for some k ≠ 0, the system is Turing-unstable.
 * -------------------------------------------------------------------------*/
typedef struct {
    int      n_modes;          /* number of Fourier modes computed */
    double  *wavenumbers;      /* k values, size n_modes */
    double  *growth_rates;     /* Re(λ(k)), size n_modes */
    double   k_critical;       /* k with maximum growth rate */
    double   lambda_max;       /* max Re(λ(k)) */
    double   k_min;            /* lower bound of unstable band */
    double   k_max;            /* upper bound of unstable band */
    int      turing_unstable;  /* 1 if Turing instability detected */
    double   predicted_wavelength; /* 2π/k_critical */
} DispersionRelation;

/* ---------------------------------------------------------------------------
 * L2: Pattern Metrics — quantitative characterization
 * -------------------------------------------------------------------------*/
typedef struct {
    double   mean;             /* spatial mean */
    double   variance;         /* spatial variance */
    double   skewness;         /* distribution skew */
    double   kurtosis;         /* distribution kurtosis */
    double   dominant_wavelength; /* from FFT peak */
    double   pattern_amplitude;   /* max - min */
    double   spot_density;     /* spots per unit area */
    double   stripe_alignment; /* orientation coherence (0-1) */
    PatternType classified_type;
} PatternMetrics;

/* ---------------------------------------------------------------------------
 * L5: Solver Configuration
 * -------------------------------------------------------------------------*/
typedef struct {
    double   dt;               /* time step */
    double   t_max;            /* final simulation time */
    int      max_steps;        /* maximum number of steps */
    int      save_interval;    /* steps between data output */
    double   tolerance;        /* steady-state tolerance */
    LaplacianStencil stencil;  /* discretization type */
    BoundaryCondition bc_x;    /* x-boundary condition */
    BoundaryCondition bc_y;    /* y-boundary condition */
    int      adaptive_dt;      /* 1 = enabled */
    double   dt_min;           /* minimum allowed time step */
    double   dt_max;           /* maximum allowed time step */
    double   safety_factor;    /* CFL safety factor for adaptive dt */
} SolverConfig;

/* ---------------------------------------------------------------------------
 * Core API: Field Lifecycle
 * -------------------------------------------------------------------------*/

/** Allocate a 2D scalar field with ghost cells for boundary handling. */
ScalarField2D *field_alloc(int Nx, int Ny, double Lx, double Ly, int ghost);

/** Free a 2D scalar field and all associated memory. */
void field_free(ScalarField2D *f);

/** Access element (i,j) of the field, accounting for ghost cells. */
static inline double field_get(const ScalarField2D *f, int i, int j) {
    int gi = i + f->ghost;
    int gj = j + f->ghost;
    return f->data[gi * (f->Ny + 2 * f->ghost) + gj];
}

/** Set element (i,j) of the field. */
static inline void field_set(ScalarField2D *f, int i, int j, double val) {
    int gi = i + f->ghost;
    int gj = j + f->ghost;
    f->data[gi * (f->Ny + 2 * f->ghost) + gj] = val;
}

/** Get raw pointer to interior data (without ghost cells) for fast access. */
double *field_interior(ScalarField2D *f);

/** Total number of interior grid points. */
static inline int field_size(const ScalarField2D *f) {
    return f->Nx * f->Ny;
}

/* ---------------------------------------------------------------------------
 * Core API: Reaction-Diffusion System Lifecycle
 * -------------------------------------------------------------------------*/

/** Allocate and initialize a reaction-diffusion system. */
ReactionDiffusion2D *rd_alloc(int Nx, int Ny, double Lx, double Ly,
                               double Du, double Dv, double dt);

/** Free the reaction-diffusion system. */
void rd_free(ReactionDiffusion2D *rd);

/* ---------------------------------------------------------------------------
 * Core API: Field Operations
 * -------------------------------------------------------------------------*/

/** Fill the entire field with a constant value (including ghost cells). */
void field_fill(ScalarField2D *f, double value);

/** Fill the interior of the field with a constant value. */
void field_fill_interior(ScalarField2D *f, double value);

/** Copy src field to dst (both must have identical dimensions). */
int field_copy(const ScalarField2D *src, ScalarField2D *dst);

/** Compute the L2 norm of the field difference: ||a - b||₂. */
double field_diff_norm(const ScalarField2D *a, const ScalarField2D *b);

/** Apply boundary conditions to the ghost cells of the field. */
void field_apply_bc(ScalarField2D *f, BoundaryCondition bc_x,
                    BoundaryCondition bc_y);

/* ---------------------------------------------------------------------------
 * Core API: Laplacian Computation
 *
 * Computes ∇²u using the specified stencil.
 * Returns 0 on success.
 * -------------------------------------------------------------------------*/

/**
 * Compute ∇² field on a 2D grid.
 * Stencil choices:
 *   5-point:   ∂²u/∂x² + ∂²u/∂y² with O(h²) accuracy
 *   7-point:   adds diagonal neighbors for improved isotropy
 *   9-point:   Mehrstellen 4th-order compact scheme
 *
 * Complexity: O(Nx·Ny) per evaluation
 * Reference: Thomas, J.W. (1995) Numerical Partial Differential Equations
 */
int laplacian_compute(const ScalarField2D *f, ScalarField2D *lap,
                      LaplacianStencil stencil);

/**
 * Compute the spectral Laplacian (via FFT) for dispersion analysis.
 * Assumes periodic boundary conditions.
 * Complexity: O(N log N)
 */
int laplacian_spectral(const ScalarField2D *f, ScalarField2D *lap);

/* ---------------------------------------------------------------------------
 * Core API: Initial Conditions
 * -------------------------------------------------------------------------*/

/**
 * Initialize field with random noise around a mean value.
 * amplitude controls the noise magnitude.
 */
void field_init_random(ScalarField2D *f, double mean, double amplitude,
                       uint64_t seed);

/** Initialize with a Gaussian bump centered at (cx, cy). */
void field_init_gaussian(ScalarField2D *f, double cx, double cy,
                         double sigma, double amplitude);

/** Initialize with a sinusoidal pattern (for testing modes). */
void field_init_sinusoidal(ScalarField2D *f, double kx, double ky,
                           double amplitude, double offset);

/** Initialize with a random dot pattern (for testing spots). */
void field_init_random_dots(ScalarField2D *f, int num_dots, double radius,
                            double amplitude, uint64_t seed);

/** Initialize with a gradient along x-axis. */
void field_init_gradient_x(ScalarField2D *f, double left_val, double right_val);

/* ---------------------------------------------------------------------------
 * Core API: Model Parameters
 * -------------------------------------------------------------------------*/

/** Set default model parameters for a given model type. */
void model_params_default(ModelType type, ModelParams *p);

/** Validate that model parameters are physically meaningful. */
int model_params_validate(ModelType type, const ModelParams *p);

/** Print model parameters to stdout for diagnostics. */
void model_params_print(ModelType type, const ModelParams *p);

/* ---------------------------------------------------------------------------
 * Core API: Dispersion Analysis
 * -------------------------------------------------------------------------*/

/** Compute the dispersion relation λ(k) for a given HSS and diffusion ratio. */
DispersionRelation *dispersion_alloc(int n_modes);

/** Free dispersion relation structure. */
void dispersion_free(DispersionRelation *dr);

/* ---------------------------------------------------------------------------
 * Utility: Safe math operations with domain checking.
 * -------------------------------------------------------------------------*/

/** Safe exponential: clamps argument to avoid overflow. */
double safe_exp(double x);

/** Safe division: returns 0.0 if denominator is near zero. */
double safe_div(double num, double den);

/** Check if a double value is a valid finite number. */
int is_finite_val(double x);

#endif /* TURING_COMMON_H */
