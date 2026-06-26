#ifndef HAKEN_DYNAMICS_H
#define HAKEN_DYNAMICS_H
#include "haken_core.h"

/* ==============================================================
 * haken_dynamics.h — Nonlinear Dynamics & Numerical Integration
 *
 * Implements time evolution of synergetic systems:
 *   dq/dt = L(α)q + N(q, α)
 *
 * with proper handling of timescale separation for stiff systems.
 * Includes projection-based splitting methods that respect the
 * slow/fast decomposition.
 *
 * References:
 *   Haken (1983) Advanced Synergetics, §1.3-1.5
 *   Press et al. (2007) Numerical Recipes, §17.5 (stiff ODEs)
 *   Gear (1971) Numerical Initial Value Problems
 * ============================================================== */

/* ==============================================================
 * L1: Integration method types
 * ============================================================== */

typedef enum {
    HAKEN_INTEGRATOR_EULER       = 0,
    HAKEN_INTEGRATOR_RK2         = 1,
    HAKEN_INTEGRATOR_RK4         = 2,
    HAKEN_INTEGRATOR_IMPLICIT    = 3, /* Backward Euler for stiff modes */
    HAKEN_INTEGRATOR_SPLIT       = 4, /* Split: explicit slow, implicit fast */
    HAKEN_INTEGRATOR_ADAPTIVE    = 5  /* Adaptive RK4(5) with error control */
} Haken_IntegratorType;

/** Integration configuration */
typedef struct {
    Haken_IntegratorType type;
    double  dt;                 /* Time step (or initial step for adaptive) */
    double  t_start;
    double  t_end;
    double  atol;               /* Absolute tolerance for adaptive methods */
    double  rtol;               /* Relative tolerance */
    int     max_steps;          /* Maximum number of steps */
    bool    record_trajectory;  /* Save full trajectory */
    double* trajectory;         /* Output: (n_steps+1)×n_dim */
    int     trajectory_capacity;/* Allocated trajectory capacity */
    int     n_saved;            /* Number of saved states */
} Haken_IntegratorConfig;

/* ==============================================================
 * L3: Nonlinear coupling structures
 * ============================================================== */

/** Quadratic coupling tensor: N_i = Σ_{j,k} Q_{ijk} q_j q_k */
typedef struct {
    int      n_dim;
    double*  Q;              /* Q[i + n*(j + n*k)], flattened rank-3 tensor */
    bool     symmetric_jk;   /* Q_{ijk} = Q_{ikj} */
} Haken_QuadraticCoupling;

/** Cubic coupling tensor: N_i = Σ_{j,k,l} C_{ijkl} q_j q_k q_l */
typedef struct {
    int      n_dim;
    double*  C;              /* Flattened rank-4 tensor */
    bool     symmetric_jkl;  /* Full symmetry under permutation of j,k,l */
} Haken_CubicCoupling;

/** Generic polynomial nonlinearity: N(q) = Σ_{d=2}^D A^{(d)} q^{⊗d} */
typedef struct {
    int      n_dim;
    int      max_degree;      /* Maximum polynomial degree */
    double** coefficients;    /* coefficients[d] for degree d+2 */
    int*     n_terms;         /* Number of monomial terms per degree */
} Haken_PolynomialCoupling;

/* ==============================================================
 * L5: Numerical integration API
 * ============================================================== */

/** Integrate the full system dynamics for one time step.
 *  @param sys   System with current state
 *  @param dt    Time step
 *  @param type  Integrator type */
void haken_integrate_step(Haken_System* sys, double dt,
                           Haken_IntegratorType type);

/** Integrate the system from t_start to t_end.
 *  @param sys    System (state is updated in place)
 *  @param config Integration configuration
 *  @return       Number of steps taken, -1 on failure */
int haken_integrate(Haken_System* sys, Haken_IntegratorConfig* config);

/** Euler forward step: q(t+dt) = q(t) + dt * f(q(t)) */
void haken_step_euler(Haken_System* sys, double dt);

/** RK4 step (classical Runge-Kutta 4th order):
 *  k1 = f(q), k2 = f(q+dt*k1/2), k3 = f(q+dt*k2/2),
 *  k4 = f(q+dt*k3), q_next = q + dt*(k1+2k2+2k3+k4)/6 */
void haken_step_rk4(Haken_System* sys, double dt);

/** RK2 (midpoint) step. */
void haken_step_rk2(Haken_System* sys, double dt);

/** Backward Euler for stiff systems:
 *  q_{n+1} = q_n + dt * f(q_{n+1})
 *  Solved via Newton iteration. */
void haken_step_backward_euler(Haken_System* sys, double dt,
                                double tol, int max_iter);

/** Split-step integrator: explicit for slow modes,
 *  implicit for fast modes. Uses projection onto critical/stable
 *  subspaces to separate the dynamics.
 *
 *  Algorithm (Strang splitting):
 *    1. Half-step implicit on fast: q_f ← q_f + (dt/2)·f_f(q_s,q_f)
 *    2. Full-step explicit on slow: q_s ← q_s + dt·f_s(q_s,q_f)
 *    3. Half-step implicit on fast: q_f ← q_f + (dt/2)·f_f(q_s,q_f)
 *
 *  This respects the slaving principle numerically. */
void haken_step_split(Haken_System* sys, double dt,
                       double tol, int max_iter);

/** Adaptive RK4(5) step using the Dormand-Prince pair.
 *  Takes a step with error estimation and adjusts dt.
 *  @param sys    System (state updated)
 *  @param dt     Input: suggested step, Output: recommended next step
 *  @param err    Output: estimated local truncation error
 *  @return       0 if step accepted, 1 if rejected (needs retry) */
int haken_step_rk45_adaptive(Haken_System* sys, double* dt, double* err,
                               double atol, double rtol);

/** Evaluate the right-hand side of the ODE:
 *  f(q) = L(α)·q + N(q)
 *  @param sys  System descriptor
 *  @param q    Input state
 *  @param dq   Output: time derivative */
void haken_rhs(const Haken_System* sys, const double* q, double* dq);

/** Evaluate only the linear part: f_lin = L(α)·q */
void haken_rhs_linear(const Haken_System* sys, const double* q,
                       double* f_lin);

/** Evaluate only the nonlinear part: f_nl = N(q) */
void haken_rhs_nonlinear(const Haken_System* sys, const double* q,
                          double* f_nl);

/* ==============================================================
 * Quadratic coupling constructor API
 * ============================================================== */

/** Allocate a quadratic coupling tensor */
Haken_QuadraticCoupling* haken_quadratic_create(int n, bool symmetric);

/** Free quadratic coupling tensor */
void haken_quadratic_free(Haken_QuadraticCoupling* qc);

/** Set all coupling coefficients to zero */
void haken_quadratic_zero(Haken_QuadraticCoupling* qc);

/** Set a specific coupling coefficient Q_{ijk} */
void haken_quadratic_set(Haken_QuadraticCoupling* qc,
                          int i, int j, int k, double value);

/** Evaluate quadratic coupling: N_i = Σ_{j≤k} Q_{ijk} q_j q_k */
void haken_quadratic_eval(const Haken_QuadraticCoupling* qc,
                           const double* q, double* Nq);

/** Allocate a cubic coupling tensor */
Haken_CubicCoupling* haken_cubic_create(int n, bool symmetric);

/** Free cubic coupling */
void haken_cubic_free(Haken_CubicCoupling* cc);

/** Evaluate cubic coupling: N_i = Σ_{j,k,l} C_{ijkl} q_j q_k q_l */
void haken_cubic_eval(const Haken_CubicCoupling* cc,
                       const double* q, double* Nq);

/** Allocate a polynomial coupling */
Haken_PolynomialCoupling* haken_polynomial_create(int n, int max_degree);

/** Free polynomial coupling */
void haken_polynomial_free(Haken_PolynomialCoupling* pc);

/** Evaluate polynomial coupling up to max_degree */
void haken_polynomial_eval(const Haken_PolynomialCoupling* pc,
                            const double* q, double* Nq);

/* ==============================================================
 * L6: Canonical model constructors
 * ============================================================== */

/** Create the Haken-Kelvin model (Haken 1977, §5.1):
 *  The simplest system exhibiting the slaving principle:
 *    dq1/dt = α·q1 - a·q1·q2
 *    dq2/dt = -β·q2 + b·q1²
 *  where α ≈ 0 (critical), β ≫ 0 (stable).
 *  Adiabatic elimination: q2 ≈ (b/β)·q1²
 *  → dq1/dt = α·q1 - (a·b/β)·q1³  (Ginzburg-Landau form) */
Haken_System* haken_create_prototype_model(double alpha, double beta,
                                             double a, double b);

/** Create the Lorenz system (reduced to Haken form):
 *  Near the first instability, the Lorenz equations reduce to
 *  a single-order-parameter system. Parameters: r (Rayleigh),
 *  σ (Prandtl), b_geo (geometry factor).
 *
 *  Reference: Haken (1975) Phys. Lett. 53A, 77 */
Haken_System* haken_create_lorenz_instability(double r, double sigma,
                                                double b_geo);

#endif /* HAKEN_DYNAMICS_H */
