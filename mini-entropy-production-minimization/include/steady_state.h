/**
 * @file steady_state.h
 * @brief Steady-state analysis and minimum entropy production principle.
 *
 * The Minimum Entropy Production Principle (MEPP), formulated by
 * Ilya Prigogine (1945), states that for a system with linear
 * constitutive relations and time-independent boundary conditions,
 * the steady non-equilibrium state minimizes the total entropy
 * production rate sigma.
 *
 * This header defines data structures and algorithms for:
 *   - Finding steady states that minimize entropy production
 *   - Stability analysis of non-equilibrium steady states
 *   - Lyapunov function construction for dissipative systems
 *   - Bifurcation analysis from equilibrium to dissipative structures
 *
 * References:
 *   - Prigogine, I. (1945) Acad. Roy. Belg., Bull. Classe Sci., 31, 600
 *   - Glansdorff & Prigogine (1971) "Thermodynamic Theory of Structure,
 *     Stability and Fluctuations"
 *   - Nicolis & Prigogine (1977) "Self-Organization in Nonequilibrium Systems"
 *   - Kondepudi & Prigogine (2014) "Modern Thermodynamics", Ch. 17-19
 */

#ifndef STEADY_STATE_H
#define STEADY_STATE_H

#include "entropy_production.h"
#include "onsager_relations.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---------------------------------------------------------------------------
 * L1: Definitions -- Steady state descriptors
 * ------------------------------------------------------------------------- */

/**
 * @brief Non-equilibrium steady state descriptor.
 *
 * A steady state is characterized by constant fluxes and forces
 * despite ongoing entropy production. The system is held away from
 * equilibrium by external constraints (fixed boundary conditions,
 * continuous matter/energy flows through the system).
 *
 * In the linear regime, the steady-state forces X^* minimize the
 * entropy production sigma(X) = X^T L X subject to constraints.
 */
typedef struct {
    ThermodynamicState state;        /**< Thermodynamic state at steady state */
    double forces[EP_MAX_FORCES];    /**< Steady-state force values X_i^* */
    double fluxes[EP_MAX_FORCES];    /**< Steady-state flux values J_i^* */
    double sigma;                    /**< Entropy production rate at steady state [W/K] */
    int    num_forces;              /**< Number of active thermodynamic forces */
    int    is_mep;                  /**< 1 if verified as minimum entropy production state */
    int    is_stable;               /**< 1 if linear stability verified */
    double lyapunov_exponent;       /**< Largest Lyapunov exponent (negative = stable) */
} SteadyState;

/**
 * @brief Stability analysis result for a non-equilibrium steady state.
 *
 * Linear stability is determined by the eigenvalues of the Jacobian
 * matrix of the evolution equations evaluated at the steady state.
 * For systems obeying MEPP, the entropy production sigma serves as
 * a Lyapunov function, guaranteeing stability.
 *
 * The spectral abscissa (maximum real part of eigenvalues) determines
 * stability: if < 0, the state is asymptotically stable.
 */
typedef struct {
    double eigenvalues[EP_MAX_FORCES]; /**< Eigenvalues of the Jacobian */
    int    dim;                        /**< System dimension */
    double max_real_part;              /**< Maximum real part of eigenvalues */
    int    is_linearly_stable;         /**< 1 if all Re(lambda) <= 0 */
    int    is_asymptotically_stable;   /**< 1 if all Re(lambda) < 0 */
    double spectral_abscissa;          /**< max_i Re(lambda_i) */
    double damping_ratio;              /**< min |Re(lambda)| / max |Im(lambda)| */
} StabilityResult;

/**
 * @brief Constraint specification for MEP optimization.
 *
 * Boundary conditions fix some thermodynamic forces to constant values
 * (e.g., fixed temperature gradient, fixed concentration differences),
 * while other forces are free to adjust to minimize entropy production.
 */
typedef struct {
    int    constrained[EP_MAX_FORCES]; /**< constrained[i]=1 if force i is fixed */
    double fixed_values[EP_MAX_FORCES]; /**< Value of each fixed force */
    int    num_forces;                  /**< Total number of forces */
    int    num_constrained;             /**< Count of constrained forces */
} SteadyStateConstraints;

/**
 * @brief Evolution trajectory of a system approaching steady state.
 *
 * Records the time evolution of forces, fluxes, and entropy production
 * as the system relaxes from an initial non-steady state toward its
 * minimum-entropy-production steady state.
 *
 * The trajectory provides empirical verification of Prigogine's theorem:
 * d(sigma)/dt <= 0 at all times during the relaxation.
 */
typedef struct {
    double *time;                  /**< Time points [s] (length n_steps) */
    double *forces_traj;           /**< Forces at each time (n_steps x dim, row-major) */
    double *sigma_traj;            /**< Entropy production at each time [W/K] */
    double *dsigma_dt_traj;        /**< Time derivative of sigma at each time */
    int     n_steps;               /**< Number of time steps recorded */
    int     dim;                   /**< System dimension */
    double  relaxation_time;       /**< Estimated relaxation time to steady state [s] */
    int     is_monotonic;          /**< 1 if sigma decreases monotonically */
} EvolutionTrajectory;

/* ---------------------------------------------------------------------------
 * L2: Core Concepts -- MEP steady state solvers
 * ------------------------------------------------------------------------- */

/**
 * @brief Find the minimum entropy production steady state.
 *
 * Solves the constrained optimization problem:
 *   minimize  sigma(X) = X^T L X
 *   subject to X_k = X_k^fixed for all constrained k
 *
 * This is a convex quadratic program. The function solves the KKT
 * system using LDL^T decomposition of the Onsager matrix.
 *
 * The KKT conditions are:
 *   partial(sigma)/partial(X_u) = 2 (L_{uu} X_u + L_{uc} X_c) = 0
 * where subscript u = unconstrained, c = constrained.
 * Solution: X_u = -L_{uu}^{-1} L_{uc} X_c
 *
 * @param L              Onsager coefficient matrix (n x n)
 * @param constraints    Constraint specification
 * @param forces         [in/out] Initial guess, overwritten with solution
 * @param n              System dimension
 * @param tol            Convergence tolerance
 * @return               0 on success, -1 on failure
 *
 * Theorem: Prigogine (1945) -- the solution is the unique steady state
 *          when L is positive definite and constrained forces are consistent.
 * Complexity: O(n^3) for LDL^T decomposition + O(n^2) for solve
 */
int ss_find_mep_state(const OnsagerMatrix *L,
                      const SteadyStateConstraints *constraints,
                      double forces[], int n, double tol);

/**
 * @brief Compute fluxes at the MEP steady state.
 *
 * Given the steady-state forces X^*, compute:
 *   J_i = sum_{j=1}^n L_{ij} X^*_j
 *
 * @param L       Onsager matrix
 * @param forces  Steady-state force vector X^* (length n)
 * @param fluxes  [out] Steady-state flux vector J^* (length n)
 * @param n       System dimension
 * @return        0 on success, -1 on null pointer
 */
int ss_compute_steady_fluxes(const OnsagerMatrix *L,
                             const double forces[], double fluxes[], int n);

/**
 * @brief Compute the minimum entropy production rate.
 *
 * sigma_min = (X^*)^T L X^* = sum_{i=1}^n sum_{j=1}^n L_{ij} X^*_i X^*_j
 *
 * This is the minimum possible entropy production rate given the
 * boundary constraints -- the value achieved at the steady state.
 *
 * @param L       Onsager matrix
 * @param forces  Steady-state force vector X^* (length n)
 * @param n       System dimension
 * @return        Minimum entropy production rate sigma_min [W/K]
 */
double ss_minimum_sigma(const OnsagerMatrix *L, const double forces[], int n);

/* ---------------------------------------------------------------------------
 * L3: Mathematical Structures -- Evolution equations (gradient flow)
 * ------------------------------------------------------------------------- */

/**
 * @brief Compute time derivatives of forces from gradient flow dynamics.
 *
 * For systems with linear constitutive relations, force evolution
 * follows from the balance equations as a gradient flow:
 *   dX_i/dt = -gamma * partial(sigma) / partial(X_i)
 *           = -gamma * 2 * sum_j L_{ij} X_j
 *
 * where gamma > 0 is a relaxation rate. Forces not fixed by boundary
 * conditions evolve to minimize sigma. Fixed forces have dX/dt = 0.
 *
 * This gradient structure is the dynamical foundation of the MEPP:
 * the system naturally evolves to minimize entropy production.
 *
 * @param forces      Current force values (length n)
 * @param dforces_dt  [out] Time derivatives dX_i/dt (length n)
 * @param L           Onsager matrix (n x n, row-major)
 * @param gamma       Relaxation rate coefficient [1/s]
 * @param fixed_forces Fixed boundary forces (NaN for free forces)
 * @param n           System dimension
 *
 * Reference: Glansdorff & Prigogine (1971), Ch. 5
 */
void ss_force_evolution(const double forces[], double dforces_dt[],
                        const OnsagerMatrix *L, double gamma,
                        const double fixed_forces[], int n);

/**
 * @brief Integrate the force evolution equations using RK4.
 *
 * Advances the system from time t to t+dt using the classic
 * fourth-order Runge-Kutta method applied to the gradient flow
 * dynamics defined by ss_force_evolution().
 *
 * RK4 provides 4th-order accuracy in dt, making it suitable for
 * accurate trajectory computation with moderate step sizes.
 *
 * @param forces      [in/out] Force values, updated from X(t) to X(t+dt)
 * @param L           Onsager matrix
 * @param gamma       Relaxation rate [1/s]
 * @param fixed_forces Fixed boundary forces (NaN for free)
 * @param n           System dimension
 * @param dt          Time step [s]
 */
void ss_rk4_step(double forces[], const OnsagerMatrix *L, double gamma,
                 const double fixed_forces[], int n, double dt);

/**
 * @brief Simulate the full relaxation to steady state.
 *
 * Integrates the gradient flow dynamics from initial conditions
 * until convergence, recording the complete trajectory. The
 * integration stops when either:
 *   - |d(sigma)/dt| < tol (steady state reached)
 *   - max_steps is reached
 *
 * The trajectory can be used to verify Prigogine's theorem by
 * checking that sigma(t) decreases monotonically.
 *
 * @param trajectory  [out] Evolution trajectory (pre-allocated)
 * @param L           Onsager matrix
 * @param gamma       Relaxation rate [1/s]
 * @param fixed_forces Fixed boundary forces
 * @param forces_init Initial force values at t=0
 * @param n           System dimension
 * @param dt          Time step [s]
 * @param max_steps   Maximum number of integration steps
 * @param tol         Convergence tolerance on |d(sigma)/dt|
 * @return            Number of steps taken, or -1 on error
 */
int ss_simulate_relaxation(EvolutionTrajectory *trajectory,
                           const OnsagerMatrix *L, double gamma,
                           const double fixed_forces[],
                           const double forces_init[], int n,
                           double dt, int max_steps, double tol);

/* ---------------------------------------------------------------------------
 * L4: Fundamental Laws -- Stability theorems for steady states
 * ------------------------------------------------------------------------- */

/**
 * @brief Compute the linear stability of a steady state.
 *
 * Constructs the Jacobian matrix:
 *   J_{ij} = partial(dX_i/dt) / partial(X_j) |_{X = X^*}
 *
 * For the gradient flow dX_i/dt = -2 gamma sum_k L_{ik} X_k,
 * the Jacobian is J = -2 gamma L for unconstrained forces.
 *
 * Since L is positive definite, all eigenvalues of J are negative,
 * guaranteeing asymptotic stability of the MEP steady state.
 *
 * @param result  [out] Stability analysis result
 * @param forces  Steady-state force values X^*
 * @param L       Onsager matrix
 * @param gamma   Relaxation rate [1/s]
 * @param n       System dimension
 * @return        0 on success, -1 on error
 */
int ss_linear_stability(StabilityResult *result, const double forces[],
                        const OnsagerMatrix *L, double gamma, int n);

/**
 * @brief Compute the excess entropy production (Lyapunov function).
 *
 * The excess entropy production delta^2 sigma is a Lyapunov function
 * for the stability of non-equilibrium steady states:
 *
 *   delta^2 sigma = sum_i (J_i - J_i^*) (X_i - X_i^*)
 *                 = (delta X)^T L (delta X) >= 0
 *
 * This quantity is:
 *   - Positive definite: delta^2 sigma > 0 for any deviation
 *   - Its time derivative is non-positive: d(delta^2 sigma)/dt <= 0
 *
 * satisfying Lyapunov's conditions for asymptotic stability.
 *
 * @param forces         Current force values X
 * @param forces_steady  Steady-state force values X^*
 * @param L              Onsager matrix (n x n)
 * @param n              System dimension
 * @return               Excess entropy production >= 0
 *
 * Reference: Glansdorff & Prigogine (1971), Ch. 6, Sec. 4
 */
double ss_excess_entropy_production(const double forces[],
                                    const double forces_steady[],
                                    const OnsagerMatrix *L, int n);

/**
 * @brief Evaluate the Glansdorff-Prigogine universal evolution criterion.
 *
 * The criterion states that for any non-equilibrium system:
 *   d_X sigma = sum_i J_i * dX_i/dt <= 0
 *
 * where d_X sigma is the change in entropy production due to changes
 * in forces alone (keeping the Onsager matrix fixed).
 *
 * At the steady state, d_X sigma = 0. Away from steady state,
 * d_X sigma < 0, indicating evolution toward the MEP state.
 *
 * @param forces         Current force values
 * @param dforces_dt     Force time derivatives
 * @param L              Onsager matrix
 * @param n              System dimension
 * @return               d_X sigma = sum_i J_i dX_i/dt [W^2/(K^2*s)]
 *
 * Reference: Glansdorff & Prigogine (1971), Ch. 6, Eq. (6.47)
 */
double ss_glansdorff_prigogine_criterion(const double forces[],
                                          const double dforces_dt[],
                                          const OnsagerMatrix *L, int n);

/* ---------------------------------------------------------------------------
 * L5: Algorithms -- Constrained optimization methods
 * ------------------------------------------------------------------------- */

/**
 * @brief Solve the MEP problem using the Lagrange multiplier method.
 *
 * Minimize sigma = X^T L X subject to A X = b.
 * Constructs and solves the saddle-point (KKT) system:
 *   [ 2L   A^T ] [ X     ]   [ 0 ]
 *   [ A    0   ] [ lambda ] = [ b ]
 *
 * This is a direct method suitable for small to medium systems.
 * For singular or near-singular systems, use ss_projected_gradient().
 *
 * @param L         Onsager matrix (n x n, row-major)
 * @param A         Constraint matrix (m x n, row-major)
 * @param b         Constraint right-hand side (length m)
 * @param X         [out] Solution forces (length n)
 * @param n         Number of forces
 * @param m         Number of linear constraints
 * @return          0 on success, -1 on singular system
 *
 * Complexity: O((n+m)^3)
 */
int ss_lagrange_mep(const double *L, const double *A, const double *b,
                    double *X, int n, int m);

/**
 * @brief Solve the MEP problem using projected gradient descent.
 *
 * For large systems or ill-conditioned Onsager matrices, iteratively
 * minimizes sigma by gradient descent with projection onto the
 * constraint manifold.
 *
 * Algorithm:
 *   1. Compute gradient: g = 2 L X
 *   2. Step: X_new = X - alpha * g
 *   3. Project: X_new[constrained_i] = fixed_values[i]
 *   4. Repeat until |g_unconstrained| < tol
 *
 * @param forces      [in/out] Force values (updated iteratively)
 * @param L           Onsager matrix
 * @param constrained Constraint flags array
 * @param fixed_vals  Fixed force values
 * @param n           System dimension
 * @param step_size   Gradient step size alpha
 * @param max_iter    Maximum iterations
 * @param tol         Convergence tolerance on gradient norm
 * @return            Number of iterations taken, or -1 on failure
 *
 * Complexity: O(n^2) per iteration
 */
int ss_projected_gradient(double forces[], const OnsagerMatrix *L,
                          const int constrained[], const double fixed_vals[],
                          int n, double step_size, int max_iter, double tol);

/* ---------------------------------------------------------------------------
 * Utility functions
 * ------------------------------------------------------------------------- */

void ss_init_constraints(SteadyStateConstraints *c, int num_forces);
int  ss_add_constraint(SteadyStateConstraints *c, int force_idx, double value);
void ss_init_steady_state(SteadyState *ss, int num_forces);
void ss_free_steady_state(SteadyState *ss);
void ss_init_trajectory(EvolutionTrajectory *traj, int dim, int max_steps);
void ss_free_trajectory(EvolutionTrajectory *traj);

#ifdef __cplusplus
}
#endif

#endif /* STEADY_STATE_H */