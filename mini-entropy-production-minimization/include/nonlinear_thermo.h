/**
 * @file nonlinear_thermo.h
 * @brief Nonlinear irreversible thermodynamics and far-from-equilibrium systems.
 *
 * While the linear regime (Onsager relations, MEPP) describes systems
 * close to equilibrium, many interesting phenomena -- pattern formation,
 * chemical oscillations, turbulence, life itself -- occur far from
 * equilibrium where nonlinearities dominate.
 *
 * This header covers:
 *   - Nonlinear constitutive relations (flux as nonlinear function of forces)
 *   - Generalized entropy production beyond bilinear form
 *   - Maximum Entropy Production Principle (MEPP debate)
 *   - Dissipative adaptation and self-organization
 *   - Thermodynamic metrics far from equilibrium
 *
 * References:
 *   - Glansdorff & Prigogine (1971)
 *   - Martyushev & Seleznev (2006) "Maximum entropy production principle"
 *   - Kleidon & Lorenz (2005) "Non-equilibrium Thermodynamics and the
 *     Production of Entropy"
 *   - England (2015) "Dissipative adaptation in driven self-assembly"
 *   - Maes, C. (2018) "Non-Dissipative Effects in Nonequilibrium Systems"
 */

#ifndef NONLINEAR_THERMO_H
#define NONLINEAR_THERMO_H

#include "entropy_production.h"
#include "onsager_relations.h"

#ifdef __cplusplus
extern "C" {
#endif

#define NT_MAX_TERMS 32

/* ---------------------------------------------------------------------------
 * L1: Definitions -- Nonlinear extension structures
 * ------------------------------------------------------------------------- */

/**
 * @brief Nonlinear constitutive relation.
 *
 * Beyond the linear regime, fluxes can be expanded as:
 *   J_i = sum_j L_{ij} X_j + sum_{j,k} M_{ijk} X_j X_k + ...
 *
 * This structure holds the expansion coefficients up to a specified
 * order. The linear term corresponds to the Onsager matrix;
 * higher-order terms represent nonlinear corrections.
 */
typedef struct {
    double *linear;              /**< Linear coefficients L_{ij} (n x n) */
    double *quadratic;           /**< Quadratic coefficients M_{ijk} (n x n x n) */
    double *cubic;               /**< Cubic coefficients N_{ijkl} (n x n x n x n) */
    int     dim;                 /**< System dimension n */
    int     order;               /**< Maximum order (2 or 3) */
    int     num_quadratic_terms; /**< Number of nonzero quadratic terms */
    int     num_cubic_terms;     /**< Number of nonzero cubic terms */
} NonlinearOnsager;

/**
 * @brief Generalized entropy production rate.
 *
 * For nonlinear constitutive relations, the entropy production
 * is not simply the bilinear form sigma = sum J_i X_i. Additional
 * contributions arise from higher-order coupling:
 *
 *   sigma = sum_i J_i X_i
 *         = sum_{i,j} L_{ij} X_i X_j
 *         + sum_{i,j,k} M_{ijk} X_i X_j X_k
 *         + sum_{i,j,k,l} N_{ijkl} X_i X_j X_k X_l
 *         + ...
 *
 * For the Second Law to hold, sigma(X) >= 0 for all physically
 * admissible X, which imposes constraints on the coefficients.
 */
typedef struct {
    double sigma_linear;         /**< Linear contribution sigma^(2) */
    double sigma_quadratic;      /**< Quadratic contribution sigma^(3) */
    double sigma_cubic;          /**< Cubic contribution sigma^(4) */
    double sigma_total;          /**< Total sigma = sum of all orders */
    int    is_positive;          /**< 1 if sigma >= 0 */
} NonlinearEntropyProduction;

/**
 * @brief Far-from-equilibrium thermodynamic state descriptor.
 *
 * Extends the basic thermodynamic state with information about
 * distance from equilibrium and nonlinear response regime.
 */
typedef struct {
    ThermodynamicState state;
    double distance_from_equilibrium; /**< Norm of force vector |X| */
    double linearity_measure;         /**< |J - L X| / |J| -- deviation from linear */
    int    regime;                    /**< 0=equilibrium, 1=linear, 2=nonlinear, 3=far */
    double max_entropy_prod_estimate; /**< Upper bound on sigma [W/K] */
} FarFromEquilibriumState;

/* ---------------------------------------------------------------------------
 * L2: Core Concepts -- Nonlinear extensions of MEPP
 * ------------------------------------------------------------------------- */

/**
 * @brief Compute fluxes from nonlinear constitutive relations.
 *
 * J_i = sum_j L_{ij} X_j + sum_{j,k} M_{ijk} X_j X_k
 *     + sum_{j,k,l} N_{ijkl} X_j X_k X_l
 *
 * @param nl     Nonlinear Onsager coefficients
 * @param forces Force vector (length nl->dim)
 * @param fluxes [out] Flux vector (length nl->dim)
 * @return       0 on success
 */
int nt_compute_nonlinear_fluxes(const NonlinearOnsager *nl,
                                 const double forces[], double fluxes[]);

/**
 * @brief Compute generalized entropy production for nonlinear systems.
 *
 * sigma = sum_i J_i X_i with J_i computed from nonlinear relations.
 *
 * @param nl     Nonlinear Onsager coefficients
 * @param forces Force vector
 * @param ep     [out] Nonlinear entropy production breakdown
 * @return       0 on success
 */
int nt_nonlinear_sigma(const NonlinearOnsager *nl,
                       const double forces[], NonlinearEntropyProduction *ep);

/**
 * @brief Check if the Second Law holds for a nonlinear constitutive relation.
 *
 * Tests sigma(X) >= 0 for a grid of force values within a specified
 * range. Since sigma is not a quadratic form, positivity must be
 * verified numerically.
 *
 * @param nl       Nonlinear Onsager coefficients
 * @param X_max    Maximum force magnitude to test
 * @param n_grid   Number of grid points per dimension
 * @return         1 if sigma >= 0 everywhere tested, 0 if violation found
 */
int nt_check_nonlinear_second_law(const NonlinearOnsager *nl,
                                   double X_max, int n_grid);

/* ---------------------------------------------------------------------------
 * L3: Mathematical Structures -- Nonlinear stability and evolution
 * ------------------------------------------------------------------------- */

/**
 * @brief Compute the Jacobian of the nonlinear flux-force relation.
 *
 * For nonlinear relations, the Jacobian is state-dependent:
 *   J_{ij}(X) = partial(J_i) / partial(X_j)
 *            = L_{ij} + sum_k (M_{ijk} + M_{ikj}) X_k + ...
 *
 * This state-dependent Jacobian determines local stability properties
 * that can change as the system moves farther from equilibrium.
 *
 * @param nl       Nonlinear Onsager coefficients
 * @param forces   Force vector at which to evaluate Jacobian
 * @param jacobian [out] Jacobian matrix (n x n, row-major)
 * @return         0 on success
 */
int nt_nonlinear_jacobian(const NonlinearOnsager *nl,
                          const double forces[], double *jacobian);

/**
 * @brief Compute the evolution of forces under nonlinear dynamics.
 *
 * dX_i/dt = -gamma * partial(sigma) / partial(X_i) + boundary_terms
 *
 * where sigma includes nonlinear terms. The gradient structure is
 * preserved even with nonlinear Onsager relations, but the evolution
 * may no longer be a simple exponential relaxation.
 *
 * @param forces      Current force values
 * @param dforces_dt  [out] Time derivatives
 * @param nl          Nonlinear Onsager coefficients
 * @param gamma       Relaxation rate [1/s]
 * @param fixed       Fixed boundary forces (NaN for free)
 */
void nt_nonlinear_force_evolution(const double forces[], double dforces_dt[],
                                  const NonlinearOnsager *nl, double gamma,
                                  const double fixed[]);

/**
 * @brief Find stationary states of the nonlinear evolution equations.
 *
 * Solves dX_i/dt = 0 for the nonlinear gradient flow, which may
 * have multiple solutions (multiple steady states) far from
 * equilibrium, in contrast to the unique MEP state in the linear regime.
 *
 * @param forces   [in/out] Initial guess, updated to solution
 * @param nl       Nonlinear Onsager coefficients
 * @param fixed    Fixed forces (NaN for free)
 * @param max_iter Maximum Newton iterations
 * @param tol      Convergence tolerance
 * @return         Number of solutions found (1 or more), -1 on failure
 */
int nt_find_nonlinear_steady_states(double forces[], const NonlinearOnsager *nl,
                                    const double fixed[], int max_iter, double tol);

/* ---------------------------------------------------------------------------
 * L4: Fundamental Laws -- Maximum Entropy Production Principle debate
 * ------------------------------------------------------------------------- */

/**
 * @brief Evaluate the local potential for nonlinear systems.
 *
 * Glansdorff and Prigogine introduced the "local potential" as a
 * generalization of MEPP to nonlinear systems. The local potential
 * is not a true potential in the sense of being globally minimized,
 * but its variation satisfies an evolution criterion:
 *
 *   d_X sigma <= 0
 *
 * This function computes d_X sigma for a given state and its
 * time derivative.
 *
 * @param forces      Current forces
 * @param dforces_dt  Force time derivatives
 * @param nl          Nonlinear Onsager coefficients
 * @return            d_X sigma (should be <= 0)
 *
 * Reference: Glansdorff & Prigogine (1971), Ch. 7
 */
double nt_local_potential_evolution(const double forces[],
                                     const double dforces_dt[],
                                     const NonlinearOnsager *nl);

/**
 * @brief Compute the maximum entropy production state (MEP state).
 *
 * In contrast to MEPP (minimum), the Maximum Entropy Production
 * Principle proposes that some far-from-equilibrium systems select
 * the steady state that MAXIMIZES entropy production among the
 * available stable states.
 *
 * This function searches for local maxima of sigma(X) subject to
 * constraints, using gradient ascent.
 *
 * @param forces      [in/out] Force values
 * @param nl          Nonlinear Onsager coefficients
 * @param constrained Constraint flags
 * @param fixed_vals  Fixed force values
 * @param step_size   Gradient step size
 * @param max_iter    Maximum iterations
 * @param tol         Convergence tolerance
 * @return            0 if a local maximum found, -1 otherwise
 *
 * Reference: Martyushev & Seleznev (2006), Phys. Rep. 426, 1-45
 */
int nt_find_max_entropy_state(double forces[], const NonlinearOnsager *nl,
                               const int constrained[], const double fixed_vals[],
                               double step_size, int max_iter, double tol);

/* ---------------------------------------------------------------------------
 * L5: Algorithms -- Numerical methods for nonlinear thermodynamics
 * ------------------------------------------------------------------------- */

/**
 * @brief Fit a nonlinear Onsager model to force-flux data.
 *
 * Given m observations of (X, J), estimates coefficients L, M, N
 * by polynomial regression with symmetry constraints where applicable.
 *
 * @param forces  Force observations (m x n, row-major)
 * @param fluxes  Flux observations (m x n, row-major)
 * @param m       Number of observations
 * @param n       System dimension
 * @param order   Maximum polynomial order (2 or 3)
 * @param nl      [out] Fitted nonlinear Onsager coefficients
 * @return        0 on success
 */
int nt_fit_nonlinear_onsager(const double *forces, const double *fluxes,
                              int m, int n, int order, NonlinearOnsager *nl);

/**
 * @brief Compute the thermodynamic force magnitude for regime classification.
 *
 * |X| = sqrt(sum_i X_i^2)
 *
 * This metric determines whether the system is in the linear regime
 * (|X| small), nonlinear near-equilibrium (|X| moderate), or far
 * from equilibrium (|X| large).
 *
 * @param state      [in/out] State, distance_from_equilibrium updated
 * @param forces     Force vector
 * @param L          Onsager matrix (for linearity check)
 */
void nt_classify_regime(FarFromEquilibriumState *state,
                        const double forces[], const OnsagerMatrix *L);

/* Utility */

void nt_init_nonlinear_onsager(NonlinearOnsager *nl, int dim, int order);
void nt_free_nonlinear_onsager(NonlinearOnsager *nl);

#ifdef __cplusplus
}
#endif

#endif /* NONLINEAR_THERMO_H */