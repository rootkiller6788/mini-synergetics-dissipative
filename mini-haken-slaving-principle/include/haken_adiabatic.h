#ifndef HAKEN_ADIABATIC_H
#define HAKEN_ADIABATIC_H
#include "haken_core.h"

/* ==============================================================
 * haken_adiabatic.h — Adiabatic Elimination (Slaving Computation)
 *
 * The slaving principle in practice: eliminate fast variables
 * to obtain closed equations for the order parameters.
 *
 * Given the decomposition:
 *   dq_s/dt = Λ_s q_s + N_s(q_s, q_f)    (slow)
 *   dq_f/dt = Λ_f q_f + N_f(q_s, q_f)    (fast)
 *
 * with ||Λ_s|| ≪ ||Λ_f|| (spectral gap).
 *
 * Adiabatic elimination: assume dq_f/dt ≈ 0, solve:
 *   0 = Λ_f q_f + N_f(q_s, q_f)
 *   →  q_f = -Λ_f^{-1} N_f(q_s, 0) + O(||q_s||³)
 *
 * This yields the slaving relation q_f = h(q_s).
 *
 * Substituting back gives the order parameter equation:
 *   dq_s/dt = Λ_s q_s + N_s(q_s, h(q_s)) = f_eff(q_s)
 *
 * Error of adiabatic approximation is O(1/|spectral_gap|).
 *
 * References:
 *   Haken (1977) Synergetics, §7.2-7.4
 *   Gardiner (2004) Handbook of Stochastic Methods, §6.6
 *   Van Kampen (1985) "Elimination of Fast Variables", Phys. Rep.
 * ============================================================== */

/* -----------------------------------------------------------------
 * L1: Adiabatic elimination types
 * ----------------------------------------------------------------- */

/** Order of the adiabatic expansion */
typedef enum {
    HAKEN_ADIABATIC_LINEAR   = 0,  /* q_f ≈ -Λ_f^{-1} N_f(q_s,0), first order */
    HAKEN_ADIABATIC_QUADRATIC = 1, /* Include O(q_s²) corrections */
    HAKEN_ADIABATIC_CUBIC     = 2, /* Include O(q_s³) corrections */
    HAKEN_ADIABATIC_EXACT     = 3  /* Solve fixed-point equation iteratively */
} Haken_AdiabaticOrder;

/** Configuration for adiabatic elimination */
typedef struct {
    Haken_AdiabaticOrder order;       /* Expansion order */
    double  tolerance;                 /* Convergence tol for exact solve */
    int     max_iterations;            /* Max Newton iterations for exact */
    bool    check_validity;            /* Verify small-error condition */
    double  max_allowed_error;         /* Threshold for validity check */
} Haken_AdiabaticConfig;

/** Result of adiabatic elimination for one fast mode */
typedef struct {
    int     fast_index;         /* Index of fast mode */
    double* slaving_coeffs;     /* Coefficients of expansion in q_s */
    int     n_coeffs;           /* Number of coefficients */
    double  residual;           /* Residual of the elimination */
    double  estimate_error;     /* Estimated error from truncation */
    bool    converged;          /* Did the elimination converge? */
} Haken_AdiabaticMode;

/* -----------------------------------------------------------------
 * L4: Slaving computation API
 * ----------------------------------------------------------------- */

/** Perform adiabatic elimination on a Haken system.
 *  Given critical and stable mode separation, compute the slaving
 *  relation q_f = h(q_s) for all stable modes.
 *
 *  Algorithm:
 *    1. Project dynamics onto critical/stable subspaces
 *    2. Set dq_f/dt = 0
 *    3. Solve for q_f as function of q_s using perturbation expansion
 *    4. Build the reduced order parameter equation
 *
 *  Theorem (Haken 1977, §7.2):
 *    If there exists a spectral gap γ > 0 such that
 *    max Re(λ_crit) - max Re(λ_stable) = γ,
 *    then the adiabatic elimination error is O(1/γ).
 *
 *  @param sys       System with classified modes
 *  @param config     Adiabatic elimination configuration
 *  @param result     Output: slaving relation and reduced dynamics
 *  @return          0 on success, nonzero on failure */
int haken_adiabatic_eliminate(const Haken_System* sys,
                               const Haken_AdiabaticConfig* config,
                               Haken_SlavingResult* result);

/** Compute the slaving relation at a given order parameter value.
 *  Evaluates q_f = h(q_s) using the computed slaving coefficients.
 *
 *  @param slaving  Slaving relation (from elimination)
 *  @param q_s      Current order parameter values
 *  @param q_f      Output: fast variable values */
void haken_evaluate_slaving(const Haken_SlavingRelation* slaving,
                             const double* q_s, int n_slow,
                             double* q_f, int n_fast);

/** Compute the reduced (effective) dynamics for order parameters:
 *  dξ/dt = Λ_s ξ + N_s(ξ, h(ξ))
 *
 *  @param sys      System with decomposition
 *  @param slaving  Computed slaving relation
 *  @param xi       Current order parameter
 *  @param dxi      Output: time derivative */
void haken_reduced_dynamics(const Haken_System* sys,
                             const Haken_SlavingRelation* slaving,
                             const double* xi, double* dxi);

/** Validate the adiatic approximation by computing the residual:
 *  R = ||q_f - h(q_s)|| under the full dynamics.
 *
 *  @param sys        System
 *  @param slaving    Computed slaving relation
 *  @param dt         Time step for residual check
 *  @return           Max residual over a short integration */
double haken_adiabatic_residual(const Haken_System* sys,
                                 const Haken_SlavingRelation* slaving,
                                 double dt);

/* -----------------------------------------------------------------
 * L5: Perturbation expansion methods
 * ----------------------------------------------------------------- */

/** Compute the first-order slaving coefficients:
 *  q_f^{(1)} = -[Λ_f]^{-1} · N_f^{(2)}(q_s, q_s)
 *  where N_f^{(2)} is the quadratic part of the fast-mode nonlinearity.
 *
 *  This is the leading-order adiabatic elimination (Haken 1977, eq 7.2.14).
 *
 *  @param lambda_fast   Diagonal of fast eigenvalues
 *  @param n_fast        Number of fast modes
 *  @param n_slow        Number of slow modes
 *  @param coeff_matrix  Coupling coefficients N_f(a,b)
 *  @param slaving_coeff Output: first-order slaving coefficients */
void haken_first_order_slaving(const double* lambda_fast,
                                int n_fast, int n_slow,
                                const double* coeff_matrix,
                                double* slaving_coeff);

/** Compute second-order correction to slaving:
 *  q_f^{(2)} includes terms from N_f^{(3)}(q_s,q_s,q_s) and
 *  N_f^{(2)}(q_s, q_f^{(1)}).
 *
 *  Reference: Kuramoto (1984) Chemical Oscillations, §5.2 */
void haken_second_order_slaving(const double* lambda_fast,
                                 int n_fast, int n_slow,
                                 const double* coeff_quad,
                                 const double* coeff_cubic,
                                 const double* first_order,
                                 double* second_order);

/** Compute the effective cubic coefficient (saturation) in the
 *  order parameter equation after eliminating fast modes.
 *
 *  The cubic term determines whether the bifurcation is
 *  supercritical (saturating, stable) or subcritical (destabilizing).
 *
 *  Reference: Cross & Greenside (2009) Pattern Formation, §5 */
double haken_effective_cubic_coeff(const Haken_System* sys,
                                    const Haken_SlavingRelation* slaving,
                                    int critical_idx);

/* -----------------------------------------------------------------
 * Newton iteration for exact slaving (AdiabaticOrder = EXACT)
 * ----------------------------------------------------------------- */

/** Solve the fixed-point equation q_f = -Λ_f^{-1}·N_f(q_s, q_f)
 *  using Newton-Kantorovich iteration.
 *
 *  f(q_f) = q_f + Λ_f^{-1}·N_f(q_s, q_f) = 0
 *  J = I + Λ_f^{-1}·∂N_f/∂q_f  (Jacobian)
 *
 *  @param q_s       Fixed slow variables
 *  @param lambda_f  Fast eigenvalues (diagonal of Λ_f)
 *  @param n_f       Number of fast variables
 *  @param q_f_init  Initial guess for fast variables
 *  @param q_f_out   Output: converged fast variables
 *  @param config    Solver configuration
 *  @return          Number of iterations, -1 if failed */
int haken_newton_slaving(const double* q_s, int n_s,
                          const double* lambda_f, int n_f,
                          void (*N_func)(const double*, int, double*),
                          double* q_f_out,
                          const Haken_AdiabaticConfig* config);

/* -----------------------------------------------------------------
 * Error estimation
 * ----------------------------------------------------------------- */

/** Estimate the adiabatic elimination error:
 *  ε_adiabatic = ||q_f(exact) - q_f(adiabatic)|| / ||q_f(exact)||
 *
 *  Bounded by O(1/γ) where γ is the spectral gap.
 *  Reference: Van Kampen (1985) Phys. Rep. 124, 69 */
double haken_adiabatic_error_estimate(const Haken_System* sys,
                                       const Haken_AdiabaticConfig* config);

/** Solve linear system A·x = b via Gaussian elimination with
 *  partial pivoting. A is n×n, b and x are length n. */
void haken_solve_linear_system(const double* A, int n,
                                const double* b, double* x);

#endif /* HAKEN_ADIABATIC_H */
