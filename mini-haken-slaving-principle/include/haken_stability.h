#ifndef HAKEN_STABILITY_H
#define HAKEN_STABILITY_H
#include "haken_core.h"

/* ==============================================================
 * haken_stability.h — Linear Stability & Bifurcation Analysis
 *
 * At the core of the slaving principle lies the separation of
 * timescales near an instability. This module provides the tools
 * to analyze stability, detect bifurcations, and identify the
 * critical modes that become order parameters.
 *
 * Linear stability:
 *   d/dt δq = J(α,q*)·δq, where J = ∂f/∂q is the Jacobian
 *   λ_k(J) determine local stability
 *
 * Critical slowing down:
 *   As α → α_c, the relaxation time τ = -1/Re(λ_crit) → ∞
 *   This is the dynamical signature of an impending phase transition.
 *
 * References:
 *   Haken (1977) Synergetics, §5.1-5.3
 *   Strogatz (2015) Nonlinear Dynamics and Chaos, §3
 *   Kuznetsov (2004) Elements of Applied Bifurcation Theory
 * ============================================================== */

/* -----------------------------------------------------------------
 * L2: Stability classification
 * ----------------------------------------------------------------- */

/** Lyapunov stability type */
typedef enum {
    HAKEN_STABLE_NODE      = 0,  /* All λ < 0, real */
    HAKEN_STABLE_FOCUS     = 1,  /* All Re(λ) < 0, some complex */
    HAKEN_UNSTABLE_NODE    = 2,  /* All λ > 0, real */
    HAKEN_UNSTABLE_FOCUS   = 3,  /* Some Re(λ) > 0, complex */
    HAKEN_SADDLE           = 4,  /* Both positive and negative λ */
    HAKEN_CENTER           = 5,  /* All Re(λ) = 0 */
    HAKEN_DEGENERATE       = 6   /* Defective/non-diagonalizable */
} Haken_StabilityType;

/** Fixed point with stability information */
typedef struct {
    double*  position;      /* Fixed point coordinates, length n_dim */
    int      n_dim;
    double*  jacobian;      /* Jacobian at fixed point, n_dim×n_dim */
    double*  eigenvalues_re; /* Eigenvalues of Jacobian */
    double*  eigenvalues_im;
    Haken_StabilityType type;
    bool     is_hyperbolic;  /* No eigenvalues on imaginary axis */
    int      unstable_dim;   /* Number of eigenvalues with Re(λ)>0 */
    int      center_dim;     /* Number of eigenvalues with Re(λ)=0 */
} Haken_FixedPoint;

/** Bifurcation detection result */
typedef struct {
    Haken_BifurcationType type;
    double   critical_value;  /* Control parameter at bifurcation */
    int     *critical_modes;  /* Indices of modes involved */
    int      n_critical_modes;
    double   crossing_speed;  /* d(Re λ)/dα at bifurcation */
    bool     is_degenerate;   /* Codimension > 1 */
    double   normal_form_coeff; /* a in ẋ = αx + a·x³ */
    char     description[256]; /* Human-readable description */
} Haken_Bifurcation;

/* -----------------------------------------------------------------
 * L4: Stability analysis API
 * ----------------------------------------------------------------- */

/** Compute the Jacobian matrix J_{ij} = ∂f_i/∂q_j using central
 *  finite differences at a given state.
 *  J_{ij} ≈ [f_i(q + h·e_j) - f_i(q - h·e_j)] / (2h)
 *  O(n²·eval) time. */
void haken_compute_jacobian(const Haken_System* sys, const double* q,
                             double* J, double h);

/** Find fixed points of the dynamics: f(q*) = 0
 *  using Newton-Raphson with line search.
 *  @param sys        System descriptor
 *  @param q_initial  Initial guess for fixed point
 *  @param q_fixed    Output: fixed point coordinates
 *  @param max_iter   Maximum Newton iterations
 *  @param tol        Convergence tolerance
 *  @return           0 on success, number of iterations on failure (> max_iter) */
int haken_find_fixed_point(const Haken_System* sys,
                            const double* q_initial,
                            double* q_fixed,
                            int max_iter, double tol);

/** Analyze stability of a fixed point.
 *  Computes Jacobian, eigenvalues, and classifies stability type.
 *  @param sys    System
 *  @param q_star Fixed point position
 *  @param fp     Output: filled FixedPoint structure */
void haken_analyze_fixed_point(const Haken_System* sys,
                                const double* q_star,
                                Haken_FixedPoint* fp);

/** Free a FixedPoint structure */
void haken_fixed_point_free(Haken_FixedPoint* fp);

/** Detect the type of bifurcation at a given control parameter.
 *  Analyzes eigenvalue crossing pattern.
 *  @param sys      System (with computed modes)
 *  @param prev_sys Previous system state (slightly different α)
 *  @param result   Output: bifurcation type and details */
void haken_detect_bifurcation(const Haken_System* sys,
                               const Haken_System* prev_sys,
                               Haken_Bifurcation* result);

/** Compute the normal form coefficient for a pitchfork bifurcation:
 *  ẋ = α·x + a·x³ + O(x⁵)
 *  a < 0 → supercritical (continuous), a > 0 → subcritical (discontinuous)
 *
 *  Uses center manifold reduction (Kuznetsov 2004, §5.4). */
double haken_pitchfork_normal_form(const Haken_System* sys,
                                    int critical_mode_idx);

/** Compute the first Lyapunov coefficient for a Hopf bifurcation:
 *  Determines supercritical (l₁<0) vs subcritical (l₁>0).
 *  Reference: Kuznetsov (2004) §5.4.3 */
double haken_hopf_lyapunov_coefficient(const Haken_System* sys,
                                        int mode_idx_re, int mode_idx_im);

/* -----------------------------------------------------------------
 * Critical slowing down analysis
 * ----------------------------------------------------------------- */

/** Compute the relaxation time for each mode:
 *  τ_k = -1/Re(λ_k)   (if Re(λ_k) < 0)
 *  τ_k = ∞            (if Re(λ_k) ≥ 0)
 *
 *  Critical slowing down: τ_crit → ∞ as α → α_c. */
void haken_relaxation_times(const Haken_System* sys, double* tau);

/** Estimate the critical parameter α_c by extrapolation:
 *  Fit τ(α) = τ₀·|α - α_c|^{-1} to relaxation time data.
 *  Uses linear regression in log-log space. */
double haken_estimate_critical_parameter(const double* alpha_values,
                                          const double* tau_values,
                                          int n_points,
                                          double* alpha_c_out);

/** Compute the dynamic critical exponent z:
 *  τ ~ |α - α_c|^{-νz},  ξ ~ |α - α_c|^{-ν}
 *  Estimates z from relaxation time data. */
double haken_dynamic_exponent(const double* alpha_values,
                               const double* tau_values,
                               int n_points, double alpha_c);

/* -----------------------------------------------------------------
 * Parameter sweep and tracking
 * ----------------------------------------------------------------- */

/** Allocate a parameter sweep structure */
Haken_ParameterSweep* haken_sweep_create(double start, double end,
                                          int n_points);

/** Free a parameter sweep */
void haken_sweep_free(Haken_ParameterSweep* sweep);

/** Perform a full parameter sweep: for each α, compute eigenvalues,
 *  mode classification, spectral gap, and order parameter amplitude.
 *  Stores results in the sweep structure. */
int haken_sweep_execute(Haken_System* sys, Haken_ParameterSweep* sweep);

/** Find the α that maximizes the spectral gap (optimal slaving regime). */
double haken_sweep_max_spectral_gap(const Haken_ParameterSweep* sweep);

/** Find the α at which the first eigenvalue crosses zero.
 *  Uses linear interpolation between sign changes. */
double haken_sweep_critical_threshold(const Haken_ParameterSweep* sweep,
                                       int mode_index);

/** Compute the susceptibility: χ = ∂ξ/∂h where h is an external field.
 *  For GL theory: χ ~ |α - α_c|^{-1} (mean-field).
 *  Estimates via finite-difference response to perturbation. */
double haken_susceptibility(const Haken_System* sys,
                             const Haken_OrderParam* op,
                             double perturbation_strength);

/** Compute the mode softening rate:
 *  S_k = d(Re λ_k)/dα at α = α_c
 *  Faster softening → larger rate → mode becomes critical more abruptly. */
double haken_mode_softening_rate(Haken_System* sys, int mode_idx,
                                  double dalpha);

/** Control parameter sweep with state tracking.
 *  Continuously varying α while following the equilibrium state.
 *  Uses predictor-corrector continuation (Allgower & Georg 1990). */
int haken_continuation(Haken_System* sys,
                        double alpha_start, double alpha_end,
                        int n_steps, double* states_out);

#endif /* HAKEN_STABILITY_H */
