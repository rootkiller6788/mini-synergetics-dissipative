/**
 * turing_analysis.h — Linear Stability & Bifurcation Analysis
 *
 * Reference: Turing, A.M. (1952) "The Chemical Basis of Morphogenesis"
 *            Murray, J.D. (2003) Mathematical Biology II, Ch. 2
 *            Cross, M.C. & Hohenberg, P.C. (1993) Rev. Mod. Phys. 65(3):851
 *            Maini, P.K. et al. (1997) Bull. Math. Biol. 59(1):61-84
 *
 * Provides the complete toolkit for analyzing Turing instability:
 *   1. Homogeneous steady state (HSS) computation via Newton's method
 *   2. Linear stability analysis of HSS against Fourier modes
 *   3. Dispersion relation λ(k) computation
 *   4. Turing space construction (parameter regions with instability)
 *   5. Bifurcation detection (Turing, Hopf, Turing-Hopf interaction)
 *
 * Knowledge Coverage:
 *   L2: Diffusion-driven instability, Turing conditions
 *   L3: Jacobian eigenvalues, dispersion relation, Fourier mode growth
 *   L4: Turing instability theorem (4 necessary conditions)
 *   L5: Newton's method for HSS, eigenvalue computation, parameter scanning
 */

#ifndef TURING_ANALYSIS_H
#define TURING_ANALYSIS_H

#include "turing_common.h"
#include "turing_models.h"

/* ---------------------------------------------------------------------------
 * L4: Turing Instability Conditions
 *
 * For a reaction-diffusion system ∂u/∂t = f(u,v) + D_u∇²u, ∂v/∂t = g(u,v) + D_v∇²v,
 * the homogeneous steady state (u*,v*) is Turing-unstable iff:
 *
 * Condition 1: f_u + g_v < 0                    (stable to homogeneous perturbations)
 * Condition 2: f_u·g_v - f_v·g_u > 0            (det(J) > 0)
 * Condition 3: D_v·f_u + D_u·g_v > 0            (unstable to some k)
 * Condition 4: (D_v·f_u + D_u·g_v)² > 4D_u·D_v·det(J)  (real eigenvalues at max)
 *
 * where J = [[f_u, f_v], [g_u, g_v]] at (u*,v*).
 *
 * These 4 conditions are necessary and sufficient for Turing instability.
 * Reference: Murray (2003), §2.3
 * -------------------------------------------------------------------------*/

/** Structure holding the four Turing conditions for diagnostics. */
typedef struct {
    double  f_u;            /* ∂f/∂u at HSS */
    double  f_v;            /* ∂f/∂v at HSS */
    double  g_u;            /* ∂g/∂u at HSS */
    double  g_v;            /* ∂g/∂v at HSS */
    double  trace_J;        /* f_u + g_v */
    double  det_J;          /* f_u*g_v - f_v*g_u */
    int     cond1;          /* trace_J < 0 ? 1 : 0 */
    int     cond2;          /* det_J > 0 ? 1 : 0 */
    int     cond3;          /* Dv*f_u + Du*g_v > 0 ? 1 : 0 */
    int     cond4;          /* discriminant > 0 ? 1 : 0 */
    double  discriminant;   /* (Dv*f_u + Du*g_v)² - 4*Du*Dv*det_J */
    int     turing_unstable; /* all 4 conditions satisfied */
} TuringConditions;

/* ---------------------------------------------------------------------------
 * L2: Turing Space
 *
 * A 2D parameter sweep showing regions of Turing instability.
 * Typical axes: (F, k) for Gray-Scott, (A, B) for Brusselator, etc.
 *
 * Each grid point (p1[i], p2[j]) is classified as:
 *   STABLE:       homogeneous steady state is stable to all perturbations
 *   TURING:       Turing-unstable (spatial pattern forms)
 *   HOPF:         Hopf-unstable (temporal oscillations)
 *   TURING_HOPF:  both instabilities (spatiotemporal complexity)
 *   NO_HSS:       no valid homogeneous steady state found
 * -------------------------------------------------------------------------*/
typedef enum {
    REGION_STABLE       = 0,
    REGION_TURING       = 1,
    REGION_HOPF         = 2,
    REGION_TURING_HOPF  = 3,
    REGION_NO_HSS       = 4
} TuringRegion;

/** 2D parameter scan result. */
typedef struct {
    double         p1_min, p1_max;   /* parameter 1 range */
    double         p2_min, p2_max;   /* parameter 2 range */
    int            n1, n2;           /* grid size */
    TuringRegion  *regions;          /* n1*n2 array, row-major */
    double         Du, Dv;           /* diffusion coefficients */
    ModelType      model;
    ModelParams    base_params;      /* all other parameters fixed */
} TuringSpace;

/* ---------------------------------------------------------------------------
 * L2: Bifurcation Type
 * -------------------------------------------------------------------------*/
typedef enum {
    BIF_NONE         = 0,
    BIF_TURING       = 1,
    BIF_HOPF         = 2,
    BIF_TURING_HOPF  = 3,   /* co-dimension 2 point */
    BIF_SADDLE_NODE  = 4
} BifurcationType;

/** Bifurcation point descriptor. */
typedef struct {
    BifurcationType  type;
    double           param1_value;
    double           param2_value;
    double           k_critical;
    double           omega_critical;  /* Hopf frequency if applicable */
} BifurcationPoint;

/* ---------------------------------------------------------------------------
 * Analysis Functions
 * -------------------------------------------------------------------------*/

/**
 * Compute the 4 Turing conditions for a given HSS, Jacobian, and diffusion coefficients.
 *
 * Input:
 *   J[2][2] — Jacobian matrix at HSS
 *   Du, Dv  — diffusion coefficients
 *
 * Output:
 *   conditions — filled with condition values
 *
 * Returns 0 on success.
 *
 * Reference: Murray (2003) §2.3, Eqs. (2.27)-(2.31)
 */
int turing_conditions_compute(const double J[2][2], double Du, double Dv,
                               TuringConditions *conditions);

/**
 * Print Turing conditions in human-readable format.
 */
void turing_conditions_print(const TuringConditions *tc);

/**
 * Compute the dispersion relation λ(k) for a given J, Du, Dv.
 *
 * For each wavenumber k:
 *   det(J - k²D - λI) = 0
 *
 * The eigenvalue with larger real part is λ(k).
 * The system is Turing-unstable if Re(λ(k)) > 0 for some k > 0.
 *
 * Input:
 *   J[2][2] — Jacobian at HSS
 *   Du, Dv  — diffusion coefficients
 *   n_modes — number of Fourier modes
 *   k_min, k_max — wavenumber range
 *
 * Returns newly allocated DispersionRelation.
 *
 * Complexity: O(n_modes) per evaluation
 * Reference: Cross & Greenside (2009) Pattern Formation and Dynamics, Ch. 5
 */
DispersionRelation *turing_dispersion_compute(const double J[2][2],
                                              double Du, double Dv,
                                              int n_modes,
                                              double k_min, double k_max);

/**
 * Print dispersion relation summary.
 */
void turing_dispersion_print(const DispersionRelation *dr);

/**
 * Find critical wavenumber k_c where Re(λ(k)) is maximized (fastest growing mode).
 * For a 2-component system, analytical formula exists:
 *   k_c² = sqrt(det_J / (Du * Dv))
 *
 * Reference: Murray (2003) §2.3, Eq. (2.28)
 */
double turing_critical_wavenumber(const double J[2][2], double Du, double Dv);

/**
 * Compute the predicted pattern wavelength from dispersion analysis:
 *   λ_predicted = 2π / k_critical
 */
double turing_predicted_wavelength(const DispersionRelation *dr);

/* ---------------------------------------------------------------------------
 * Homogeneous Steady State Solver
 * -------------------------------------------------------------------------*/

/**
 * Newton's method for finding homogeneous steady state of a reaction-diffusion system.
 *
 * Solves: f(u,v) = 0, g(u,v) = 0
 *
 * Iteration:
 *   [u_{n+1}]   [u_n]           [f(u_n, v_n)]
 *   [v_{n+1}] = [v_n] - J⁻¹(u_n,v_n) * [g(u_n, v_n)]
 *
 * Parameters:
 *   reaction_f — reaction rate function
 *   jacobian_f — Jacobian function
 *   p          — model parameters
 *   u0, v0     — initial guess
 *   max_iter   — maximum Newton iterations
 *   tol        — convergence tolerance
 *
 * Returns:
 *   Pointer to HomogeneousSteadyState (caller frees with free())
 *   NULL if method fails to converge.
 *
 * Complexity: O(max_iter) Newton steps, each O(1)
 * Reference: Kelley, C.T. (2003) Solving Nonlinear Equations with Newton's Method
 */
HomogeneousSteadyState *newton_find_hss(ReactionFunc reaction_f,
                                         JacobianFunc jacobian_f,
                                         const ModelParams *p,
                                         double u0, double v0,
                                         int max_iter, double tol);

/**
 * Find all real steady states by systematic multiple-start Newton search.
 * Returns array of HSS pointers, with *n_hss set to the count.
 * Caller must free each HSS and the array.
 */
HomogeneousSteadyState **newton_find_all_hss(ReactionFunc reaction_f,
                                              JacobianFunc jacobian_f,
                                              const ModelParams *p,
                                              int *n_hss);

/**
 * Compute eigenvalues of a 2x2 real matrix.
 * Stores real parts in eigen_real[2], imaginary parts in eigen_imag[2].
 * Returns 0 on success.
 */
int matrix_2x2_eigenvalues(const double M[2][2],
                           double eigen_real[2], double eigen_imag[2]);

/**
 * Check Hopf instability: trace(J) > 0 (with det(J) > 0).
 * Hopf bifurcation occurs when a pair of complex conjugate eigenvalues
 * crosses the imaginary axis.
 *
 * Returns 1 if Hopf-unstable, 0 otherwise.
 */
int is_hopf_unstable(const double J[2][2]);

/* ---------------------------------------------------------------------------
 * Turing Space Construction
 * -------------------------------------------------------------------------*/

/**
 * Allocate and initialize a Turing space for parameter scanning.
 *
 * Parameters:
 *   model      — kinetic model
 *   base       — base parameters (will be varied in p1, p2 axis)
 *   Du, Dv     — diffusion coefficients
 *   n1, n2     — grid resolution for each parameter
 *
 * p1_range/p2_range specify which parameters to vary and their range.
 * param1_idx: 0=F, 1=k, 2=A, 3=B, etc. (model-dependent interpretation)
 */
TuringSpace *turing_space_alloc(ModelType model, const ModelParams *base,
                                double Du, double Dv,
                                int param1_idx, double p1_min, double p1_max, int n1,
                                int param2_idx, double p2_min, double p2_max, int n2);

/** Free Turing space structure. */
void turing_space_free(TuringSpace *ts);

/**
 * Compute the Turing space by scanning parameters and classifying stability
 * at each grid point.
 *
 * For each (p1[i], p2[j]):
 *   1. Find HSS
 *   2. Compute Jacobian at HSS
 *   3. Check Turing and Hopf conditions
 *   4. Classify into region type
 *
 * Complexity: O(n1*n2) — each grid point requires Newton + eigenvalue analysis
 */
int turing_space_compute(TuringSpace *ts, const SolverConfig *ignored_var);

/**
 * Count grid points in each region type.
 * Output: counts[REGION_STABLE], counts[REGION_TURING], etc.
 */
void turing_space_count_regions(const TuringSpace *ts, int counts[5]);

/**
 * Print Turing space summary with region counts.
 */
void turing_space_print_summary(const TuringSpace *ts);

/* ---------------------------------------------------------------------------
 * Bifurcation Analysis
 * -------------------------------------------------------------------------*/

/**
 * Find Turing bifurcation points along a 1D parameter sweep.
 *
 * A Turing bifurcation occurs where max_k Re(λ(k)) crosses 0 from below.
 * Uses bisection to locate the critical parameter value precisely.
 *
 * Parameters:
 *   model          — the kinetic model
 *   base_params    — base parameters (one will be varied)
 *   vary_idx       — index of parameter to vary
 *   p_min, p_max   — search range
 *   Du, Dv         — diffusion coefficients
 *   max_bifurcations — max number to find
 *   n_found         — output: number of bifurcations found
 *
 * Returns array of BifurcationPoint (caller frees).
 */
BifurcationPoint *find_turing_bifurcations(ModelType model,
                                            const ModelParams *base_params,
                                            int vary_idx,
                                            double p_min, double p_max,
                                            double Du, double Dv,
                                            int max_bifurcations,
                                            int *n_found);

/**
 * Print bifurcation point details.
 */
void bifurcation_point_print(const BifurcationPoint *bp);

/** Check for Turing-Hopf co-dimension 2 point. */
int detect_turing_hopf_interaction(const double J[2][2], double Du, double Dv);

/* ---------------------------------------------------------------------------
 * Sensitivity Analysis
 *
 * How does pattern wavelength / growth rate vary with parameters?
 * -------------------------------------------------------------------------*/

/**
 * Compute ∂λ_max/∂param for each parameter via finite differences.
 *
 * Parameters:
 *   model        — kinetic model
 *   p_nominal    — nominal parameter set
 *   Du, Dv       — diffusion coefficients
 *   sensitivities — output array, size n_params (model-dependent)
 *   n_params     — number of parameters
 *
 * Returns 0 on success.
 */
int turing_sensitivity_analysis(ModelType model, const ModelParams *p_nominal,
                                double Du, double Dv,
                                double *sensitivities, int n_params);

/**
 * Compute pattern wavelength sensitivity to domain size changes.
 * Useful for understanding scale-dependence of patterns.
 */
double turing_wavelength_sensitivity(double Lx, double Ly,
                                     const double J[2][2],
                                     double Du, double Dv);

#endif /* TURING_ANALYSIS_H */
