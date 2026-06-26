/**
 * @file bz_stability.h
 * @brief Linear and nonlinear stability analysis for BZ reaction models.
 *        Jacobian eigenvalues, Hopf bifurcation, limit cycle stability.
 *
 * Knowledge coverage:
 *   L2: Linear stability, asymptotic stability, instability
 *   L3: Eigenvalue problem, characteristic polynomial
 *   L4: Hopf bifurcation theorem, Lyapunov indirect method
 *   L5: Eigenvalue computation, Floquet analysis
 *   L8: Bifurcation tracking, stability exchange
 */

#ifndef BZ_STABILITY_H
#define BZ_STABILITY_H

#include "bz_reaction.h"

/* ---------------------------------------------------------------------------
 * L2-L3: Linear stability analysis
 * --------------------------------------------------------------------------- */

/** Stability classification of a steady state based on eigenvalues. */
typedef enum {
    BZ_STABILITY_STABLE_NODE       = 0,
    BZ_STABILITY_STABLE_FOCUS      = 1,
    BZ_STABILITY_UNSTABLE_NODE     = 2,
    BZ_STABILITY_UNSTABLE_FOCUS    = 3,
    BZ_STABILITY_SADDLE            = 4,
    BZ_STABILITY_CENTER            = 5,
    BZ_STABILITY_DEGENERATE        = 6
} bz_stability_type_t;

/** Classify stability of a 3x3 Jacobian from its eigenvalues.
 *  @param ev_re  Real parts of eigenvalues
 *  @param ev_im  Imaginary parts of eigenvalues
 *  @return       Stability classification */
bz_stability_type_t bz_classify_stability(const double *ev_re,
                                           const double *ev_im);

/** Compute spectral radius of a 3x3 matrix: rho(J) = max_i |lambda_i|
 *  @param J  3x3 matrix (flat, row-major)
 *  @return   Spectral radius */
double bz_spectral_radius_3x3(const double *J);

/** Check if steady state is linearly stable (all eigenvalues Re < 0).
 *  @param J  3x3 Jacobian matrix
 *  @return   1 if stable, 0 otherwise */
int bz_is_linearly_stable(const double *J);

/* ---------------------------------------------------------------------------
 * L4: Hopf bifurcation analysis
 * --------------------------------------------------------------------------- */

/** Compute first Lyapunov coefficient for Oregonator Hopf bifurcation.
 *  Determines supercritical (l1 < 0) vs subcritical (l1 > 0).
 *  Uses projection method (Kuznetsov formula).
 *  @param p       Parameters at bifurcation
 *  @param ss      Steady state at bifurcation
 *  @param omega   Frequency at bifurcation
 *  @return        First Lyapunov coefficient l1 */
double bz_first_lyapunov_coeff(const bz_oregonator_params_t *p,
                                const bz_oregonator_state_t *ss,
                                double omega);

/** Determine Hopf bifurcation type.
 *  @param p       Parameters near bifurcation
 *  @param ss      Steady state at bifurcation
 *  @return        1 = supercritical, -1 = subcritical, 0 = degenerate */
int bz_hopf_type(const bz_oregonator_params_t *p,
                  const bz_oregonator_state_t *ss);

/* ---------------------------------------------------------------------------
 * L5: Limit cycle stability (Floquet analysis)
 * --------------------------------------------------------------------------- */

/** Compute Floquet multipliers for a limit cycle via Poincare map.
 *  @param multipliers  Output: Floquet multipliers (excluding trivial 1)
 *  @param n_mult       Output: number of non-trivial multipliers
 *  @param crossings_z  Poincare section crossing z-values
 *  @param n_cross      Number of crossings
 *  @return             1 if limit cycle is stable (all |mult| < 1) */
int bz_floquet_stability(double *multipliers, int *n_mult,
                          const double *crossings_z, int n_cross);

/** Estimate largest Floquet exponent from trajectory divergence.
 *  lambda_F = (1/T) * log(||delta(T)|| / ||delta(0)||)
 *  @param traj1     Reference trajectory
 *  @param traj2     Perturbed trajectory
 *  @param n         Number of points
 *  @param dt        Time step
 *  @param exponent  Output: largest Floquet exponent
 *  @return          1 if computed successfully */
int bz_floquet_exponent(double *exponent,
                         const double *traj1, const double *traj2,
                         int n, double dt);

/* ---------------------------------------------------------------------------
 * L8: Advanced: Bifurcation tracking
 * --------------------------------------------------------------------------- */

/** Track Hopf bifurcation curve in (f, q) parameter space.
 *  Uses secant continuation.
 *  @param f_vals    Output f values along curve
 *  @param q_vals    Output q values along curve
 *  @param n_pts     Number of points on curve
 *  @param q_start   Starting q value */
void bz_track_hopf_curve(double *f_vals, double *q_vals, int n_pts,
                          double q_start);

/** Detect Bogdanov-Takens bifurcation point (double-zero eigenvalue).
 *  @param J         3x3 Jacobian at candidate point
 *  @return          1 if Bogdanov-Takens, 0 otherwise */
int bz_detect_bogdanov_takens(const double *J);

/** Compute the criticality of the Hopf bifurcation.
 *  @param p    Parameters
 *  @param ss   Steady state
 *  @return     Supercritical (1) or subcritical (-1) */
int bz_hopf_criticality(const bz_oregonator_params_t *p,
                         const bz_oregonator_state_t *ss);

#endif /* BZ_STABILITY_H */
