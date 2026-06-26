/**
 * @file bz_oregonator.h
 * @brief Oregonator model: Field-Noyes reduced model of the BZ reaction.
 *        Covers rigorous analysis of the 3-variable Oregonator ODE system.
 *
 * Knowledge coverage:
 *   L2: Time-scale separation, excitability, limit cycle
 *   L3: Dimensionless ODE system, Jacobian matrix
 *   L4: Hopf bifurcation theorem applied to Oregonator
 *   L5: Parameter sweeps, bifurcation tracking
 *
 * The Oregonator (Field & Noyes, 1974) is a 3-variable ODE model obtained
 * by reducing the full FKN mechanism using time-scale separation.
 *
 * Dimensionless form:
 *   eps   * dx/dt = q*y - x*y + x*(1 - x)
 *   eps'  * dy/dt = -q*y - x*y + f*z
 *           dz/dt = x - z
 *
 * References:
 *   - Field & Noyes (1974) J. Chem. Phys. 60, 1877-1884
 *   - Tyson (1985) "Oscillations and Traveling Waves in Chemical Systems"
 *   - Murray (2002) "Mathematical Biology" Vol I, Ch 7
 */

#ifndef BZ_OREGONATOR_H
#define BZ_OREGONATOR_H

#include "bz_reaction.h"

/* ---------------------------------------------------------------------------
 * L2: Oregonator dynamical analysis -- limit cycle detection
 * --------------------------------------------------------------------------- */

/** Oregonator oscillation period approximation.
 *  Uses scaling estimate: T ~ (1/epsilon_p)*log(f/(f-1)) for f > 1.
 *  More accurate for small epsilon. Returns period in dimensionless time.
 *  @param p   Oregonator parameters
 *  @return    Estimated oscillation period */
double bz_oregonator_period_approx(const bz_oregonator_params_t *p);

/** Detect if parameters correspond to an oscillatory regime.
 *  Uses the condition from Hopf bifurcation analysis:
 *    oscillatory if f_min(q) < f < f_max(q), where
 *    f_min = (1/2)*(1+sqrt(2)) when q << 1.
 *  @param p   Oregonator parameters
 *  @return    1 if oscillatory, 0 otherwise */
int bz_oregonator_is_oscillatory(const bz_oregonator_params_t *p);

/** Compute the critical stoichiometric factor f_c at which the Hopf
 *  bifurcation occurs (transition from steady state to oscillation).
 *  For Oregonator with q small:
 *    f_c ¡Ö 0.5 + 0.5*sqrt(1 + 4*q)
 *  @param q   Rate constant ratio
 *  @return    Critical stoichiometric factor f_c */
double bz_oregonator_hopf_critical_f(double q);

/** Characterize oscillation shape: relaxation vs sinusoidal.
 *  Returns ratio of fast to slow time scale.
 *  Ratio >> 1 => relaxation oscillations (sharp spikes).
 *  Ratio ~ 1  => near-sinusoidal (supercritical Hopf).
 *  @param p   Oregonator parameters
 *  @return    Relaxation ratio (>1 means relaxation) */
double bz_oregonator_relaxation_ratio(const bz_oregonator_params_t *p);

/* ---------------------------------------------------------------------------
 * L3: Jacobian matrix and linear stability
 * --------------------------------------------------------------------------- */

/** Oregonator Jacobian matrix at a given state.
 *  J[0..8] row-major:
 *    row0: d(dx/dt)/dx, d(dx/dt)/dy, d(dx/dt)/dz
 *    row1: d(dy/dt)/dx, d(dy/dt)/dy, d(dy/dt)/dz
 *    row2: d(dz/dt)/dx, d(dz/dt)/dy, d(dz/dt)/dz
 *  @param J    Output: 3x3 Jacobian (flat, 9 doubles)
 *  @param s    State at which Jacobian is evaluated
 *  @param p    Parameters */
void bz_oregonator_jacobian(double *J,
                             const bz_oregonator_state_t *s,
                             const bz_oregonator_params_t *p);

/** Compute eigenvalues of a 3x3 real matrix.
 *  Solves the cubic characteristic equation analytically.
 *  @param ev_re  Output: real parts (3 doubles)
 *  @param ev_im  Output: imaginary parts (3 doubles)
 *  @param J      3x3 matrix (flat, row-major)
 *  @return       Number of eigenvalues with Re(lambda) > 0 */
int bz_eigenvalues_3x3(double *ev_re, double *ev_im, const double *J);

/** Compute characteristic polynomial coefficients of a 3x3 matrix.
 *  det(J - lambda*I) = -lambda^3 + tr*lambda^2 - m2*lambda + det
 *  @param trace  Output: trace(J)
 *  @param m2     Output: sum of principal minors
 *  @param det    Output: determinant
 *  @param J      3x3 matrix (flat, row-major) */
void bz_char_poly_3x3(double *trace, double *m2, double *det, const double *J);

/* ---------------------------------------------------------------------------
 * L4: Hopf bifurcation analysis
 * --------------------------------------------------------------------------- */

/** Assess whether a Hopf bifurcation occurs at the given steady state.
 *  Checks: (1) pair of complex conjugate eigenvalues crosses imaginary axis,
 *  (2) transversality condition d(Re lambda)/df != 0.
 *  @param p      Oregonator parameters
 *  @param ss     Steady state (from bz_oregonator_steady_states)
 *  @return       1 if Hopf bifurcation detected, 0 otherwise */
int bz_oregonator_hopf_test(const bz_oregonator_params_t *p,
                             const bz_oregonator_state_t *ss);

/** Compute the Hopf bifurcation curve f = f_H(q) in parameter space.
 *  For the Oregonator with epsilon << 1:
 *    f_H = (1 + q)^2 / (3*q - 1)   for q > 1/3
 *  @param q   Rate constant ratio
 *  @return    f value at Hopf bifurcation */
double bz_oregonator_hopf_curve(double q);

/* ---------------------------------------------------------------------------
 * L5: Parameter sweep and bifurcation tracking
 * --------------------------------------------------------------------------- */

/** Perform a parameter sweep of f to detect oscillation onset.
 *  For each f, integrates Oregonator and checks for threshold crossings.
 *  @param f_min     Minimum f
 *  @param f_max     Maximum f
 *  @param n_steps   Number of sweep points
 *  @param results   Output: results[i] = 1 if oscillatory at step i
 *  @param p_base    Base parameters (f will be overridden) */
void bz_oregonator_sweep_f(double f_min, double f_max, int n_steps,
                            int *results, const bz_oregonator_params_t *p_base);

/* ---------------------------------------------------------------------------
 * L6: Canonical problem -- Oregonator oscillation integration
 * --------------------------------------------------------------------------- */

/** Integrate Oregonator and analyze oscillation.
 *  Uses RK4 with fixed step size.
 *  @param s            Initial state (modified to final state)
 *  @param p            Parameters
 *  @param duration     Total integration time
 *  @param n_steps      Number of integration steps
 *  @param n_peaks      Output: number of peaks detected
 *  @param period_avg   Output: average period
 *  @param amplitude    Output: peak-to-trough amplitude of x */
void bz_oregonator_integrate_and_analyze(bz_oregonator_state_t *s,
                                          const bz_oregonator_params_t *p,
                                          double duration, int n_steps,
                                          int *n_peaks,
                                          double *period_avg,
                                          double *amplitude);

#endif /* BZ_OREGONATOR_H */
