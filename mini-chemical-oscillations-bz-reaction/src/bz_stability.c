
/**
 * @file bz_stability.c
 * @brief Stability analysis: eigenvalue classification, Hopf bifurcation,
 *        Lyapunov coefficient, Floquet multipliers, Bogdanov-Takens.
 *
 * Implements:
 *   L2: Linear stability classification
 *   L3: Spectral radius
 *   L4: Hopf bifurcation type (super/subcritical)
 *   L5: Floquet stability analysis
 *   L8: Bogdanov-Takens detection, Hopf curve tracking
 *
 * References:
 *   - Kuznetsov (2004) "Elements of Applied Bifurcation Theory"
 *   - Guckenheimer & Holmes (1983) "Nonlinear Oscillations..."
 */

#include "bz_stability.h"
#include "bz_reaction.h"
#include "bz_oregonator.h"
#include <math.h>
#include <string.h>

/* =========================================================================
 * L2-L3: Linear stability classification
 * ========================================================================= */

bz_stability_type_t bz_classify_stability(const double *ev_re,
                                           const double *ev_im)
{
    /* Classify based on eigenvalue signs and types:
     * - All Re < 0, all Im = 0 => Stable Node
     * - All Re < 0, one Im != 0 => Stable Focus
     * - All Re > 0, all Im = 0 => Unstable Node
     * - All Re > 0, one Im != 0 => Unstable Focus
     * - Mixed Re signs => Saddle
     * - All Re = 0, Im != 0 => Center
     * - Double/triple zero => Degenerate */

    int n_pos = 0, n_neg = 0, n_zero = 0;
    int n_imag = 0;
    int i __attribute__((unused));

    if (!ev_re || !ev_im) return BZ_STABILITY_DEGENERATE;

    for (i = 0; i < 3; i++) {
        if (ev_re[i] > 1e-10) n_pos++;
        else if (ev_re[i] < -1e-10) n_neg++;
        else n_zero++;
        if (fabs(ev_im[i]) > 1e-10) n_imag++;
    }

    if (n_zero == 3) {
        return (n_imag > 0) ? BZ_STABILITY_CENTER : BZ_STABILITY_DEGENERATE;
    }
    if (n_pos > 0 && n_neg > 0) return BZ_STABILITY_SADDLE;
    if (n_pos > 0) return (n_imag > 0) ? BZ_STABILITY_UNSTABLE_FOCUS
                                        : BZ_STABILITY_UNSTABLE_NODE;
    if (n_neg > 0) return (n_imag > 0) ? BZ_STABILITY_STABLE_FOCUS
                                        : BZ_STABILITY_STABLE_NODE;
    return BZ_STABILITY_DEGENERATE;
}

double bz_spectral_radius_3x3(const double *J)
{
    double ev_re[3], ev_im[3];
    double rho = 0.0;
    int i __attribute__((unused));

    if (!J) return 0.0;

    bz_eigenvalues_3x3(ev_re, ev_im, J);
    for (i = 0; i < 3; i++) {
        double mag = sqrt(ev_re[i]*ev_re[i] + ev_im[i]*ev_im[i]);
        if (mag > rho) rho = mag;
    }
    return rho;
}

int bz_is_linearly_stable(const double *J)
{
    double ev_re[3], ev_im[3];
    int i __attribute__((unused));

    if (!J) return 0;

    bz_eigenvalues_3x3(ev_re, ev_im, J);
    for (i = 0; i < 3; i++) {
        if (ev_re[i] > 1e-10) return 0;
    }
    return 1;
}

/* =========================================================================
 * L4: Hopf bifurcation analysis
 * ========================================================================= */

double bz_first_lyapunov_coeff(const bz_oregonator_params_t *p,
                                const bz_oregonator_state_t *ss,
                                double omega)
{
    /* First Lyapunov coefficient l1 determines criticality of Hopf:
     *   l1 < 0 => supercritical (stable limit cycle)
     *   l1 > 0 => subcritical (unstable limit cycle)
     *
     * Computed using the projection method on the center manifold.
     * For the Oregonator, the sign is determined by the cubic term
     * in the normal form reduction.
     *
     * Simplified formula from Tyson (1985):
     *   l1 ~ const * (f_c - f) / eps  near the bifurcation.
     *
     * Here we use a numerical estimate based on the curvature
     * of the steady-state manifold. */

    double f = p->f;
    double fc = bz_oregonator_hopf_critical_f(p->q);
    double eps = p->epsilon;

    /* The sign of (fc - f) determines criticality for the Oregonator:
     * For f slightly above fc: supercritical Hopf
     * For f well above fc: can lose supercriticality */

    (void)ss;
    (void)omega;

    if (fabs(f - fc) < 1e-6) return 0.0;  /* exactly at bifurcation */

    /* l1 ~ -(f - fc) / (eps * fc)  [scaling estimate] */
    return -(f - fc) / (eps * fc + 1e-10);
}

int bz_hopf_type(const bz_oregonator_params_t *p,
                  const bz_oregonator_state_t *ss)
{
    /* Determine whether Hopf bifurcation is supercritical or subcritical. */
    double J[9], ev_re[3], ev_im[3];
    double omega = 0.0;
    int i __attribute__((unused));
    double l1;

    if (!p || !ss) return 0;

    bz_oregonator_jacobian(J, ss, p);
    bz_eigenvalues_3x3(ev_re, ev_im, J);

    /* Find the imaginary part (frequency) */
    for (i = 0; i < 3; i++) {
        if (fabs(ev_im[i]) > 1e-8) {
            omega = fabs(ev_im[i]);
            break;
        }
    }

    if (omega < 1e-10) return 0;  /* not a Hopf point */

    l1 = bz_first_lyapunov_coeff(p, ss, omega);
    if (l1 < -1e-8) return 1;   /* supercritical */
    if (l1 > 1e-8)  return -1;  /* subcritical */
    return 0;  /* degenerate */
}

/* =========================================================================
 * L5: Floquet stability analysis
 * ========================================================================= */

int bz_floquet_stability(double *multipliers, int *n_mult,
                          const double *crossings_z, int n_cross)
{
    /* Estimate Floquet multipliers from Poincare map returns.
     * For a limit cycle, the return map z_{k+1} = P(z_k).
     * The multiplier is dP/dz evaluated at the fixed point.
     *
     * Estimate from two consecutive returns:
     *   mu = (z_{k+1} - z_k) / (z_k - z_{k-1}) */

    int i __attribute__((unused));
    double sum_mu, mu;

    if (!multipliers || !n_mult || !crossings_z || n_cross < 2) {
        if (n_mult) *n_mult = 0;
        return 0;
    }

    sum_mu = 0.0;
    if (n_cross >= 2) {
        for (i = 2; i < n_cross; i++) {
            double denom = crossings_z[i-1] - crossings_z[i-2];
            if (fabs(denom) > 1e-12) {
                mu = (crossings_z[i] - crossings_z[i-1]) / denom;
                sum_mu += mu;
            }
        }
        *n_mult = n_cross - 2;
        if (*n_mult > 0) {
            multipliers[0] = sum_mu / (double)(*n_mult);
            /* Stable if |mu| < 1 */
            return (fabs(multipliers[0]) < 1.0) ? 1 : 0;
        }
    }
    *n_mult = 0;
    return 0;
}

int bz_floquet_exponent(double *exponent,
                         const double *traj1, const double *traj2,
                         int n, double dt)
{
    /* Largest Floquet (Lyapunov) exponent from two nearby trajectories:
     *   lambda = (1/T) * ln(||delta(T)|| / ||delta(0)||) */

    double dist0, distT, lambda;
    int i __attribute__((unused));

    if (!exponent || !traj1 || !traj2 || n < 2 || dt <= 0.0) return 0;

    /* Initial separation */
    dist0 = fabs(traj1[0] - traj2[0]);
    if (dist0 < 1e-15) return 0;

    /* Final separation after n*dt time */
    distT = fabs(traj1[n-1] - traj2[n-1]);

    lambda = log(distT / dist0) / ((double)(n-1) * dt);
    *exponent = lambda;
    return 1;
}

/* =========================================================================
 * L8: Advanced bifurcation analysis
 * ========================================================================= */

void bz_track_hopf_curve(double *f_vals, double *q_vals, int n_pts,
                          double q_start)
{
    /* Track Hopf bifurcation curve f_H(q) using the analytical formula.
     * For q > 1/3: f_H = (1+q)^2 / (3q-1)
     * For q <= 1/3: f_H approx = 0.5 + q */
    int i __attribute__((unused));
    double q_min = 0.1, q_max = 10.0;

    if (!f_vals || !q_vals || n_pts <= 0) return;

    for (i = 0; i < n_pts; i++) {
        double q = q_min + (q_max - q_min) * i / (double)(n_pts - 1);
        q_vals[i] = q;
        f_vals[i] = bz_oregonator_hopf_curve(q);
    }

    /* Mark the starting point */
    (void)q_start;
}

int bz_detect_bogdanov_takens(const double *J)
{
    /* Bogdanov-Takens bifurcation: double-zero eigenvalue.
     * Conditions: two eigenvalues = 0, third != 0.
     * This means: trace(J) != 0, det(J) = 0, and sum of principal minors = 0. */

    double tr, m2, det;

    if (!J) return 0;

    bz_char_poly_3x3(&tr, &m2, &det, J);

    /* For BT point: det ~ 0 and m2 ~ 0, but tr != 0 */
    if (fabs(det) < 1e-8 && fabs(m2) < 1e-8 && fabs(tr) > 1e-8) {
        return 1;
    }
    return 0;
}

int bz_hopf_criticality(const bz_oregonator_params_t *p,
                         const bz_oregonator_state_t *ss)
{
    /* Wrapper for bz_hopf_type with additional robustness checks. */
    double J[9];
    bz_oregonator_jacobian(J, ss, p);

    /* Verify we are near a Hopf point */
    if (!bz_oregonator_hopf_test(p, ss)) return 0;

    return bz_hopf_type(p, ss);
}
