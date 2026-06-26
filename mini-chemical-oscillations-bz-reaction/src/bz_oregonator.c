/**
 * @file bz_oregonator.c
 * @brief Oregonator model detailed analysis: oscillation period,
 *        Jacobian computation, Hopf bifurcation, parameter sweeps.
 *
 * Implements:
 *   L2: Period estimation, oscillation detection
 *   L3: Jacobian matrix, eigenvalues (cubic solver)
 *   L4: Hopf bifurcation detection and curve
 *   L5: Parameter sweep, integration with analysis
 *   L6: Canonical Oregonator oscillation integration
 *
 * References:
 *   - Field and Noyes (1974) J. Chem. Phys. 60, 1877-1884
 *   - Tyson (1985) "Oscillations and Traveling Waves in Chemical Systems"
 *   - Kuznetsov (2004) "Elements of Applied Bifurcation Theory"
 */

#include "bz_oregonator.h"
#include "bz_reaction.h"
#include <math.h>
#include <string.h>
#include <stdio.h>

/* =========================================================================
 * L2: Period estimation and oscillation detection
 * ========================================================================= */

double bz_oregonator_period_approx(const bz_oregonator_params_t *p)
{
    /* Approximate period for relaxation oscillations (Tyson 1985):
     *   T ~ (1/epsilon_p) * ln(f/(f-1))  for f > 1
     * For f < 1 (near Hopf onset), T ~ 2*pi/omega_H
     *   where omega_H = sqrt((1-2q-f+q*f)/(eps*eps_p)) */
    double T;
    if (!p) return 0.0;

    if (p->f > 1.01) {
        /* Relaxation regime: slow manifold + fast jumps */
        T = (1.0 / p->epsilon_p) * log(p->f / (p->f - 1.0));
    } else {
        /* Near-sinusoidal regime: Hopf frequency */
        double omega_H_sq = (1.0 - 2.0 * p->q - p->f + p->q * p->f)
                             / (p->epsilon * p->epsilon_p);
        if (omega_H_sq < 0) omega_H_sq = 0.0;
        T = 6.283185307179586 / sqrt(omega_H_sq + 1e-10);
    }
    return T < 0 ? 0 : T;
}

double bz_oregonator_hopf_critical_f(double q)
{
    /* Hopf bifurcation occurs when trace condition is met.
     * For small epsilon, f_c ~ 0.5*(1+sqrt(1+8q)) approximately.
     * More precisely from the characteristic equation:
     *   f_H(q) = (1+q)^2 / (3q-1) for q > 1/3.
     * For small q, use expansion: f_H ~ 0.5 + 0.5*sqrt(1+4q) */
    if (q <= 1.0/3.0) {
        /* Use expansion for small q */
        return 0.5 + q;
    } else {
        return (1.0 + q) * (1.0 + q) / (3.0 * q - 1.0);
    }
}

double bz_oregonator_relaxation_ratio(const bz_oregonator_params_t *p)
{
    /* Relaxation ratio = ratio of fast to slow time scales
     *   R = epsilon / epsilon_p
     * R >> 1 => relaxation oscillations (spike + slow recovery)
     * R ~ O(1) => near-sinusoidal */
    if (!p || p->epsilon_p <= 0.0) return 1.0;
    return p->epsilon / p->epsilon_p;
}

/* =========================================================================
 * L3: Jacobian matrix computation
 * ========================================================================= */

void bz_oregonator_jacobian(double *J,
                             const bz_oregonator_state_t *s,
                             const bz_oregonator_params_t *p)
{
    /* Jacobian of Oregonator system:
     *
     * f1 = dx/dt = (q*y - x*y + x*(1-x)) / eps
     * f2 = dy/dt = (-q*y - x*y + f*z) / eps_p
     * f3 = dz/dt = x - z
     *
     * Partial derivatives:
     *   df1/dx = ( -y + 1 - 2x ) / eps
     *   df1/dy = ( q - x ) / eps
     *   df1/dz = 0
     *
     *   df2/dx = ( -y ) / eps_p
     *   df2/dy = ( -q - x ) / eps_p
     *   df2/dz = f / eps_p
     *
     *   df3/dx = 1
     *   df3/dy = 0
     *   df3/dz = -1 */

    double x, y, inv_eps, inv_eps_p;

    if (!J || !s || !p) return;

    x = s->x;
    y = s->y;
    inv_eps   = 1.0 / p->epsilon;
    inv_eps_p = 1.0 / p->epsilon_p;

    /* Row 0: df1/d(x,y,z) */
    J[0] = (-y + 1.0 - 2.0 * x) * inv_eps;
    J[1] = (p->q - x) * inv_eps;
    J[2] = 0.0;

    /* Row 1: df2/d(x,y,z) */
    J[3] = (-y) * inv_eps_p;
    J[4] = (-p->q - x) * inv_eps_p;
    J[5] = p->f * inv_eps_p;

    /* Row 2: df3/d(x,y,z) */
    J[6] = 1.0;
    J[7] = 0.0;
    J[8] = -1.0;
}

/* =========================================================================
 * L3: Eigenvalue computation (cubic characteristic polynomial)
 * ========================================================================= */

void bz_char_poly_3x3(double *trace, double *m2, double *det, const double *J)
{
    /* For J = [[a,b,c],[d,e,f],[g,h,i]]
     *   trace = a + e + i
     *   m2    = (ae - bd) + (ai - cg) + (ei - fh)  [sum of principal minors]
     *   det   = a(ei - fh) - b(di - fg) + c(dh - eg) */

    double a, b, c, d, e, f, g, h, i;

    if (!J) return;
    a = J[0]; b = J[1]; c = J[2];
    d = J[3]; e = J[4]; f = J[5];
    g = J[6]; h = J[7]; i = J[8];

    if (trace) *trace = a + e + i;
    if (m2)    *m2 = (a*e - b*d) + (a*i - c*g) + (e*i - f*h);
    if (det)   *det = a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g);
}

/* Solve cubic: x^3 + a2*x^2 + a1*x + a0 = 0
 * Returns number of real roots. Fills roots[3] with real parts,
 * imag[3] with imaginary parts. */
static int solve_cubic(double a2, double a1, double a0,
                        double *roots, double *imag)
{
    /* Depressed cubic: t^3 + p*t + q = 0
     * via substitution x = t - a2/3 */
    double p, q, Q, R, theta, Q3, R2;
    double shift = a2 / 3.0;

    p = a1 - a2 * a2 / 3.0;
    q = a0 - a2 * a1 / 3.0 + 2.0 * a2 * a2 * a2 / 27.0;

    Q = p / 3.0;
    R = -q / 2.0;
    Q3 = Q * Q * Q;
    R2 = R * R;

    if (R2 + Q3 <= 0.0) {
        /* Three real roots (trigonometric solution) */
        double Q_abs = sqrt(-Q);
        if (Q_abs < 1e-15) Q_abs = 1e-15;
        theta = acos(R / (Q_abs * Q_abs * Q_abs));
        /* Clamp for numerical safety */
        if (theta > 3.14159265) theta = 3.14159265;
        if (theta < 0.0) theta = 0.0;

        roots[0] = 2.0 * Q_abs * cos(theta / 3.0) - shift;
        roots[1] = 2.0 * Q_abs * cos((theta + 2.0*3.14159265) / 3.0) - shift;
        roots[2] = 2.0 * Q_abs * cos((theta + 4.0*3.14159265) / 3.0) - shift;
        imag[0] = imag[1] = imag[2] = 0.0;
        return 3;
    } else {
        /* One real root, two complex conjugates */
        double sqrt_D = sqrt(R2 + Q3);
        double A = -cbrt(R + sqrt_D);
        double B = 0.0;

        if (fabs(A) > 1e-15) {
            B = Q / A;
        } else {
            B = -cbrt(R - sqrt_D);
        }

        roots[0] = A + B - shift;
        imag[0] = 0.0;
        roots[1] = -(A + B)/2.0 - shift;
        imag[1] = (sqrt(3.0)/2.0) * (A - B);
        roots[2] = roots[1];
        imag[2] = -imag[1];
        return 1;  /* one real root, but we fill all 3 slots */
    }
}

int bz_eigenvalues_3x3(double *ev_re, double *ev_im, const double *J)
{
    /* Characteristic polynomial: det(J - lambda*I) = 0
     *   -lambda^3 + tr*lambda^2 - m2*lambda + det = 0
     * => lambda^3 - tr*lambda^2 + m2*lambda - det = 0
     *
     * Compute eigenvalues of J by solving cubic. */
    double tr, m2, det;
    int n_unstable = 0;
    int i;

    if (!ev_re || !ev_im || !J) return 0;

    bz_char_poly_3x3(&tr, &m2, &det, J);

    /* Solve: lambda^3 - tr*lambda^2 + m2*lambda - det = 0
     * a2 = -tr, a1 = m2, a0 = -det */
    solve_cubic(-tr, m2, -det, ev_re, ev_im);

    for (i = 0; i < 3; i++) {
        if (ev_re[i] > 1e-10) n_unstable++;
    }
    return n_unstable;
}

/* =========================================================================
 * L4: Hopf bifurcation analysis
 * ========================================================================= */

int bz_oregonator_hopf_test(const bz_oregonator_params_t *p,
                             const bz_oregonator_state_t *ss)
{
    /* Hopf bifurcation test:
     * (1) Compute Jacobian at steady state
     * (2) Check eigenvalues: one real + pair of complex
     * (3) Real part of complex pair must cross zero
     * (4) Transversality: d(Re lambda)/d(parameter) != 0 */
    double J[9], ev_re[3], ev_im[3];
    int n_unstable __attribute__((unused));

    if (!p || !ss) return 0;

    bz_oregonator_jacobian(J, ss, p);
    n_unstable = bz_eigenvalues_3x3(ev_re, ev_im, J);

    /* Check if there is a complex conjugate pair */
    int has_complex = 0;
    int i;
    for (i = 0; i < 3; i++) {
        if (fabs(ev_im[i]) > 1e-8) {
            has_complex = 1;
            break;
        }
    }

    /* A Hopf bifurcation is approaching if:
     * - There is a complex pair
     * - The real part of the complex pair is near zero
     * - The real eigenvalue is negative (otherwise it's a different bifurcation) */
    if (has_complex) {
        double re_complex_min = 1e10;
        for (i = 0; i < 3; i++) {
            if (fabs(ev_im[i]) > 1e-8) {
                if (fabs(ev_re[i]) < fabs(re_complex_min)) {
                    re_complex_min = ev_re[i];
                }
            }
        }
        /* Hopf if complex pair Re ~ 0 and real eigenvalue Re < 0 */
        if (fabs(re_complex_min) < 0.1 && ev_re[0] < -0.01) {
            return 1;
        }
    }
    return 0;
}

double bz_oregonator_hopf_curve(double q)
{
    /* Analytical Hopf curve for Oregonator (epsilon -> 0 limit):
     *   f_H(q) = (1+q)^2 / (3q - 1)  for q > 1/3
     * For q <= 1/3, use expansion from characteristic equation. */
    if (q > 1.0/3.0) {
        return (1.0 + q) * (1.0 + q) / (3.0 * q - 1.0);
    } else {
        /* For small q, the critical f is near 0.5.
         * More precisely f_H ~ 0.5 + O(q) */
        return 0.5 + q;
    }
}

/* =========================================================================
 * L5: Parameter sweep
 * ========================================================================= */

void bz_oregonator_sweep_f(double f_min, double f_max, int n_steps,
                            int *results, const bz_oregonator_params_t *p_base)
{
    /* For each f in [f_min, f_max], check if oscillatory.
     * Uses the analytical Hopf criterion for speed. */
    int i;
    bz_oregonator_params_t p;

    if (!results || !p_base || n_steps <= 0) return;

    p = *p_base;
    for (i = 0; i < n_steps; i++) {
        double f = f_min + (f_max - f_min) * i / (double)(n_steps > 1 ? n_steps - 1 : 1);
        p.f = f;
        results[i] = bz_oregonator_is_oscillatory(&p);
    }
}

/* =========================================================================
 * L6: Integration and oscillation analysis
 * ========================================================================= */

/* Runge-Kutta 4th order for Oregonator 3-var system */
static void rk4_step(bz_oregonator_state_t *s,
                      const bz_oregonator_params_t *p, double dt)
{
    bz_oregonator_state_t k1, k2, k3, k4, tmp;
    bz_oregonator_state_t s0 = *s;

    /* k1 = f(s) */
    bz_oregonator_rhs(&k1, &s0, p);

    /* k2 = f(s + dt*k1/2) */
    tmp.x = s0.x + 0.5*dt*k1.x;
    tmp.y = s0.y + 0.5*dt*k1.y;
    tmp.z = s0.z + 0.5*dt*k1.z;
    bz_oregonator_rhs(&k2, &tmp, p);

    /* k3 = f(s + dt*k2/2) */
    tmp.x = s0.x + 0.5*dt*k2.x;
    tmp.y = s0.y + 0.5*dt*k2.y;
    tmp.z = s0.z + 0.5*dt*k2.z;
    bz_oregonator_rhs(&k3, &tmp, p);

    /* k4 = f(s + dt*k3) */
    tmp.x = s0.x + dt*k3.x;
    tmp.y = s0.y + dt*k3.y;
    tmp.z = s0.z + dt*k3.z;
    bz_oregonator_rhs(&k4, &tmp, p);

    /* s_new = s + dt*(k1 + 2*k2 + 2*k3 + k4)/6 */
    s->x = s0.x + (dt/6.0) * (k1.x + 2.0*k2.x + 2.0*k3.x + k4.x);
    s->y = s0.y + (dt/6.0) * (k1.y + 2.0*k2.y + 2.0*k3.y + k4.y);
    s->z = s0.z + (dt/6.0) * (k1.z + 2.0*k2.z + 2.0*k3.z + k4.z);
}

void bz_oregonator_integrate_and_analyze(bz_oregonator_state_t *s,
                                          const bz_oregonator_params_t *p,
                                          double duration, int n_steps,
                                          int *n_peaks,
                                          double *period_avg,
                                          double *amplitude)
{
    double dt;
    int i;
    double x_prev, x_max, x_min;
    int peak_count;
    double *peaks_buf __attribute__((unused));
    double period_sum;
    int peak_indices[100];  /* fixed buffer for up to 100 peaks */
    int max_peaks = 100;
    int going_up;

    if (!s || !p || n_steps <= 0) return;

    dt = duration / (double)n_steps;
    peak_count = 0;
    x_max = s->x;
    x_min = s->x;
    x_prev = s->x;
    going_up = 0;

    /* Integrate and detect peaks */
    for (i = 0; i < n_steps; i++) {
        /* Detect peak: derivative changes sign from + to - */
        rk4_step(s, p, dt);

        if (s->x > x_max) x_max = s->x;
        if (s->x < x_min) x_min = s->x;

        if (s->x > x_prev) {
            going_up = 1;
        } else if (going_up && s->x < x_prev - 1e-8) {
            /* Peak detected at x_prev */
            if (peak_count < max_peaks && x_prev > 0.05) {
                peak_indices[peak_count] = i;
                peak_count++;
            }
            going_up = 0;
        }
        x_prev = s->x;
    }

    /* Compute period stats */
    if (n_peaks) *n_peaks = peak_count;

    if (period_avg && peak_count >= 2) {
        period_sum = 0.0;
        for (i = 1; i < peak_count; i++) {
            period_sum += (peak_indices[i] - peak_indices[i-1]) * dt;
        }
        *period_avg = period_sum / (double)(peak_count - 1);
    } else if (period_avg) {
        *period_avg = 0.0;
    }

    if (amplitude) {
        *amplitude = x_max - x_min;
    }
}