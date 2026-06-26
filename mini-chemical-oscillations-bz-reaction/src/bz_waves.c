
/**
 * @file bz_waves.c
 * @brief Chemical wave analysis: front detection, wave speed,
 *        eikonal equation, gradient, target patterns, spiral arms.
 *
 * Implements:
 *   L2: Wave front detection, wave speed computation
 *   L3: Eikonal equation (fast marching), front curvature
 *   L5: Gradient computation, cross-correlation, probe analysis
 *   L6: Target pattern generation, spiral arm counting
 *
 * References:
 *   - Winfree (1987) "When Time Breaks Down"
 *   - Keener & Sneyd (1998) "Mathematical Physiology"
 *   - Sethian (1999) "Level Set Methods and Fast Marching Methods"
 */

#include "bz_waves.h"
#include "bz_reaction_diffusion.h"
#include "bz_oscillation.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* =========================================================================
 * L2: Wave front detection
 * ========================================================================= */

int bz_wave_front_2d(int *front, const double *u, int nx, int ny,
                      double threshold)
{
    /* Detect wave front as points where u crosses threshold from below.
     * A front point satisfies: u[i] >= threshold AND neighbor u[prev] < threshold */
    int i, j, idx, count = 0;

    if (!front || !u || nx < 2 || ny < 2) return 0;

    memset(front, 0, nx * ny * sizeof(int));

    for (j = 0; j < ny; j++) {
        for (i = 1; i < nx; i++) {
            idx = j * nx + i;
            if (u[idx] >= threshold && u[idx-1] < threshold) {
                front[idx] = 1;
                count++;
            }
        }
    }
    return count;
}

double bz_wave_speed_field(double *speed, const double *u1, const double *u2,
                            int nx, int ny, double dx, double dy __attribute__((unused)), double dt)
{
    /* Compute wave speed via level-set method:
     *   c = (u2 - u1) / (dt * |grad u|)
     * Average speed is returned. */
    int i, j, idx, count = 0;
    double *grad_x, *grad_y;
    double total_speed = 0.0;
    double eps = 1e-8;

    if (!speed || !u1 || !u2 || nx < 2 || ny < 2 || dt <= 0.0) return 0.0;

    grad_x = (double*)malloc(nx * ny * sizeof(double));
    grad_y = (double*)malloc(nx * ny * sizeof(double));

    if (!grad_x || !grad_y) {
        free(grad_x);
        free(grad_y);
        return 0.0;
    }

    bz_gradient_2d(grad_x, grad_y, u1, nx, ny, dx, dy);

    for (j = 0; j < ny; j++) {
        for (i = 0; i < nx; i++) {
            idx = j * nx + i;
            double grad_mag = sqrt(grad_x[idx]*grad_x[idx]
                                    + grad_y[idx]*grad_y[idx]);

            if (grad_mag > eps) {
                double dtdu = (u2[idx] - u1[idx]) / dt;
                speed[idx] = -dtdu / grad_mag;
            } else {
                speed[idx] = 0.0;
            }

            if (fabs(speed[idx]) > eps) {
                total_speed += speed[idx];
                count++;
            }
        }
    }

    free(grad_x);
    free(grad_y);

    return (count > 0) ? (total_speed / (double)count) : 0.0;
}

/* =========================================================================
 * L3: Eikonal equation (fast marching method)
 * ========================================================================= */

void bz_eikonal_fast_marching(double *T, const double *speed,
                               int nx, int ny, double dx, double dy __attribute__((unused)),
                               int start_x, int start_y)
{
    /* Fast marching method for eikonal equation |grad T| = 1/c.
     * Simplified version: Dijkstra-like front propagation.
     *
     * For each neighbor of known points, solve the quadratic:
     *   max(T - T_left, 0)^2 + max(T - T_down, 0)^2 = (dx/c)^2 */

    int i, j, n;
    double *known;
    double h = dx;  /* assume dx = dy */

    if (!T || !speed || nx < 3 || ny < 3) return;

    n = nx * ny;

    /* Initialize */
    for (i = 0; i < n; i++) T[i] = 1e10;
    known = (double*)calloc(n, sizeof(double));
    if (!known) return;

    T[start_y * nx + start_x] = 0.0;
    known[start_y * nx + start_x] = 1.0;

    /* Simple iterative propagation (not true FMM, but captures the idea) */
    {
        int iter, changed;
        for (iter = 0; iter < 100; iter++) {
            changed = 0;
            for (j = 1; j < ny - 1; j++) {
                for (i = 1; i < nx - 1; i++) {
                    int idx = j * nx + i;
                    double T_left  = T[idx - 1];
                    double T_right = T[idx + 1];
                    double T_up    = T[idx - nx];
                    double T_down  = T[idx + nx];
                    double T_min_x = (T_left < T_right) ? T_left : T_right;
                    double T_min_y = (T_up < T_down)   ? T_up : T_down;
                    double c = speed[idx];
                    double T_new;

                    if (c < 1e-10) continue;

                    if (fabs(T_min_x - T_min_y) > 1e10) {
                        /* One-sided update */
                        double T_min = (T_min_x < T_min_y) ? T_min_x : T_min_y;
                        T_new = T_min + h / c;
                    } else {
                        /* Two-sided update */
                        double a = 2.0;
                        double b = -2.0 * (T_min_x + T_min_y);
                        double cc = T_min_x*T_min_x + T_min_y*T_min_y
                                    - (h*h)/(c*c);
                        double disc = b*b - 4.0*a*cc;
                        if (disc >= 0.0) {
                            T_new = (-b + sqrt(disc)) / (2.0*a);
                        } else {
                            T_new = T_min_x + h / c;
                        }
                    }

                    if (T_new < T[idx] - 1e-12) {
                        T[idx] = T_new;
                        changed = 1;
                    }
                }
            }
            if (!changed) break;
        }
    }

    free(known);
}

void bz_front_curvature(double *curvature, const double *u,
                         int nx, int ny, double dx, double dy)
{
    /* Curvature of level sets: kappa = div(grad u / |grad u|)
     * Computed via central differences. */
    int i, j, idx;
    double eps = 1e-12;

    if (!curvature || !u || nx < 3 || ny < 3) return;

    for (j = 1; j < ny - 1; j++) {
        for (i = 1; i < nx - 1; i++) {
            idx = j * nx + i;

            /* First derivatives */
            double ux = (u[idx+1] - u[idx-1]) / (2.0 * dx);
            double uy = (u[idx+nx] - u[idx-nx]) / (2.0 * dy);

            /* Second derivatives */
            double uxx = (u[idx+1] - 2.0*u[idx] + u[idx-1]) / (dx * dx);
            double uyy = (u[idx+nx] - 2.0*u[idx] + u[idx-nx]) / (dy * dy);
            double uxy = (u[idx+1+nx] - u[idx-1+nx] - u[idx+1-nx] + u[idx-1-nx])
                         / (4.0 * dx * dy);

            double grad_mag_sq = ux*ux + uy*uy + eps;

            /* kappa = (uxx*uy^2 - 2*ux*uy*uxy + uyy*ux^2) / |grad u|^3 */
            curvature[idx] = (uxx*uy*uy - 2.0*ux*uy*uxy + uyy*ux*ux)
                             / (grad_mag_sq * sqrt(grad_mag_sq));
        }
    }
}

/* =========================================================================
 * L5: Gradient and cross-correlation
 * ========================================================================= */

void bz_gradient_2d(double *grad_x, double *grad_y, const double *u,
                     int nx, int ny, double dx, double dy)
{
    /* Central difference gradient */
    int i, j, idx;

    if (!grad_x || !grad_y || !u || nx < 3 || ny < 3) return;

    for (j = 1; j < ny - 1; j++) {
        for (i = 1; i < nx - 1; i++) {
            idx = j * nx + i;
            grad_x[idx] = (u[idx+1] - u[idx-1]) / (2.0 * dx);
            grad_y[idx] = (u[idx+nx] - u[idx-nx]) / (2.0 * dy);
        }
    }
}

double bz_probe_period(const double *times, int n, double dt)
{
    /* Compute dominant period at a spatial probe point. */
    int peaks[100], n_peaks;
    double period;

    if (!times || n < 5) return 0.0;

    n_peaks = bz_count_peaks(times, n, peaks, 100);

    if (n_peaks < 2) return 0.0;

    bz_estimate_period(peaks, n_peaks, dt, &period, NULL);
    return period;
}

double bz_cross_correlation(int *delay,
                             const double *ts1, const double *ts2,
                             int n, int max_lag)
{
    /* Compute cross-correlation and find lag with maximum correlation.
     * C[lag] = sum_i ts1[i] * ts2[i+lag] */
    int lag, i, best_lag = 0;
    double max_corr = -1e10;

    if (!delay || !ts1 || !ts2 || n <= max_lag) {
        if (delay) *delay = 0;
        return 0.0;
    }

    for (lag = -max_lag; lag <= max_lag; lag++) {
        double corr = 0.0;
        int count = 0;
        for (i = 0; i < n; i++) {
            int j = i + lag;
            if (j >= 0 && j < n) {
                corr += ts1[i] * ts2[j];
                count++;
            }
        }
        if (count > 0) corr /= (double)count;
        if (corr > max_corr) {
            max_corr = corr;
            best_lag = lag;
        }
    }

    if (delay) *delay = best_lag;
    return max_corr;
}

/* =========================================================================
 * L6: Pattern generation
 * ========================================================================= */

void bz_pacemaker_target(double *u, double *v, int nx, int ny,
                          int cx, int cy, int radius,
                          int period, int step)
{
    /* Create target pattern from a periodic pacemaker.
     * Every 'period' steps, set a disc of u=1, v=0 around (cx,cy). */
    int i, j, idx;
    int r2 = radius * radius;

    if (!u || !v) return;

    if (step % period == 0) {
        for (j = -radius; j <= radius; j++) {
            for (i = -radius; i <= radius; i++) {
                if (i*i + j*j <= r2) {
                    int px = cx + i, py = cy + j;
                    if (px >= 0 && px < nx && py >= 0 && py < ny) {
                        idx = py * nx + px;
                        u[idx] = 1.0;
                        v[idx] = 0.0;
                    }
                }
            }
        }
    }
}

int bz_spiral_arm_count(double tip_x, double tip_y,
                         const double *u, int nx, int ny)
{
    /* Count spiral arms by computing winding number around the tip.
     * Winding number = (1/2*pi) * sum(delta_theta) along a closed loop.
     *
     * Simplified: count the number of phase = 0 -> 2*pi transitions
     * along a circle centered at the tip. */
    int n_samples = 36;  /* sample points around circle */
    double radius = 5.0;  /* circle radius in grid units */
    int k;
    double pi = 3.141592653589793;
    double prev_phase = 0.0;
    int arms = 0;

    if (!u) return 0;

    for (k = 0; k < n_samples; k++) {
        double angle = 2.0 * pi * k / n_samples;
        int sx = (int)(tip_x + radius * cos(angle));
        int sy = (int)(tip_y + radius * sin(angle));

        if (sx < 0 || sx >= nx || sy < 0 || sy >= ny) continue;

        {
            int idx = sy * nx + sx;
            double phase = atan2(u[idx] - 0.5, u[idx] - 0.3);

            if (k > 0) {
                double dphase = phase - prev_phase;
                /* Normalize to [-pi, pi] */
                while (dphase > pi)  dphase -= 2.0 * pi;
                while (dphase < -pi) dphase += 2.0 * pi;
                if (fabs(dphase) > pi * 0.8) arms++;
            }
            prev_phase = phase;
        }
    }

    return arms;
}
