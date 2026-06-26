#include "pds_brusselator.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ============================================================================
 * Brusselator Model Implementation
 * Prigogine-Lefever (1968) trimolecular reaction model
 * ============================================================================ */

void pds_brusselator_rhs(const PDSBrusselatorState* state,
                          const PDSBrusselatorParams* params,
                          PDSBrusselatorState* deriv) {
    if (!state || !params || !deriv) return;

    double x = state->x;
    double y = state->y;
    double a = params->a;
    double b = params->b;

    /* dx/dt = a - (b+1)*x + x^2 * y */
    deriv->x = a - (b + 1.0) * x + x * x * y;

    /* dy/dt = b*x - x^2 * y */
    deriv->y = b * x - x * x * y;
}

void pds_brusselator_steady_state(const PDSBrusselatorParams* params,
                                   double* x_s, double* y_s) {
    if (!params || !x_s || !y_s) return;

    /* From d/dt = 0:
     * 0 = a - (b+1)*x_s + x_s^2 * y_s
     * 0 = b*x_s - x_s^2 * y_s
     * => x_s^2 * y_s = b*x_s => y_s = b/x_s (for x_s != 0)
     * => 0 = a - (b+1)*x_s + x_s^2 * (b/x_s) = a - (b+1)*x_s + b*x_s
     * => a - x_s = 0 => x_s = a
     * => y_s = b / a
     */
    *x_s = params->a;
    if (fabs(params->a) > 1e-15) {
        *y_s = params->b / params->a;
    } else {
        *y_s = 0.0;
    }
}

void pds_brusselator_jacobian(const PDSBrusselatorState* state,
                               const PDSBrusselatorParams* params,
                               double J[2][2]) {
    if (!state || !params || !J) return;

    double x = state->x;
    double y = state->y;

    /* J = [ df1/dx,  df1/dy  ]
     *     [ df2/dx,  df2/dy  ]
     *
     * f1 = a - (b+1)*x + x^2 * y
     *   df1/dx = -(b+1) + 2*x*y
     *   df1/dy = x^2
     *
     * f2 = b*x - x^2 * y
     *   df2/dx = b - 2*x*y
     *   df2/dy = -x^2
     */

    J[0][0] = -(params->b + 1.0) + 2.0 * x * y;
    J[0][1] = x * x;
    J[1][0] = params->b - 2.0 * x * y;
    J[1][1] = -(x * x);
}

double pds_brusselator_bifurcation_point(const PDSBrusselatorParams* params) {
    if (!params) return 0.0;

    /* At steady state (x_s=a, y_s=b/a):
     * J_ss = [ b-1,  a^2  ]
     *        [ -b,   -a^2 ]
     *
     * Trace(J_ss) = b - 1 - a^2
     * Det(J_ss) = -(b-1)*a^2 + b*a^2 = -b*a^2 + a^2 + b*a^2 = a^2 > 0
     *
     * Hopf bifurcation when Trace = 0:
     * b_c = 1 + a^2
     *
     * For b < b_c: Trace < 0 => stable focus
     * For b > b_c: Trace > 0 => unstable focus
     * For b = b_c: Hopf bifurcation (purely imaginary eigenvalues)
     */
    return 1.0 + params->a * params->a;
}

/* ============================================================================
 * RK4 Integration ? L5 Numerical Algorithm
 *
 * The classic 4th-order Runge-Kutta method for ODE integration.
 * For the Brusselator, this is the workhorse for simulating
 * limit cycle oscillations and bifurcation sequences.
 *
 * Butcher tableau for RK4:
 *   0   |
 *   1/2 | 1/2
 *   1/2 |  0   1/2
 *   1   |  0    0    1
 *  -------------------------
 *       | 1/6  1/3  1/3  1/6
 *
 * Global error: O(h^4)
 * ============================================================================ */

void pds_brusselator_rk4_step(PDSBrusselatorState* state,
                               const PDSBrusselatorParams* params,
                               double dt) {
    if (!state || !params || dt <= 0.0) return;

    PDSBrusselatorState k1, k2, k3, k4, temp;

    /* k1 = f(state) */
    pds_brusselator_rhs(state, params, &k1);

    /* k2 = f(state + dt/2 * k1) */
    temp.x = state->x + 0.5 * dt * k1.x;
    temp.y = state->y + 0.5 * dt * k1.y;
    pds_brusselator_rhs(&temp, params, &k2);

    /* k3 = f(state + dt/2 * k2) */
    temp.x = state->x + 0.5 * dt * k2.x;
    temp.y = state->y + 0.5 * dt * k2.y;
    pds_brusselator_rhs(&temp, params, &k3);

    /* k4 = f(state + dt * k3) */
    temp.x = state->x + dt * k3.x;
    temp.y = state->y + dt * k3.y;
    pds_brusselator_rhs(&temp, params, &k4);

    /* Update: state += dt/6 * (k1 + 2k2 + 2k3 + k4) */
    state->x += (dt / 6.0) * (k1.x + 2.0 * k2.x + 2.0 * k3.x + k4.x);
    state->y += (dt / 6.0) * (k1.y + 2.0 * k2.y + 2.0 * k3.y + k4.y);

    /* Ensure concentrations stay non-negative */
    if (state->x < 0.0) state->x = 0.0;
    if (state->y < 0.0) state->y = 0.0;
}

PDSBrusselatorResult* pds_brusselator_simulate(
    const PDSBrusselatorState* initial,
    const PDSBrusselatorParams* params,
    double duration, double dt,
    int save_every) {
    if (!initial || !params || duration <= 0.0 || dt <= 0.0 || save_every <= 0) {
        return NULL;
    }

    PDSBrusselatorResult* result =
        (PDSBrusselatorResult*)calloc(1, sizeof(PDSBrusselatorResult));
    if (!result) return NULL;

    int n_steps = (int)(duration / dt) + 1;
    int n_saved = n_steps / save_every + 2;

    result->n_points = n_saved;
    result->trajectory = (PDSBrusselatorState*)calloc(n_saved,
        sizeof(PDSBrusselatorState));
    result->time = (double*)calloc(n_saved, sizeof(double));

    if (!result->trajectory || !result->time) {
        pds_brusselator_result_free(result);
        return NULL;
    }

    /* Compute steady state and bifurcation point */
    pds_brusselator_steady_state(params, &result->steady_x, &result->steady_y);
    result->bifurcation_parameter_b_c = pds_brusselator_bifurcation_point(params);

    /* Check if beyond Hopf bifurcation */
    if (params->b > result->bifurcation_parameter_b_c) {
        result->is_oscillatory = true; /* Expect limit cycle */
    }

    /* Initialize simulation */
    PDSBrusselatorState state;
    state.x = initial->x;
    state.y = initial->y;

    /* Save initial condition */
    result->trajectory[0] = state;
    result->time[0] = 0.0;

    int save_idx = 1;
    double t = 0.0;

    /* Main integration loop */
    for (int step = 1; step <= n_steps && save_idx < n_saved; step++) {
        pds_brusselator_rk4_step(&state, params, dt);
        t = step * dt;

        if (step % save_every == 0) {
            result->trajectory[save_idx] = state;
            result->time[save_idx] = t;
            save_idx++;
        }
    }

    result->n_points = save_idx;

    /* Detect limit cycle if oscillatory */
    if (result->is_oscillatory) {
        pds_brusselator_detect_limit_cycle(result);
    }

    return result;
}

void pds_brusselator_result_free(PDSBrusselatorResult* result) {
    if (!result) return;
    free(result->trajectory);
    free(result->time);
    free(result);
}

/* ============================================================================
 * Limit Cycle Detection ? L6 Classical Problem
 *
 * Detects sustained chemical oscillations in the Brusselator.
 * Method:
 *   1. Wait for transient decay (skip first ~20% of trajectory)
 *   2. Find zero-crossings of x - x_s to determine period
 *   3. Compute amplitude from maxima/minima
 *   4. Check if successive periods are consistent (closed orbit)
 *
 * Theorem: For b > b_c = 1 + a^2, the Brusselator exhibits a stable
 * limit cycle. The period T ~ 2*pi / omega_0 where omega_0 = sqrt(a)
 * at the Hopf bifurcation point (from J_ss eigenvalues).
 * ============================================================================ */

void pds_brusselator_detect_limit_cycle(PDSBrusselatorResult* result) {
    if (!result || result->n_points < 10) return;

    /* Skip initial transient: analyze last 60% of trajectory */
    int start_idx = result->n_points * 40 / 100;
    if (start_idx < 0) start_idx = 0;
    int n_analyze = result->n_points - start_idx;
    if (n_analyze < 5) return;

    /* Find extrema of x(t) */
    double x_min = result->trajectory[start_idx].x;
    double x_max = x_min;
    int n_crossings = 0;
    double* crossing_times = (double*)malloc(n_analyze * sizeof(double));

    if (!crossing_times) return;

    for (int i = start_idx + 1; i < result->n_points; i++) {
        double x = result->trajectory[i].x;
        if (x < x_min) x_min = x;
        if (x > x_max) x_max = x;

        /* Detect zero-crossing of (x - x_s): upward crossing */
        double prev = result->trajectory[i - 1].x - result->steady_x;
        double curr = x - result->steady_x;
        if (prev < 0.0 && curr >= 0.0) {
            crossing_times[n_crossings] = result->time[i];
            n_crossings++;
        }
    }

    result->oscillation_amplitude = (x_max - x_min) / 2.0;

    /* Estimate period from crossing intervals */
    if (n_crossings >= 3) {
        double sum_periods = 0.0;
        int n_periods = 0;
        for (int i = 1; i < n_crossings; i++) {
            double period = crossing_times[i] - crossing_times[i - 1];
            if (period > 0.0) {
                sum_periods += period;
                n_periods++;
            }
        }
        if (n_periods > 0) {
            result->oscillation_period = sum_periods / n_periods;
            /* Check if periods are consistent (limit cycle = closed orbit) */
            double variance = 0.0;
            for (int i = 1; i < n_crossings; i++) {
                double period = crossing_times[i] - crossing_times[i - 1];
                if (period > 0.0) {
                    double diff = period - result->oscillation_period;
                    variance += diff * diff;
                }
            }
            variance /= n_periods;
            result->has_limit_cycle = (variance < 0.01 * result->oscillation_period *
                                        result->oscillation_period);
        }
    }

    /* Simple Lyapunov exponent estimate from trajectory divergence.
     * For a limit cycle, the maximum Lyapunov exponent should be ~0
     * (neutral stability along the cycle). */
    if (result->n_points >= 4) {
        double d0 = sqrt(
            pow(result->trajectory[start_idx].x - result->trajectory[start_idx + 1].x, 2) +
            pow(result->trajectory[start_idx].y - result->trajectory[start_idx + 1].y, 2));
        double d1 = sqrt(
            pow(result->trajectory[result->n_points - 2].x -
                result->trajectory[result->n_points - 1].x, 2) +
            pow(result->trajectory[result->n_points - 2].y -
                result->trajectory[result->n_points - 1].y, 2));
        if (d0 > 1e-10 && d1 > 1e-10) {
            double dt_total = result->time[result->n_points - 1] -
                              result->time[start_idx];
            if (dt_total > 0.0) {
                result->max_lyapunov_exponent = log(d1 / d0) / dt_total;
            }
        }
    }

    free(crossing_times);
    result->is_oscillatory = result->has_limit_cycle;
}

/* ============================================================================
 * Spatial Brusselator ? Reaction-Diffusion Extension
 *
 * L5: Extends the well-mixed Brusselator to spatially extended systems
 * with diffusion. The reaction-diffusion equations are:
 *
 *   dX/dt = a - (b+1)*X + X^2*Y + Dx * nabla^2 X
 *   dY/dt = b*X - X^2*Y + Dy * nabla^2 Y
 *
 * Discretization: finite differences on a 2D grid.
 * Laplacian: 5-point stencil
 *   nabla^2 u_{i,j} = (u_{i+1,j} + u_{i-1,j} + u_{i,j+1} + u_{i,j-1}
 *                       - 4*u_{i,j}) / dx^2
 *
 * Turing Instability Condition (L4):
 * For Dx << Dy, the homogeneous steady state can become unstable to
 * spatial perturbations with a characteristic wavenumber k_c.
 * This produces stationary spatial patterns ? Turing's morphogenesis.
 * ============================================================================ */

PDSBrusselatorSpatial* pds_brusselator_spatial_create(int nx, int ny,
                                                        double Dx, double Dy) {
    if (nx <= 0 || ny <= 0) return NULL;

    PDSBrusselatorSpatial* spat = (PDSBrusselatorSpatial*)calloc(1,
        sizeof(PDSBrusselatorSpatial));
    if (!spat) return NULL;

    spat->n_x = nx;
    spat->n_y = ny;
    spat->dx = 1.0 / (double)((nx < ny) ? nx : ny); /* Grid spacing */
    spat->Dx = Dx;
    spat->Dy = Dy;
    spat->has_turing_pattern = false;
    spat->dominant_wavenumber = 0.0;

    spat->spatial_profile = (PDSBrusselatorState*)calloc(nx * ny,
        sizeof(PDSBrusselatorState));
    if (!spat->spatial_profile) {
        free(spat);
        return NULL;
    }

    return spat;
}

void pds_brusselator_spatial_free(PDSBrusselatorSpatial* spat) {
    if (!spat) return;
    free(spat->spatial_profile);
    free(spat);
}

/* Helper: get index in flattened 2D array */
static int _idx(int i, int j, int ny) {
    return i * ny + j;
}

void pds_brusselator_spatial_rhs(const PDSBrusselatorSpatial* spat,
                                  const PDSBrusselatorParams* params,
                                  PDSBrusselatorSpatial* deriv) {
    if (!spat || !params || !deriv) return;

    int nx = spat->n_x;
    int ny = spat->n_y;
    double dx2 = spat->dx * spat->dx;
    double a = params->a;
    double b = params->b;

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int idx = _idx(i, j, ny);
            double X = spat->spatial_profile[idx].x;
            double Y = spat->spatial_profile[idx].y;

            /* Reaction terms */
            double rx = a - (b + 1.0) * X + X * X * Y;
            double ry = b * X - X * X * Y;

            /* Diffusion terms (5-point Laplacian) */
            double lapX = 0.0, lapY = 0.0;

            double X_c = X;
            double Y_c = Y;

            /* i-1 */
            if (i > 0) {
                double Xm = spat->spatial_profile[_idx(i-1, j, ny)].x;
                double Ym = spat->spatial_profile[_idx(i-1, j, ny)].y;
                lapX += Xm - X_c;
                lapY += Ym - Y_c;
            }
            /* i+1 */
            if (i < nx - 1) {
                double Xp = spat->spatial_profile[_idx(i+1, j, ny)].x;
                double Yp = spat->spatial_profile[_idx(i+1, j, ny)].y;
                lapX += Xp - X_c;
                lapY += Yp - Y_c;
            }
            /* j-1 */
            if (j > 0) {
                double Xm = spat->spatial_profile[_idx(i, j-1, ny)].x;
                double Ym = spat->spatial_profile[_idx(i, j-1, ny)].y;
                lapX += Xm - X_c;
                lapY += Ym - Y_c;
            }
            /* j+1 */
            if (j < ny - 1) {
                double Xp = spat->spatial_profile[_idx(i, j+1, ny)].x;
                double Yp = spat->spatial_profile[_idx(i, j+1, ny)].y;
                lapX += Xp - X_c;
                lapY += Yp - Y_c;
            }

            lapX /= dx2;
            lapY /= dx2;

            deriv->spatial_profile[idx].x = rx + spat->Dx * lapX;
            deriv->spatial_profile[idx].y = ry + spat->Dy * lapY;
        }
    }
}

void pds_brusselator_euler_step_spatial(PDSBrusselatorSpatial* spat,
                                         const PDSBrusselatorParams* params,
                                         double dt) {
    if (!spat || !params || dt <= 0.0) return;

    /* Allocate temporary derivative storage */
    PDSBrusselatorSpatial* deriv = pds_brusselator_spatial_create(
        spat->n_x, spat->n_y, spat->Dx, spat->Dy);
    if (!deriv) return;

    pds_brusselator_spatial_rhs(spat, params, deriv);

    int n_total = spat->n_x * spat->n_y;
    for (int i = 0; i < n_total; i++) {
        spat->spatial_profile[i].x += dt * deriv->spatial_profile[i].x;
        spat->spatial_profile[i].y += dt * deriv->spatial_profile[i].y;
        /* Non-negative constraint */
        if (spat->spatial_profile[i].x < 0.0) spat->spatial_profile[i].x = 0.0;
        if (spat->spatial_profile[i].y < 0.0) spat->spatial_profile[i].y = 0.0;
    }

    pds_brusselator_spatial_free(deriv);
}

/* ============================================================================
 * Display Functions
 * ============================================================================ */

void pds_brusselator_print_result(const PDSBrusselatorResult* result) {
    if (!result) { printf("BrusselatorResult: NULL\n"); return; }

    printf("=== Brusselator Simulation Results ===\n");
    printf("  Steady State: (x_s=%.4f, y_s=%.4f)\n",
           result->steady_x, result->steady_y);
    printf("  Bifurcation b_c = 1 + a^2 = %.4f\n",
           result->bifurcation_parameter_b_c);
    printf("  Trajectory points: %d\n", result->n_points);

    if (result->is_oscillatory) {
        printf("  OSCILLATORY (Dissipative Structure):\n");
        printf("    Period: %.4f\n", result->oscillation_period);
        printf("    Amplitude: %.4f\n", result->oscillation_amplitude);
        printf("    Limit Cycle: %s\n",
               result->has_limit_cycle ? "YES" : "not confirmed");
        printf("    Max Lyapunov Exponent: %.6f\n",
               result->max_lyapunov_exponent);
    } else {
        printf("  STABLE (Thermodynamic Branch)\n");
    }

    /* Print first and last few points */
    printf("  Trajectory preview:\n");
    int show = (result->n_points < 8) ? result->n_points : 8;
    for (int i = 0; i < show; i++) {
        printf("    t=%.2f: (x=%.4f, y=%.4f)\n",
               result->time[i],
               result->trajectory[i].x,
               result->trajectory[i].y);
    }
    if (result->n_points > 8) {
        printf("    ... (%d more points)\n", result->n_points - 8);
    }
}