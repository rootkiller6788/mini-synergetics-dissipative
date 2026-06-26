/**
 * @file bz_waves.h
 * @brief Chemical wave propagation in BZ reaction: spiral waves,
 *        target patterns, eikonal equation, wave speed.
 *
 * Knowledge coverage:
 *   L2: Excitable media, wave propagation, curvature effects
 *   L3: Eikonal equation, wave front kinematics
 *   L5: Wave speed computation, curvature-velocity relation
 *   L6: Spiral and target wave characterization
 *   L7: Chemical wave applications (cardiac modeling)
 *
 * References:
 *   - Winfree (1987) "When Time Breaks Down"
 *   - Keener & Sneyd (1998) "Mathematical Physiology"
 */

#ifndef BZ_WAVES_H
#define BZ_WAVES_H

#include "bz_reaction_diffusion.h"

/* ---------------------------------------------------------------------------
 * L2: Wave propagation concepts
 * --------------------------------------------------------------------------- */

/** Wave front detection by threshold crossing.
 *  Identifies points where u crosses threshold from below to above.
 *  @param front     Output binary mask (1=front, size nx*ny)
 *  @param u         Activator field
 *  @param nx, ny    Grid dimensions
 *  @param threshold Detection threshold
 *  @return          Number of front points found */
int bz_wave_front_2d(int *front, const double *u, int nx, int ny,
                      double threshold);

/** Compute wave propagation speed from two time-separated fields.
 *  Uses level-set method: c = (u2 - u1) / (dt * |grad u|)
 *  @param u1        Field at time t
 *  @param u2        Field at time t+dt
 *  @param nx, ny    Grid dimensions
 *  @param dx, dy    Grid spacing
 *  @param dt        Time between fields
 *  @param speed     Output speed field (size nx*ny)
 *  @return          Average wave speed */
double bz_wave_speed_field(double *speed, const double *u1, const double *u2,
                            int nx, int ny, double dx, double dy, double dt);

/* ---------------------------------------------------------------------------
 * L3: Eikonal equation for wave propagation
 * --------------------------------------------------------------------------- */

/** Solve eikonal equation |grad T| = 1/c for wave arrival time T.
 *  Uses fast marching method (first-order upwind).
 *  @param T         Output arrival time field (size nx*ny)
 *  @param speed     Wave speed field (size nx*ny)
 *  @param nx, ny    Grid dimensions
 *  @param dx, dy    Grid spacing
 *  @param start_x   x-index of wave source
 *  @param start_y   y-index of wave source */
void bz_eikonal_fast_marching(double *T, const double *speed,
                               int nx, int ny, double dx, double dy,
                               int start_x, int start_y);

/** Compute curvature of wave front.
 *  kappa = div(grad u / |grad u|)
 *  @param curvature Output curvature field (size nx*ny)
 *  @param u         Activator field
 *  @param nx, ny    Grid dimensions
 *  @param dx, dy    Grid spacing */
void bz_front_curvature(double *curvature, const double *u,
                         int nx, int ny, double dx, double dy);

/* ---------------------------------------------------------------------------
 * L5: Wave analysis algorithms
 * --------------------------------------------------------------------------- */

/** Compute gradient field of activator concentration.
 *  @param grad_x    Output x-gradient (size nx*ny)
 *  @param grad_y    Output y-gradient (size nx*ny)
 *  @param u         Input field
 *  @param nx, ny    Grid dimensions
 *  @param dx, dy    Grid spacing */
void bz_gradient_2d(double *grad_x, double *grad_y, const double *u,
                     int nx, int ny, double dx, double dy);

/** Compute wave period at a probe point from time series.
 *  @param times     Probe recording at one spatial point
 *  @param n         Number of time samples
 *  @param dt        Sampling interval
 *  @return          Dominant period (0 if no oscillation) */
double bz_probe_period(const double *times, int n, double dt);

/** Cross-correlation of two probe time series for wave delay estimation.
 *  @param delay     Output: lag with maximum correlation
 *  @param ts1, ts2  Time series from two probes
 *  @param n         Number of samples
 *  @param max_lag   Maximum lag to test
 *  @return          Maximum cross-correlation value */
double bz_cross_correlation(int *delay,
                             const double *ts1, const double *ts2,
                             int n, int max_lag);

/* ---------------------------------------------------------------------------
 * L6: Target pattern and spiral wave generation
 * --------------------------------------------------------------------------- */

/** Create a target pattern by setting a periodic pacemaker.
 *  @param u, v      Fields (modified in-place)
 *  @param nx, ny    Grid dimensions
 *  @param cx, cy    Pacemaker center
 *  @param radius    Pacemaker radius
 *  @param period    Oscillation period (in time steps)
 *  @param step      Current time step index */
void bz_pacemaker_target(double *u, double *v, int nx, int ny,
                          int cx, int cy, int radius,
                          int period, int step);

/** Count spiral wave arms by winding number around tip.
 *  @param tip_x, tip_y  Spiral tip coordinates
 *  @param u             Activator field
 *  @param nx, ny        Grid dimensions
 *  @return              Number of spiral arms */
int bz_spiral_arm_count(double tip_x, double tip_y,
                         const double *u, int nx, int ny);

#endif /* BZ_WAVES_H */
