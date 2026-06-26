/**
 * @file bz_bistability.h
 * @brief Bistability and multistability analysis in BZ-type reactions.
 *        Hysteresis, basin of attraction, saddle-node bifurcation.
 *
 * Knowledge coverage:
 *   L2: Bistability, hysteresis, multistability
 *   L3: Multiple steady states, saddle-node bifurcation
 *   L4: Bistability conditions, exchange of stability
 *   L5: Basin of attraction computation
 *   L7: Bistable chemical switch applications
 *
 * Bistability in BZ arises when two stable steady states coexist,
 * separated by an unstable saddle, enabling hysteresis.
 */

#ifndef BZ_BISTABILITY_H
#define BZ_BISTABILITY_H

#include "bz_reaction.h"
#include "bz_stability.h"

/* ---------------------------------------------------------------------------
 * L2: Bistability characterization
 * --------------------------------------------------------------------------- */

/** Hysteresis loop description for bistable BZ reaction. */
typedef struct {
    double  f_low;       /* Low-f transition (saddle-node 1) */
    double  f_high;      /* High-f transition (saddle-node 2) */
    double  width;       /* Hysteresis width: f_high - f_low */
    int     n_steady;    /* Number of steady states in bistable region */
} bz_hysteresis_t;

/** Detect if system is bistable for given Oregonator parameters.
 *  Bistable when bz_oregonator_steady_states returns 3 real solutions,
 *  two stable and one saddle.
 *  @param p    Oregonator parameters
 *  @return     1 if bistable, 0 otherwise */
int bz_is_bistable(const bz_oregonator_params_t *p);

/** Compute hysteresis loop boundaries for Oregonator.
 *  Finds f values where saddle-node bifurcations occur.
 *  @param hysteresis  Output hysteresis description
 *  @param q           Rate constant ratio (fixed)
 *  @return            1 if hysteresis found, 0 otherwise */
int bz_compute_hysteresis(bz_hysteresis_t *hysteresis, double q);

/* ---------------------------------------------------------------------------
 * L3: Multiple steady state analysis
 * --------------------------------------------------------------------------- */

/** Classify stability of all steady states.
 *  @param stability  Output stability type for each steady state
 *  @param max_ss     Maximum number of steady states
 *  @param ss         Steady state array
 *  @param n_ss       Number of steady states
 *  @param p          Parameters */
void bz_classify_steady_states(bz_stability_type_t *stability, int max_ss,
                                const bz_oregonator_state_t *ss, int n_ss,
                                const bz_oregonator_params_t *p);

/** Find saddle-node bifurcation points by bisection.
 *  At saddle-node, Jacobian has a zero eigenvalue (det = 0).
 *  @param f_sn   Output f value at saddle-node
 *  @param f_lo   Lower bound for bisection
 *  @param f_hi   Upper bound for bisection
 *  @param q      Fixed rate ratio
 *  @param tol    Tolerance for f
 *  @return       1 if found, 0 otherwise */
int bz_find_saddle_node(double *f_sn, double f_lo, double f_hi,
                         double q, double tol);

/* ---------------------------------------------------------------------------
 * L5: Basin of attraction computation
 * --------------------------------------------------------------------------- */

/** Determine which steady state an initial condition converges to.
 *  Integrates Oregonator and checks proximity to each steady state.
 *  @param s0     Initial condition
 *  @param p      Parameters (in bistable regime)
 *  @param ss     Array of steady states
 *  @param n_ss   Number of steady states
 *  @return       Index of attracting steady state (-1 if no convergence) */
int bz_basin_attractor(const bz_oregonator_state_t *s0,
                        const bz_oregonator_params_t *p,
                        const bz_oregonator_state_t *ss, int n_ss);

/** Compute basin of attraction via grid sampling.
 *  Samples a 2D slice of state space, maps each point to its attractor.
 *  @param basin      Output basin labels for each grid point (size nx*nz)
 *  @param x_range    [x_min, x_max]
 *  @param z_range    [z_min, z_max]
 *  @param nx, nz     Grid size
 *  @param y_fixed    Fixed y value (2D slice)
 *  @param p          Parameters
 *  @param ss         Steady states
 *  @param n_ss       Number of steady states */
void bz_basin_grid(int *basin,
                    double x_min, double x_max, int nx,
                    double z_min, double z_max, int nz,
                    double y_fixed,
                    const bz_oregonator_params_t *p,
                    const bz_oregonator_state_t *ss, int n_ss);

/* ---------------------------------------------------------------------------
 * L7: Applications: chemical switch and memory
 * --------------------------------------------------------------------------- */

/** Bistable chemical switch: toggle state by perturbation.
 *  Simulates a flip-flop using BZ bistability.
 *  @param state      Current state (modified)
 *  @param p          Parameters (in bistable regime)
 *  @param perturb_x  x perturbation magnitude
 *  @param target     0 = switch to lower state, 1 = switch to upper state
 *  @return           1 if switch successful */
int bz_chemical_switch(bz_oregonator_state_t *state,
                        const bz_oregonator_params_t *p,
                        double perturb_x, int target);

#endif /* BZ_BISTABILITY_H */
