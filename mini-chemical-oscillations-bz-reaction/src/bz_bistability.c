
/**
 * @file bz_bistability.c
 * @brief Bistability analysis: hysteresis, saddle-node detection,
 *        basin of attraction, chemical switch.
 *
 * Implements:
 *   L2: Bistability detection, hysteresis computation
 *   L3: Stability classification of multiple steady states
 *   L5: Saddle-node bifurcation finding via bisection
 *   L5: Basin of attraction sampling
 *   L7: Chemical switch (bistable flip-flop)
 *
 * References:
 *   - Epstein & Pojman (1998) "Nonlinear Chemical Dynamics"
 *   - Strogatz (2015) "Nonlinear Dynamics and Chaos"
 */

#include "bz_bistability.h"
#include "bz_reaction.h"
#include "bz_stability.h"
#include "bz_oregonator.h"
#include <math.h>
#include <string.h>

/* =========================================================================
 * L2: Bistability detection
 * ========================================================================= */

int bz_is_bistable(const bz_oregonator_params_t *p)
{
    /* Check if the Oregonator has 3 steady states (bistable).
     * This occurs when the quadratic for x has 3 real positive roots.
     *
     * Condition: discriminant of the cubic steady-state equation > 0
     * and all three roots are real and positive. */

    bz_oregonator_state_t ss[4];
    int n_ss, n_stable;
    int i;

    if (!p) return 0;

    n_ss = bz_oregonator_steady_states(ss, p);
    if (n_ss < 3) return 0;

    /* Check stability: need two stable and one saddle for bistability */
    n_stable = 0;
    for (i = 0; i < n_ss; i++) {
        double J[9];
        bz_oregonator_jacobian(J, &ss[i], p);
        if (bz_is_linearly_stable(J)) n_stable++;
    }

    return (n_stable >= 2) ? 1 : 0;
}

int bz_compute_hysteresis(bz_hysteresis_t *hysteresis, double q)
{
    /* Find hysteresis loop by detecting saddle-node bifurcations in f.
     *
     * The bistable region exists between two saddle-node points f_low and f_high.
     * Uses bisection to locate them. */

    double f_lo, f_hi;
    bz_oregonator_params_t p;

    if (!hysteresis) return 0;

    /* Initialize Oregonator params with given q */
    bz_oregonator_params_init(&p, 1.0);
    p.q = q;

    /* Find first saddle-node (lower bound of hysteresis) */
    if (!bz_find_saddle_node(&f_lo, 0.1, 1.5, q, 1e-5)) {
        hysteresis->f_low = 0.0;
        hysteresis->f_high = 0.0;
        hysteresis->width = 0.0;
        hysteresis->n_steady = 1;
        return 0;
    }

    /* Find second saddle-node (upper bound of hysteresis) */
    if (!bz_find_saddle_node(&f_hi, f_lo + 0.1, 3.0, q, 1e-5)) {
        hysteresis->f_low = f_lo;
        hysteresis->f_high = f_lo;
        hysteresis->width = 0.0;
        hysteresis->n_steady = 2;
        return 0;
    }

    hysteresis->f_low  = f_lo;
    hysteresis->f_high = f_hi;
    hysteresis->width  = f_hi - f_lo;
    hysteresis->n_steady = 3;  /* bistable region has 3 steady states */
    return 1;
}

/* =========================================================================
 * L3: Steady state classification
 * ========================================================================= */

void bz_classify_steady_states(bz_stability_type_t *stability, int max_ss,
                                const bz_oregonator_state_t *ss, int n_ss,
                                const bz_oregonator_params_t *p)
{
    int i;

    if (!stability || !ss) return;

    for (i = 0; i < n_ss && i < max_ss; i++) {
        double J[9], ev_re[3], ev_im[3];
        bz_oregonator_jacobian(J, &ss[i], p);
        bz_eigenvalues_3x3(ev_re, ev_im, J);
        stability[i] = bz_classify_stability(ev_re, ev_im);
    }
}

/* =========================================================================
 * L5: Saddle-node bifurcation finding
 * ========================================================================= */

int bz_find_saddle_node(double *f_sn, double f_lo, double f_hi,
                         double q, double tol)
{
    /* Bisection to find f where saddle-node bifurcation occurs.
     * At saddle-node, the Jacobian has a zero eigenvalue => det(J) = 0.
     *
     * For the Oregonator, this means the quadratic for x has a double root:
     *   (b^2 - 4ac) = 0 where a=1, b=q+f-1, c=-q(1+f)
     *   => (q+f-1)^2 + 4q(1+f) = 0
     *   => f^2 + 2(q-1)f + (q^2-2q+1+4q+4qf) = 0  -- no, let me redo
     *
     *   (q+f-1)^2 = -4q(1+f)
     * Actually from the steady state equation:
     *   x^2 + (q+f-1)x - q(1+f) = 0
     * The discriminant is (q+f-1)^2 + 4q(1+f).
     * Saddle-node when discriminant = 0:
     *   (q+f-1)^2 + 4q(1+f) = 0
     *
     * This is a quadratic in f:
     *   f^2 + 2(q-1)f + (q^2-2q+1+4q) + 4qf = 0
     * Better: just compute numerically via bisection. */

    bz_oregonator_params_t p;
    bz_oregonator_state_t ss[4];
    int n_ss_mid;
    double f_mid;
    int iter, max_iter = 50;

    if (!f_sn) return 0;

    /* Initialize with given q */
    bz_oregonator_params_init(&p, f_lo);
    p.q = q;

    /* Check bounds have different numbers of steady states */
    {
        p.f = f_lo;
        int n_lo = bz_oregonator_steady_states(ss, &p);
        p.f = f_hi;
        int n_hi = bz_oregonator_steady_states(ss, &p);

        if (n_lo == n_hi) return 0;  /* no bifurcation in interval */
    }

    /* Bisection */
    for (iter = 0; iter < max_iter; iter++) {
        f_mid = 0.5 * (f_lo + f_hi);
        p.f = f_mid;
        n_ss_mid = bz_oregonator_steady_states(ss, &p);

        if (f_hi - f_lo < tol) {
            *f_sn = f_mid;
            return 1;
        }

        /* Check which half contains the bifurcation */
        p.f = f_lo;
        {
            int n_lo = bz_oregonator_steady_states(ss, &p);
            if (n_lo != n_ss_mid) {
                f_hi = f_mid;
            } else {
                f_lo = f_mid;
            }
        }
    }

    *f_sn = f_mid;
    return 1;
}

/* =========================================================================
 * L5: Basin of attraction
 * ========================================================================= */

int bz_basin_attractor(const bz_oregonator_state_t *s0,
                        const bz_oregonator_params_t *p,
                        const bz_oregonator_state_t *ss, int n_ss)
{
    /* Integrate from initial condition and determine which steady state
     * it converges to (within tolerance). */

    bz_oregonator_state_t s;
    double dt = 0.01;
    int i, n_steps = 5000;
    double tol = 1e-3;

    if (!s0 || !p || !ss || n_ss <= 0) return -1;

    s = *s0;

    /* Integrate forward */
    for (i = 0; i < n_steps; i++) {
        bz_oregonator_state_t dsdt;
        bz_oregonator_rhs(&dsdt, &s, p);
        s.x += dt * dsdt.x;
        s.y += dt * dsdt.y;
        s.z += dt * dsdt.z;

        /* Check convergence to any steady state */
        {
            int j;
            for (j = 0; j < n_ss; j++) {
                double dx = s.x - ss[j].x;
                double dy = s.y - ss[j].y;
                double dz = s.z - ss[j].z;
                double dist = sqrt(dx*dx + dy*dy + dz*dz);
                if (dist < tol) return j;
            }
        }
    }

    return -1;  /* no convergence */
}

void bz_basin_grid(int *basin,
                    double x_min, double x_max, int nx,
                    double z_min, double z_max, int nz,
                    double y_fixed,
                    const bz_oregonator_params_t *p,
                    const bz_oregonator_state_t *ss, int n_ss)
{
    /* Sample a 2D grid in (x,z) space at fixed y, determine
     * basin of attraction for each grid point. */

    int i, j;

    if (!basin || !p || !ss || n_ss <= 0 || nx <= 0 || nz <= 0) return;

    for (j = 0; j < nz; j++) {
        for (i = 0; i < nx; i++) {
            bz_oregonator_state_t s0;
            double x = x_min + (x_max - x_min) * i / (double)(nx - 1);
            double z = z_min + (z_max - z_min) * j / (double)(nz - 1);
            s0.x = x;
            s0.y = y_fixed;
            s0.z = z;
            basin[j * nx + i] = bz_basin_attractor(&s0, p, ss, n_ss);
        }
    }
}

/* =========================================================================
 * L7: Chemical switch (application)
 * ========================================================================= */

int bz_chemical_switch(bz_oregonator_state_t *state,
                        const bz_oregonator_params_t *p,
                        double perturb_x, int target)
{
    /* Implement a bistable chemical switch.
     * In the bistable regime, a sufficiently large perturbation
     * can flip the system from one stable steady state to the other.
     *
     * target = 0: switch to lower steady state (decrease x)
     * target = 1: switch to upper steady state (increase x) */

    bz_oregonator_state_t ss[4];
    int n_ss, attractor_before, attractor_after;
    bz_oregonator_state_t s_before;

    if (!state || !p) return 0;

    n_ss = bz_oregonator_steady_states(ss, p);
    if (n_ss < 3) return 0;  /* not bistable */

    /* Record current attractor */
    s_before = *state;
    attractor_before = bz_basin_attractor(&s_before, p, ss, n_ss);

    /* Apply perturbation */
    if (target == 0) {
        state->x -= perturb_x;
    } else {
        state->x += perturb_x;
    }

    /* Enforce physical bounds */
    if (state->x < 0.0) state->x = 0.0;
    if (state->x > 1.0) state->x = 1.0;

    /* Check new attractor */
    attractor_after = bz_basin_attractor(state, p, ss, n_ss);

    /* Switch successful if attractor changed */
    return (attractor_after != attractor_before && attractor_after >= 0) ? 1 : 0;
}
