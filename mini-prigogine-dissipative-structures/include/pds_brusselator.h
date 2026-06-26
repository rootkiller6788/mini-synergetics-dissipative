#ifndef PDS_BRUSSELATOR_H
#define PDS_BRUSSELATOR_H

#include "pds_core.h"

/* ============================================================================
 * Prigogine-Lefever Brusselator Model
 *
 * L1 Definition: The Brusselator is a trimolecular chemical reaction model
 * proposed by Ilya Prigogine and Rene Lefever (1968) as a minimal model for
 * self-organization in chemically reacting systems far from equilibrium.
 *
 * Reaction Scheme:
 *   A -> X                    (k1)
 *   B + X -> Y + D           (k2)
 *   2X + Y -> 3X              (k3, autocatalytic)
 *   X -> E                    (k4)
 *
 * Kinetic Equations (well-mixed):
 *   dX/dt = k1*A - k2*B*X + k3*X^2*Y - k4*X
 *   dY/dt = k2*B*X - k3*X^2*Y
 *
 * Dimensionless Form (scaling):
 *   dx/dt = a - (b+1)*x + x^2*y
 *   dy/dt = b*x - x^2*y
 *
 * where:
 *   a = (k1/k4) * (k3/k4)^(1/2) * A   (normalized A concentration)
 *   b = (k2/k4) * B                     (normalized B concentration, control param)
 *
 * Key Properties:
 *   - Unique steady state: x_s = a, y_s = b/a
 *   - Hopf bifurcation at b_c = 1 + a^2
 *   - For b > b_c: stable limit cycle (chemical oscillations)
 *   - First theoretical demonstration of sustained oscillations in
 *     a simple chemical system (Nicolis & Prigogine, 1977)
 *
 * L4: The Brusselator demonstrates Prigogine's principle:
 *   "Order through fluctuations" ? beyond the critical bifurcation point,
 *   the thermodynamic branch becomes unstable and a new ordered state
 *   (the limit cycle) emerges, sustained by energy dissipation.
 *
 * References:
 *   Prigogine, I. & Lefever, R. (1968). Symmetry Breaking Instabilities
 *     in Dissipative Systems. II. J. Chem. Phys. 48, 1695-1700.
 *   Nicolis, G. & Prigogine, I. (1977). Self-Organization in
 *     Non-Equilibrium Systems. Wiley, Chapter 7.
 *   Glansdorff, P. & Prigogine, I. (1971). Thermodynamic Theory of
 *     Structure, Stability and Fluctuations. Chapter 11.
 * ============================================================================ */

/* --- L3: Brusselator Parameters --- */

typedef struct {
    double a;          /* Normalized concentration of reactant A */
    double b;          /* Normalized concentration of reactant B (bifurcation param) */
    double k1, k2, k3, k4;  /* Raw reaction rate constants */
    double A0, B0;     /* Initial concentrations (raw) */
} PDSBrusselatorParams;

/* --- L3: Brusselator State --- */

typedef struct {
    double x;          /* Concentration of species X (activator) */
    double y;          /* Concentration of species Y (inhibitor) */
} PDSBrusselatorState;

/* --- L2: Brusselator Simulation Result --- */

typedef struct {
    PDSBrusselatorState* trajectory;
    int n_points;
    double* time;
    double steady_x;
    double steady_y;
    bool is_oscillatory;
    double oscillation_period;
    double oscillation_amplitude;
    double bifurcation_parameter_b_c;  /* b_c = 1 + a^2 */
    bool has_limit_cycle;
    double max_lyapunov_exponent;
} PDSBrusselatorResult;

/* --- L2: Brusselator Spatial Extension (reaction-diffusion) --- */

typedef struct {
    PDSBrusselatorState* spatial_profile;
    int n_x;
    int n_y;
    double dx;
    double Dx;         /* Diffusion coefficient of X */
    double Dy;         /* Diffusion coefficient of Y */
    bool has_turing_pattern;
    double dominant_wavenumber;
} PDSBrusselatorSpatial;

/* ============================================================================
 * API ? Brusselator Model Functions
 * ============================================================================ */

/**
 * pds_brusselator_rhs ? Compute RHS of the ODE system.
 *
 * dx/dt = a - (b+1)*x + x^2 * y
 * dy/dt = b*x - x^2 * y
 *
 * Theorem Context: This cubic nonlinearity is the minimal form that
 * produces a Hopf bifurcation to a stable limit cycle.
 * Complexity: O(1)
 */
void pds_brusselator_rhs(const PDSBrusselatorState* state,
                          const PDSBrusselatorParams* params,
                          PDSBrusselatorState* deriv);

/**
 * pds_brusselator_steady_state ? Compute the unique steady state.
 *
 * Setting d/dt = 0:
 *   x_s = a
 *   y_s = b / a  (for a > 0)
 *
 * Theorem: The steady state is unique for all a > 0, b > 0.
 * Complexity: O(1)
 */
void pds_brusselator_steady_state(const PDSBrusselatorParams* params,
                                   double* x_s, double* y_s);

/**
 * pds_brusselator_jacobian ? Compute Jacobian at a given state.
 *
 * J = [ -(b+1) + 2*x*y,  x^2      ]
 *     [ b - 2*x*y,       -x^2     ]
 *
 * Theorem: At the steady state (x_s=a, y_s=b/a):
 *   J_ss = [ b-1,  a^2  ]
 *          [ -b,   -a^2 ]
 *   Trace(J) = b - 1 - a^2
 *   Det(J) = a^2
 *   Hopf bifurcation when Trace = 0 => b_c = 1 + a^2
 * Complexity: O(1)
 */
void pds_brusselator_jacobian(const PDSBrusselatorState* state,
                               const PDSBrusselatorParams* params,
                               double J[2][2]);

/**
 * pds_brusselator_bifurcation_point ? Compute critical b for Hopf.
 *
 * b_c = 1 + a^2
 *
 * Theorem (Nicolis & Prigogine, 1977):
 * For b < b_c: stable focus (thermodynamic branch)
 * For b = b_c: Hopf bifurcation point
 * For b > b_c: unstable focus + stable limit cycle (dissipative structure)
 * Complexity: O(1)
 */
double pds_brusselator_bifurcation_point(const PDSBrusselatorParams* params);

/**
 * pds_brusselator_rk4_step ? Single RK4 integration step.
 *
 * Uses the classic 4th-order Runge-Kutta method.
 * Complexity: O(1)
 */
void pds_brusselator_rk4_step(PDSBrusselatorState* state,
                               const PDSBrusselatorParams* params,
                               double dt);

/**
 * pds_brusselator_simulate ? Simulate the Brusselator over time interval.
 *
 * Computes trajectory from t=0 to t=duration with step dt.
 * Records state history in the result structure.
 * Complexity: O(n_steps) = O(duration/dt)
 */
PDSBrusselatorResult* pds_brusselator_simulate(
    const PDSBrusselatorState* initial,
    const PDSBrusselatorParams* params,
    double duration, double dt,
    int save_every);

/**
 * pds_brusselator_result_free ? Deallocate simulation result.
 */
void pds_brusselator_result_free(PDSBrusselatorResult* result);

/**
 * pds_brusselator_detect_limit_cycle ? Analyze trajectory for limit cycle.
 *
 * Detects oscillatory behavior by:
 * 1. Computing the period via zero-crossing of the deviation from steady state
 * 2. Computing the amplitude of oscillations
 * 3. Checking if trajectory is closed (Poincare section method)
 *
 * Theorem: A limit cycle exists if the system is beyond the Hopf bifurcation
 * (b > b_c) and the trajectory approaches a closed orbit.
 * Complexity: O(n_points)
 */
void pds_brusselator_detect_limit_cycle(PDSBrusselatorResult* result);

/**
 * pds_brusselator_spatial_alloc ? Allocate spatial Brusselator.
 */
PDSBrusselatorSpatial* pds_brusselator_spatial_create(int nx, int ny,
                                                        double Dx, double Dy);
void pds_brusselator_spatial_free(PDSBrusselatorSpatial* spat);

/**
 * pds_brusselator_spatial_rhs ? RHS for reaction-diffusion Brusselator.
 *
 * dX/dt = a - (b+1)*X + X^2*Y + Dx * nabla^2 X
 * dY/dt = b*X - X^2*Y + Dy * nabla^2 Y
 *
 * Theorem (Turing, 1952): When Dx << Dy, the homogeneous steady state
 * can become unstable to spatial perturbations, producing stationary
 * spatial patterns (Turing instability).
 * Complexity: O(nx * ny)
 */
void pds_brusselator_spatial_rhs(const PDSBrusselatorSpatial* spat,
                                  const PDSBrusselatorParams* params,
                                  PDSBrusselatorSpatial* deriv);

/**
 * pds_brusselator_euler_step_spatial ? Euler step for spatial Brusselator.
 *
 * Complexity: O(nx * ny)
 */
void pds_brusselator_euler_step_spatial(PDSBrusselatorSpatial* spat,
                                         const PDSBrusselatorParams* params,
                                         double dt);

/**
 * pds_brusselator_print_result ? Display simulation results.
 */
void pds_brusselator_print_result(const PDSBrusselatorResult* result);

#endif /* PDS_BRUSSELATOR_H */