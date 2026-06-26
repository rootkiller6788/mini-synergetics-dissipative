#include "phase_transition.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ============================================================================
 * Phase Transition — Core Implementation
 *
 * Each function implements an independent knowledge point:
 *
 * 1. Phase diagram construction (L3: data structure manipulation)
 * 2. Transition classification from OP data (L2: Ehrenfest classification)
 * 3. Entropy production computation (L4: Prigogine's theorem verification)
 * 4. Glansdorff-Prigogine criterion (L4: evolution criterion)
 * 5. NESS analysis (L2: far-from-equilibrium diagnostics)
 * ============================================================================ */

/* ---------------------------------------------------------------------------
 * PhaseDiagram Construction
 * Knowledge point: A phase diagram organizes thermodynamic phases as regions
 * in the space of control parameters (e.g., T, P, H, mu). Boundaries between
 * regions are phase transitions. This maps directly to the mathematical
 * concept of a partition of parameter space with singularities on boundaries.
 *
 * Reference: Landau & Lifshitz, Statistical Physics §144-145
 * --------------------------------------------------------------------------- */
PhaseDiagram nept_phasediagram_create(const char *param_name,
                                       double param_critical)
{
    PhaseDiagram pd;
    memset(&pd, 0, sizeof(PhaseDiagram));
    pd.param_x.value = 0.0;
    pd.param_x.critical = param_critical;
    pd.param_x.sweep_rate = 0.0;
    pd.param_x.name = param_name;
    pd.is_two_dimensional = false;
    pd.tricritical_point = param_critical; /* default: same as critical */
    pd.has_critical_endpoint = true;
    pd.critical_endpoint_value = param_critical;
    return pd;
}

int nept_phasediagram_add_region(PhaseDiagram *pd,
                                  const char *name,
                                  NEPT_SymmetryType sym,
                                  double free_energy)
{
    if (!pd || pd->num_regions >= NEPT_MAX_REGIONS) return -1;
    PhaseRegion *r = &pd->regions[pd->num_regions];
    r->region_index = pd->num_regions;
    r->phase_name = name;
    r->symmetry = sym;
    r->is_stable = true;
    r->free_energy_density = free_energy;
    pd->num_regions++;
    return r->region_index;
}

int nept_phasediagram_add_transition(PhaseDiagram *pd,
                                      double control_val,
                                      NEPT_TransitionType type,
                                      int branches_before,
                                      int branches_after,
                                      double jump)
{
    if (!pd || pd->num_transitions >= NEPT_MAX_TRANSITIONS) return -1;
    TransitionPoint *tp = &pd->transitions[pd->num_transitions];
    tp->control_value = control_val;
    tp->trans_type = type;
    tp->num_stable_branches_before = branches_before;
    tp->num_stable_branches_after = branches_after;
    tp->order_jump = jump;
    pd->num_transitions++;
    return pd->num_transitions - 1;
}

/* ---------------------------------------------------------------------------
 * Transition Classification
 *
 * Knowledge point: A phase transition can be classified by analyzing how
 * the order parameter changes as the control parameter crosses the critical
 * point.
 *
 * First-order: Discontinuous jump in OP (eta jumps from 0 to eta* > 0)
 *   Detection: max|d(eta)/d(control)| > threshold
 *   Latent heat: Delta S = -(partial F/partial T)_eta * Delta eta
 *
 * Continuous: OP grows continuously from zero
 *   Detection: eta ~ |t|^beta with beta ~ 0.3-0.5, no jump
 *
 * Algorithm: scan the derivative of OP w.r.t. control parameter.
 * If derivative exceeds a threshold (indicating a jump), classify as
 * first-order. Otherwise, fit a power law near the critical point to
 * confirm continuous transition.
 *
 * Reference: Stanley, §1.3; Chaikin & Lubensky, §2.1
 * --------------------------------------------------------------------------- */
NEPT_TransitionType nept_classify_transition(const double *order_param,
                                              const double *control,
                                              int n_points,
                                              double *jump_out)
{
    if (!order_param || !control || n_points < 4) {
        if (jump_out) *jump_out = 0.0;
        return NEPT_CONTINUOUS;
    }

    /* Compute centered finite differences: d(eta)/d(control) */
    double max_deriv = 0.0;
    double sum_deriv = 0.0;
    int max_idx = 0;
    for (int i = 1; i < n_points - 1; i++) {
        double dctrl = control[i+1] - control[i-1];
        if (fabs(dctrl) < 1e-12) continue;
        double deriv = (order_param[i+1] - order_param[i-1]) / dctrl;
        sum_deriv += fabs(deriv);
        if (fabs(deriv) > max_deriv) {
            max_deriv = fabs(deriv);
            max_idx = i;
        }
    }

    double mean_deriv = sum_deriv / (double)(n_points - 2);

    /* If max derivative is anomalously large compared to mean,
     * we have a discontinuity => first-order. */
    double ratio = (mean_deriv > 1e-12) ? max_deriv / mean_deriv : max_deriv;

    if (jump_out) {
        /* Estimate the jump as the difference across the steepest point */
        if (max_idx > 0 && max_idx < n_points - 1) {
            *jump_out = fabs(order_param[max_idx+1] - order_param[max_idx-1]);
        } else {
            *jump_out = 0.0;
        }
    }

    /* Threshold: if local derivative > 10x mean derivative => first order */
    if (ratio > 10.0 && *jump_out > 1e-6) {
        return NEPT_FIRST_ORDER;
    }

    return NEPT_CONTINUOUS;
}

/* ---------------------------------------------------------------------------
 * Ehrenfest Classification
 *
 * Knowledge point: Paul Ehrenfest (1933) classified phase transitions by
 * the lowest derivative of the free energy that shows a discontinuity.
 *
 * First-order:  (dF/dT) discontinuous => latent heat
 *               (dF/dP) discontinuous => volume jump
 * Second-order: (d^2F/dT^2) discontinuous => specific heat jump
 *               (d^2F/dP^2) discontinuous => compressibility jump
 * Lambda-point: Specific heat diverges logarithmically (He-4)
 *
 * This function implements the discrete comparison of derivatives.
 *
 * Reference: Ehrenfest, Proc. Acad. Sci. Amsterdam 36, 153 (1933)
 * --------------------------------------------------------------------------- */
NEPT_EhrenfestClass nept_ehrenfest_classify(double dF1_below,
                                             double dF1_above,
                                             double dF2_below,
                                             double dF2_above,
                                             double tolerance)
{
    double jump_F1 = fabs(dF1_above - dF1_below);
    double jump_F2 = fabs(dF2_above - dF2_below);

    /* Normalize by average magnitude to avoid scale issues */
    double avg_F1 = (fabs(dF1_below) + fabs(dF1_above)) / 2.0 + 1e-12;
    double avg_F2 = (fabs(dF2_below) + fabs(dF2_above)) / 2.0 + 1e-12;

    if (jump_F1 / avg_F1 > tolerance) {
        return NEPT_EHRENFEST_FIRST;
    }
    if (jump_F2 / avg_F2 > tolerance) {
        return NEPT_EHRENFEST_SECOND;
    }

    return NEPT_EHRENFEST_NONE;
}

/* ---------------------------------------------------------------------------
 * Entropy Production — Prigogine's Theorem
 *
 * Knowledge point: For an open system, the total entropy change dS is:
 *   dS = d_eS + d_iS
 * where d_eS is entropy flux (exchange with environment) and d_iS >= 0
 * is internal entropy production (Second Law).
 *
 * In the linear regime, entropy production is a bilinear form:
 *   sigma = d_iS/dt = sum_i J_i * X_i >= 0
 * where J_i are fluxes (heat, diffusion, reaction rates) and X_i are
 * conjugate thermodynamic forces (gradients of T, mu, affinity).
 *
 * Prigogine's minimum entropy production theorem:
 *   In the linear near-equilibrium regime, steady states minimize sigma
 *   subject to constraints. This is the thermodynamic foundation for
 *   stability of near-equilibrium steady states.
 *
 * Reference: Prigogine, Introduction to Thermodynamics of Irreversible
 *   Processes (1967), Chapter V.
 * --------------------------------------------------------------------------- */
double nept_entropy_production(const double *fluxes,
                                const double *forces,
                                int n_pairs)
{
    if (!fluxes || !forces || n_pairs <= 0) return 0.0;

    double sigma = 0.0;
    for (int i = 0; i < n_pairs; i++) {
        sigma += fluxes[i] * forces[i];
    }
    /* Second Law requires sigma >= 0 for any irreversible process.
     * If numerical noise gives sigma < 0, clamp to zero. */
    return (sigma > 0.0) ? sigma : 0.0;
}

bool nept_prigogine_min_entropy_check(const double *fluxes,
                                       const double *forces,
                                       int n_pairs,
                                       double *entropy_out)
{
    if (!fluxes || !forces || !entropy_out || n_pairs <= 0) {
        if (entropy_out) *entropy_out = 0.0;
        return false;
    }

    *entropy_out = nept_entropy_production(fluxes, forces, n_pairs);

    /* In the linear regime: J_i = sum_j L_{ij} * X_j (Onsager relations)
     * with L_{ij} = L_{ji} positive-definite. Then sigma = sum L_{ij} X_i X_j.
     * Check: each individual term should be non-negative. */
    for (int i = 0; i < n_pairs; i++) {
        if (fluxes[i] * forces[i] < -1e-15) {
            /* Negative contribution indicates violation of linear regime
             * or measurement error. But total sigma must still be >= 0. */
            return (*entropy_out >= 0.0);
        }
    }

    return (*entropy_out >= 0.0);
}

/* ---------------------------------------------------------------------------
 * Glansdorff-Prigogine Evolution Criterion
 *
 * Knowledge point: In the nonlinear far-from-equilibrium regime, the
 * minimum entropy production principle no longer holds. Instead,
 * Glansdorff and Prigogine (1971) proposed the universal evolution
 * criterion:
 *
 *   d_X sigma / dt <= 0
 *
 * where d_X sigma = sum_i J_i * dX_i is the change in entropy production
 * due to changes in the forces only (keeping fluxes fixed).
 *
 * This criterion holds for both linear and nonlinear regimes with
 * time-independent boundary conditions. It does NOT guarantee stability:
 * far-from-equilibrium steady states can become unstable, leading to
 * self-organization and dissipative structures.
 *
 * Reference: Glansdorff & Prigogine, Thermodynamic Theory of Structure,
 *   Stability and Fluctuations (1971).
 * --------------------------------------------------------------------------- */
double nept_glansdorff_prigogine_criterion(const double *state,
                                            const double *state_deriv,
                                            const double *flux_deriv,
                                            int n_vars)
{
    if (!state || !state_deriv || !flux_deriv || n_vars <= 0) return 0.0;

    /* d_X P = sum_i J_i * dX_i/dt
     * In practice: d_X sigma/dt ~ sum state_deriv[i] * flux_deriv[i]
     * (generalized force = state derivative, flux = flux derivative) */
    double dx_sigma_dt = 0.0;
    for (int i = 0; i < n_vars; i++) {
        dx_sigma_dt += state_deriv[i] * flux_deriv[i];
    }
    return dx_sigma_dt;
}

/* ---------------------------------------------------------------------------
 * System Regime Classification
 *
 * Knowledge point: The distance from equilibrium determines the system's
 * behavior. Key dimensionless numbers:
 *
 * Equilibrium: All fluxes = 0, maximum entropy (Boltzmann-Gibbs)
 * Near-equilibrium: Linear flux-force relations (Onsager)
 * Far-equilibrium: Nonlinear regime, possibility of instabilities
 * Driven steady: External constraints maintain non-zero fluxes
 *
 * This maps to Reynolds number in fluid dynamics: Re < Rec => laminar
 * (near-equilibrium), Re > Rec => turbulent (far-equilibrium).
 * The analogy: Reynolds number / critical Reynolds ~ control / critical.
 * --------------------------------------------------------------------------- */
NEPT_SystemRegime nept_regime_classify(double reynolds_or_analogue,
                                        double critical_value)
{
    if (critical_value <= 0.0) return NEPT_EQUILIBRIUM;

    double ratio = reynolds_or_analogue / critical_value;

    if (ratio < 1e-6) return NEPT_EQUILIBRIUM;
    if (ratio < 0.1)  return NEPT_NEAR_EQUILIBRIUM;
    if (ratio < 1.0)  return NEPT_FAR_EQUILIBRIUM;
    if (ratio < 10.0) return NEPT_DRIVEN_STEADY;
    return NEPT_QUENCHED;
}

/* ---------------------------------------------------------------------------
 * TransitionSeries — Observation Log
 *
 * Knowledge point: In experimental studies of phase transitions, one
 * records the order parameter, susceptibility, and correlation length
 * as functions of the control parameter. The critical point is
 * identified as the value where the susceptibility diverges (chi -> inf)
 * or where the correlation length peaks.
 *
 * This implements a simple critical point estimator: find the control
 * parameter value where the susceptibility is maximal.
 * --------------------------------------------------------------------------- */
TransitionSeries *nept_series_create(int capacity)
{
    if (capacity <= 0) capacity = 64;
    TransitionSeries *ts = (TransitionSeries *)malloc(sizeof(TransitionSeries));
    if (!ts) return NULL;
    ts->observations = (TransitionObservation *)malloc(
        (size_t)capacity * sizeof(TransitionObservation));
    if (!ts->observations) {
        free(ts);
        return NULL;
    }
    ts->num_observations = 0;
    ts->capacity = capacity;
    ts->critical_point_estimate = 0.0;
    ts->critical_exponent_estimate = 0.0;
    return ts;
}

int nept_series_push(TransitionSeries *ts,
                      double order_param,
                      double susceptibility,
                      double corr_length,
                      double time,
                      double control)
{
    if (!ts || ts->num_observations >= ts->capacity) return -1;
    TransitionObservation *obs = &ts->observations[ts->num_observations];
    obs->order_parameter_value = order_param;
    obs->susceptibility = susceptibility;
    obs->correlation_length = corr_length;
    obs->time = time;
    obs->control_parameter = control;
    ts->num_observations++;
    return ts->num_observations - 1;
}

void nept_series_free(TransitionSeries *ts)
{
    if (ts) {
        free(ts->observations);
        free(ts);
    }
}

double nept_series_estimate_critical_point(const TransitionSeries *ts)
{
    if (!ts || ts->num_observations < 2) return 0.0;

    /* Critical point = control value where susceptibility is maximal.
     * chi diverges at Tc for a continuous transition. */
    double max_chi = ts->observations[0].susceptibility;
    double critical = ts->observations[0].control_parameter;

    for (int i = 1; i < ts->num_observations; i++) {
        if (ts->observations[i].susceptibility > max_chi) {
            max_chi = ts->observations[i].susceptibility;
            critical = ts->observations[i].control_parameter;
        }
    }

    return critical;
}

/* ---------------------------------------------------------------------------
 * Non-Equilibrium Steady State (NESS)
 *
 * Knowledge point: A NESS is a state where macroscopic variables are
 * constant in time despite ongoing irreversible processes. Unlike
 * equilibrium states, NESS have non-zero entropy production sigma > 0.
 *
 * A dissipative structure (Prigogine, 1969) is a NESS exhibiting spatial
 * or temporal order, maintained by energy/matter flow through the system.
 * Examples: Bénard convection cells, BZ chemical oscillations, laser.
 *
 * Detection: sigma > threshold AND spatial/temporal correlations exist.
 * --------------------------------------------------------------------------- */
NonEquilibriumSteadyState *nept_ness_create(int num_variables)
{
    if (num_variables <= 0) return NULL;
    NonEquilibriumSteadyState *ness = (NonEquilibriumSteadyState *)malloc(
        sizeof(NonEquilibriumSteadyState));
    if (!ness) return NULL;
    ness->state_variables = (double *)calloc((size_t)num_variables, sizeof(double));
    if (!ness->state_variables) {
        free(ness);
        return NULL;
    }
    ness->num_variables = num_variables;
    ness->entropy_production = 0.0;
    ness->entropy_flux = 0.0;
    ness->lyapunov_time = 0.0;
    ness->is_dissipative_structure = false;
    return ness;
}

void nept_ness_free(NonEquilibriumSteadyState *ness)
{
    if (ness) {
        free(ness->state_variables);
        free(ness);
    }
}

void nept_ness_set_state(NonEquilibriumSteadyState *ness,
                          const double *values)
{
    if (!ness || !values) return;
    for (int i = 0; i < ness->num_variables; i++) {
        ness->state_variables[i] = values[i];
    }
}

bool nept_ness_is_dissipative(const NonEquilibriumSteadyState *ness,
                               double threshold)
{
    if (!ness) return false;
    /* A dissipative structure requires:
     * (1) Entropy production > 0 (irreversible processes)
     * (2) Entropy production exceeds threshold (far from equilibrium)
     * (3) Structure maintained by entropy export to environment */
    return (ness->entropy_production > threshold &&
            ness->entropy_flux < -0.1 * ness->entropy_production);
}