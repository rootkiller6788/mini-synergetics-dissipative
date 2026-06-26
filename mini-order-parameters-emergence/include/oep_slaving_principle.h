#ifndef OEP_SLAVING_PRINCIPLE_H
#define OEP_SLAVING_PRINCIPLE_H
#include "oep_core.h"

/* ==============================================================
 * Haken's Slaving (Subordination) Principle
 * Based on: Haken (1977, 1983), Wunderlin & Haken (1981)
 *
 * The slaving principle states that near a critical point, the
 * dynamics of a complex system with many degrees of freedom can
 * be reduced to the dynamics of a few "order parameters" which
 * dominate (enslave) all other "slaved" variables.
 *
 * Mathematical formulation:
 *   dξ/dt = Gξ(ξ)           — Order parameter dynamics (slow)
 *   dq/dt = Gq(ξ, q)         — Slaved variable dynamics (fast)
 *
 * Time-scale separation condition (adiabatic hypothesis):
 *   τ_slow / τ_fast → ∞  as p → pc
 *
 * The slaved variables can then be expressed as:
 *   q = q_s(ξ)  (adiabatic elimination)
 * ============================================================== */

/* ── Slaved Variable ───────────────────────────────────── */

typedef struct {
    OEPVector* values;       /* q: slaved variable values */
    OEPVector* damping_rates; /* γᵢ: damping rates (fast relaxation) */
    int        n_vars;
    double     timescale_ratio; /* τ_slow / τ_fast */
} OEPSlavedVars;

/* ── Slaving System ────────────────────────────────────── */

/* Complete synergetic system:
 *   Slow variables: ξ (order parameters)
 *   Fast variables: q (slaved, with damping rates)
 *   Coupling: bidirectional ξ ↔ q */
typedef struct {
    /* Order parameters (slow manifold) */
    OEPVector* order_params;
    int        n_order_params;

    /* Slaved variables (fast manifold) */
    OEPSlavedVars* slaved;

    /* Coupling matrices */
    OEPMatrix* K_slaved_by_order;    /* q_dot = ... + K_slaved_by_order * f(ξ) */
    OEPMatrix* K_order_by_slaved;    /* ξ_dot = ... + K_order_by_slaved * g(q) */

    /* Time scale ratio */
    double     timescale_ratio;

    /* Criticality measure */
    double     distance_to_critical;
    bool       is_near_critical;
} OEPSlavingSystem;

/* ── Slaving Principle API ─────────────────────────────── */

OEPSlavingSystem* oep_slaving_create(int n_order, int n_slaved);
void              oep_slaving_free(OEPSlavingSystem* sys);
void              oep_slaving_set_linear(OEPSlavingSystem* sys, int i, double alpha);
void              oep_slaving_set_cubic(OEPSlavingSystem* sys, int i, double beta);
void              oep_slaving_set_damping(OEPSlavingSystem* sys, int i, double gamma);
void              oep_slaving_set_coupling(OEPSlavingSystem* sys,
                                           int order_i, int slaved_j, double val);
void              oep_slaving_set_reverse_coupling(OEPSlavingSystem* sys,
                                                   int slaved_j, int order_i, double val);

/* Compute the full system dynamics (order parameters + slaved variables).
 *   total_state = [ξ₁...ξₙ | q₁...qₘ]
 *   total_deriv = [dξ₁/dt...dξₙ/dt | dq₁/dt...dqₘ/dt] */
void              oep_slaving_full_dynamics(const OEPSlavingSystem* sys,
                                            const OEPVector* total_state,
                                            OEPVector* total_deriv);

/* Compute the slaved variable quasi-stationary solution:
 *   Set dq/dt = 0, solve for q as function of ξ:
 *   q_s(ξ) = -K_slaved_by_order * f(ξ) / γ  (linear approximation) */
OEPVector*        oep_slaving_quasistationary_q(const OEPSlavingSystem* sys,
                                                const OEPVector* order_state);

/* Adiabatic elimination: project the full dynamics onto the slow manifold.
 *   Returns the effective order parameter equation after eliminating q. */
void              oep_slaving_adiabatic_elimination(const OEPSlavingSystem* sys,
                                                    OEPVector* effective_alpha,
                                                    OEPVector* effective_beta,
                                                    OEPMatrix* effective_coupling);

/* Compute the critical slowing down factor:
 *   τ(p) = τ₀ / |p - pc|^ν  (ν: critical exponent)
 * Returns relaxation time at current control parameter. */
double            oep_slaving_critical_slowing(const OEPSlavingSystem* sys,
                                               double control_param,
                                               double critical_point,
                                               double exponent_nu);

/* Time-scale separation check: verify that
 *   max(|λ_slow|) << min(|λ_fast|)
 * which is the condition for adiabatic elimination validity. */
bool              oep_slaving_check_timescale_separation(const OEPSlavingSystem* sys,
                                                         double* ratio_out);

/* Compute Floquet exponents for periodically driven slaved modes.
 * Used when the order parameter oscillates. */
void              oep_slaving_floquet_exponents(const OEPSlavingSystem* sys,
                                                const OEPMatrix* periodic_state,
                                                double period,
                                                OEPVector* exponents);

/* ── Center Manifold Approximation ─────────────────────── */

/* Center manifold theorem application: for a system near bifurcation,
 * the dynamics on the center manifold (dimension = #zero eigenvalues)
 * determine the full system behavior.
 *
 * Linearized system: dz/dt = J z
 *   with J block-diagonalized into stable/unstable/center eigenspaces. */

/* Compute the center manifold coordinates:
 *   The slaved variables on the center manifold:
 *   q = h(ξ) = h₂ξ² + h₃ξ³ + ...  (Taylor expansion)
 * Returns coefficients of h up to 3rd order. */
OEPVector*        oep_slaving_center_manifold_h2(const OEPSlavingSystem* sys);
OEPVector*        oep_slaving_center_manifold_h3(const OEPSlavingSystem* sys);

/* Normal form on center manifold:
 *   dξ/dt = αξ + a₃ξ³ + a₅ξ⁵ + ... (pitchfork normal form)
 *   dξ/dt = (α + iω)ξ + a|ξ|²ξ + ... (Hopf normal form) */
void              oep_slaving_normal_form(const OEPSlavingSystem* sys,
                                          double* alpha, double* omega,
                                          double* a3_re, double* a3_im);

/* ── Validation of Slaving ─────────────────────────────── */

/* Verify that the slaving assumption holds by comparing:
 *   (a) Full 2n-dimensional simulation
 *   (b) Reduced n-dimensional simulation (adiabatic elimination)
 * Returns relative L2 error between the two trajectories. */
double            oep_slaving_validate_reduction(const OEPSlavingSystem* sys,
                                                  const OEPVector* init_full,
                                                  double dt, int n_steps,
                                                  OEPMatrix* full_traj,
                                                  OEPMatrix* reduced_traj);

/* Compute the enslaving strength: how much variance of q is explained by ξ. */
double            oep_slaving_enslaving_strength(const OEPSlavingSystem* sys,
                                                  const OEPMatrix* observations);

#endif
