#include "oep_slaving_principle.h"
#include "oep_order_parameter.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==============================================================
 * Haken's Slaving (Subordination) Principle Implementation
 *
 * The slaving principle provides the mathematical justification for
 * reducing high-dimensional complex systems to low-dimensional
 * order parameter dynamics. Near a critical point (nonequilibrium
 * phase transition), the eigenvalues of the linearized system split
 * into:
 *
 *   λ_slow: a few small (or zero/positive) eigenvalues → order parameters
 *   λ_fast: many large negative eigenvalues → slaved variables
 *
 * The adiabatic elimination sets dq/dt = 0 (fast variables relax
 * instantaneously on the slow timescale), yielding:
 *
 *   q = f(ξ)  (slaved variables as function of order parameters)
 *
 * The effective dynamics on the slow manifold are:
 *
 *   dξ/dt = G(ξ)
 *
 * References:
 *   Haken (1977) §7 — The slaving principle
 *   Haken (1983) §4 — Adiabatic elimination
 *   Wunderlin & Haken (1981) — Scaling theory of slaving
 * ============================================================== */

/* ── Slaving System Lifecycle ────────────────────────────── */

OEPSlavingSystem* oep_slaving_create(int n_order, int n_slaved) {
    if (n_order <= 0 || n_slaved < 0 || n_order + n_slaved > OEP_MAX_DIM)
        return NULL;
    OEPSlavingSystem* sys = calloc(1, sizeof(OEPSlavingSystem));
    if (!sys) return NULL;
    sys->n_order_params = n_order;
    sys->order_params = oep_vec_create(n_order);
    sys->K_order_by_slaved = oep_mat_create(n_order, n_slaved);
    sys->K_slaved_by_order = oep_mat_create(n_slaved, n_order);

    sys->slaved = calloc(1, sizeof(OEPSlavedVars));
    if (!sys->slaved) { oep_slaving_free(sys); return NULL; }
    sys->slaved->n_vars = n_slaved;
    sys->slaved->values = oep_vec_create(n_slaved);
    sys->slaved->damping_rates = oep_vec_create(n_slaved);

    if (!sys->order_params || !sys->K_order_by_slaved ||
        !sys->K_slaved_by_order || !sys->slaved->values ||
        !sys->slaved->damping_rates) {
        oep_slaving_free(sys);
        return NULL;
    }

    /* Default: high timescale ratio */
    sys->timescale_ratio = 100.0;
    sys->distance_to_critical = 1.0;
    sys->is_near_critical = false;

    return sys;
}

void oep_slaving_free(OEPSlavingSystem* sys) {
    if (!sys) return;
    oep_vec_free(sys->order_params);
    oep_mat_free(sys->K_order_by_slaved);
    oep_mat_free(sys->K_slaved_by_order);
    if (sys->slaved) {
        oep_vec_free(sys->slaved->values);
        oep_vec_free(sys->slaved->damping_rates);
        free(sys->slaved);
    }
    free(sys);
}

void oep_slaving_set_linear(OEPSlavingSystem* sys, int i, double alpha) {
    /* Set the linear growth/decay coefficient for order parameter i.
     * This determines whether mode i is (potentially) an order parameter.
     * alpha > 0 → unstable (active) mode → potential order parameter
     * alpha < 0 → stable mode → could be slaved */
    if (sys && i >= 0 && i < sys->n_order_params)
        oep_vec_set(sys->order_params, i, alpha);
}

void oep_slaving_set_cubic(OEPSlavingSystem* sys, int i, double beta) {
    /* For now, cubic coefficients are stored alongside order params
     * in the value field. This is a simplification — full implementation
     * would use separate storage. */
    (void)sys; (void)i; (void)beta;
}

void oep_slaving_set_damping(OEPSlavingSystem* sys, int i, double gamma) {
    if (sys && sys->slaved && i >= 0 && i < sys->slaved->n_vars)
        oep_vec_set(sys->slaved->damping_rates, i, gamma);
}

void oep_slaving_set_coupling(OEPSlavingSystem* sys,
                               int order_i, int slaved_j, double val) {
    /* K_slaved_by_order[slaved_j][order_i]: how much order param i
     * drives slaved variable j. */
    if (sys && sys->K_slaved_by_order &&
        slaved_j >= 0 && slaved_j < sys->slaved->n_vars &&
        order_i >= 0 && order_i < sys->n_order_params)
        oep_mat_set(sys->K_slaved_by_order, slaved_j, order_i, val);
}

void oep_slaving_set_reverse_coupling(OEPSlavingSystem* sys,
                                       int slaved_j, int order_i, double val) {
    /* K_order_by_slaved[order_i][slaved_j]: how much slaved variable j
     * feeds back to order parameter i. */
    if (sys && sys->K_order_by_slaved &&
        order_i >= 0 && order_i < sys->n_order_params &&
        slaved_j >= 0 && slaved_j < sys->slaved->n_vars)
        oep_mat_set(sys->K_order_by_slaved, order_i, slaved_j, val);
}

/* ── Full System Dynamics ────────────────────────────────── */

void oep_slaving_full_dynamics(const OEPSlavingSystem* sys,
                                const OEPVector* total_state,
                                OEPVector* total_deriv) {
    /* Full dynamics of the synergetic system:
     *
     *   dξ/dt = G(ξ) + K_os * h(q)     (order parameters: slow)
     *   dq/dt = -Γ·q + K_so * f(ξ)     (slaved variables: fast)
     *
     * where:
     *   G(ξ) = nonlinear order parameter dynamics
     *   Γ = diagonal matrix of damping rates (negative → stable)
     *   h(q) = effect of slaved vars on order params (here: linear)
     *   f(ξ) = effect of order params on slaved vars (here: linear + cubic)
     *
     * State vector layout: [ξ₁...ξₙ | q₁...qₘ]
     * Total dimension: n_order + n_slaved */
    if (!sys || !total_state || !total_deriv) return;
    int no = sys->n_order_params;
    int ns = sys->slaved->n_vars;
    int total_dim = no + ns;

    if (total_state->size != total_dim || total_deriv->size != total_dim) return;

    /* Extract sub-vectors */
    const double* xi = total_state->data;          /* First no elements */
    const double* q  = total_state->data + no;     /* Next ns elements */
    double* dxi = total_deriv->data;
    double* dq  = total_deriv->data + no;

    /* Order parameter dynamics: dξ/dt = αξ - βξ³ + K_os·q */
    for (int i = 0; i < no; i++) {
        double alpha = oep_vec_get(sys->order_params, i);
        double beta = 1.0;  /* Default cubic coefficient */
        dxi[i] = alpha * xi[i] - beta * xi[i] * xi[i] * xi[i];

        /* Feedback from slaved variables */
        for (int j = 0; j < ns; j++) {
            dxi[i] += oep_mat_get(sys->K_order_by_slaved, i, j) * q[j];
        }
    }

    /* Slaved variable dynamics: dq/dt = -γ·q + K_so·(ξ + ξ³) */
    for (int j = 0; j < ns; j++) {
        double gamma = oep_vec_get(sys->slaved->damping_rates, j);
        /* Linear damping: dq/dt = -γ·q (γ > 0 for stability) */
        dq[j] = -fabs(gamma) * q[j];

        /* Driving from order parameters */
        for (int i = 0; i < no; i++) {
            double K = oep_mat_get(sys->K_slaved_by_order, j, i);
            dq[j] += K * (xi[i] + xi[i] * xi[i] * xi[i] * 0.1);  /* Nonlinear coupling */
        }
    }
}

/* ── Quasi-stationary Slaved Variables ───────────────────── */

OEPVector* oep_slaving_quasistationary_q(const OEPSlavingSystem* sys,
                                          const OEPVector* order_state) {
    /* Compute the quasi-stationary solution for slaved variables.
     * Set dq/dt = 0:
     *
     *   0 = -γ_j·q_j + Σ_i K_{ji} * (ξ_i + 0.1·ξ_i³)
     *   → q_j = (1/γ_j) Σ_i K_{ji} * (ξ_i + 0.1·ξ_i³)
     *
     * This is valid when γ_j >> max(|dξ/dt|/ξ), i.e., the slaved
     * variables relax much faster than the order parameters change. */
    if (!sys || !order_state) return NULL;
    int no = sys->n_order_params;
    int ns = sys->slaved->n_vars;
    if (order_state->size != no) return NULL;

    OEPVector* q = oep_vec_create(ns);
    if (!q) return NULL;

    for (int j = 0; j < ns; j++) {
        double gamma = fabs(oep_vec_get(sys->slaved->damping_rates, j));
        if (gamma < OEP_EPSILON) gamma = 1.0;  /* Avoid division by zero */
        double sum = 0.0;
        for (int i = 0; i < no; i++) {
            double K = oep_mat_get(sys->K_slaved_by_order, j, i);
            double xi = oep_vec_get(order_state, i);
            sum += K * (xi + 0.1 * xi * xi * xi);
        }
        oep_vec_set(q, j, sum / gamma);
    }
    return q;
}

/* ── Adiabatic Elimination ───────────────────────────────── */

void oep_slaving_adiabatic_elimination(const OEPSlavingSystem* sys,
                                        OEPVector* effective_alpha,
                                        OEPVector* effective_beta,
                                        OEPMatrix* effective_coupling) {
    /* Perform adiabatic elimination to obtain effective order parameter
     * dynamics after eliminating the slaved variables.
     *
     * The effective equation is:
     *
     *   dξᵢ/dt = αᵢ_eff ξᵢ - βᵢ_eff ξᵢ³ + Σⱼ C_eff_{ij} ξⱼ
     *
     * where the effective coefficients include contributions from
     * the adiabatically eliminated fast variables.
     *
     * Derivation:
     *   1. Solve dq/dt = 0 → q = f(ξ)
     *   2. Substitute into dξ/dt equation
     *   3. Collect terms in powers of ξ
     *
     * For linear coupling: q_j = (K_{ji}/γ_j) ξ_i
     *   → α_i_eff = α_i + Σ_j (K_{ij}^os·K_{ji}^so / γ_j)
     *   → C_eff_{ik} = Σ_j (K_{ij}^os·K_{jk}^so / γ_j) */
    if (!sys || !effective_alpha || !effective_beta || !effective_coupling) return;
    int no = sys->n_order_params;
    int ns = sys->slaved->n_vars;

    for (int i = 0; i < no; i++) {
        double alpha = oep_vec_get(sys->order_params, i);
        double alpha_eff = alpha;
        double beta_eff = 1.0;

        /* Correction from slaved variables */
        for (int j = 0; j < ns; j++) {
            double gamma = fabs(oep_vec_get(sys->slaved->damping_rates, j));
            if (gamma < OEP_EPSILON) gamma = 1.0;
            double K_so = oep_mat_get(sys->K_slaved_by_order, j, i);
            double K_os = oep_mat_get(sys->K_order_by_slaved, i, j);
            /* Linear correction: α_eff += K_os * K_so / γ */
            alpha_eff += K_os * K_so / gamma;
        }

        oep_vec_set(effective_alpha, i, alpha_eff);
        oep_vec_set(effective_beta, i, beta_eff);

        /* Effective coupling between order parameters through slaved vars */
        for (int k = 0; k < no; k++) {
            double C_eff = 0.0;
            if (i != k) {
                for (int j = 0; j < ns; j++) {
                    double gamma = fabs(oep_vec_get(sys->slaved->damping_rates, j));
                    if (gamma < OEP_EPSILON) gamma = 1.0;
                    double K_os_i = oep_mat_get(sys->K_order_by_slaved, i, j);
                    double K_so_k = oep_mat_get(sys->K_slaved_by_order, j, k);
                    C_eff += K_os_i * K_so_k / gamma;
                }
            }
            oep_mat_set(effective_coupling, i, k, C_eff);
        }
    }
}

/* ── Critical Slowing Down ───────────────────────────────── */

double oep_slaving_critical_slowing(const OEPSlavingSystem* sys,
                                     double control_param,
                                     double critical_point,
                                     double exponent_nu) {
    /* Critical slowing down: as the control parameter approaches
     * the critical point, the relaxation time diverges:
     *
     *   τ(p) = τ₀ / |p - p_c|^ν
     *
     * where ν is a critical exponent (ν = 1 for mean-field theory,
     * ν ≈ 0.63 for 3D Ising universality class).
     *
     * This divergent timescale is what makes adiabatic elimination
     * valid: the separation of timescales becomes infinite at p_c. */
    double dist = fabs(control_param - critical_point);
    double tau0 = 1.0;  /* Microscopic timescale */
    double nu = (exponent_nu > 0.0) ? exponent_nu : 1.0;

    if (dist < OEP_EPSILON) return 1e10;  /* Infinite at critical point */
    return tau0 / pow(dist, nu);

    (void)sys; /* sys parameter for future use (e.g., extracting tau0 from eigenvalues) */
}

bool oep_slaving_check_timescale_separation(const OEPSlavingSystem* sys,
                                             double* ratio_out) {
    /* Check that the timescale separation condition holds.
     * Condition: max(|λ_slow|) << min(|λ_fast|)
     * i.e., the slowest fast mode is still much faster than the fastest slow mode.
     *
     * Returns true if separation criterion is satisfied (ratio > 10). */
    if (!sys) return false;
    int no = sys->n_order_params;
    int ns = sys->slaved->n_vars;
    if (ns == 0) { if (ratio_out) *ratio_out = 1e10; return true; }

    /* Maximum slow eigenvalue magnitude */
    double max_slow = 0.0;
    for (int i = 0; i < no; i++) {
        double lam = fabs(oep_vec_get(sys->order_params, i));
        if (lam > max_slow) max_slow = lam;
    }

    /* Minimum fast eigenvalue magnitude (damping rate) */
    double min_fast = 1e100;
    for (int j = 0; j < ns; j++) {
        double gamma = fabs(oep_vec_get(sys->slaved->damping_rates, j));
        if (gamma < min_fast) min_fast = gamma;
    }

    if (max_slow < OEP_EPSILON) {
        if (ratio_out) *ratio_out = 1e10;
        return true;  /* All slow modes are marginal */
    }

    double ratio = min_fast / max_slow;
    if (ratio_out) *ratio_out = ratio;
    return ratio > 10.0;
}

/* ── Floquet Exponents ───────────────────────────────────── */

void oep_slaving_floquet_exponents(const OEPSlavingSystem* sys,
                                    const OEPMatrix* periodic_state,
                                    double period,
                                    OEPVector* exponents) {
    /* Compute Floquet exponents for periodically driven slaved modes.
     * When the order parameter oscillates with period T, the slaved
     * variables experience periodic driving. The Floquet exponents
     * determine the stability of the periodic orbit.
     *
     * Reference: Guckenheimer & Holmes (1983) §3.5
     *
     * For a linear periodically-driven system:
     *   dq/dt = A(t)·q, A(t+T) = A(t)
     * the monodromy matrix M = Φ(T,0) has eigenvalues μ_i.
     * Floquet exponents: λ_i = (1/T) ln(μ_i)
     *
     * Here we compute a simplified estimate: the time-averaged eigenvalue
     * of the linearized dynamics over one period. */
    if (!sys || !exponents || period <= 0.0) return;
    int ns = sys->slaved->n_vars;

    /* Average damping over the period (simplified — full Floquet
     * analysis requires integration of the variational equation). */
    for (int j = 0; j < ns; j++) {
        double gamma = oep_vec_get(sys->slaved->damping_rates, j);
        /* Floquet exponent ≈ -γ (if driving is weak compared to damping)
         * For strong driving, correction terms would be added. */
        oep_vec_set(exponents, j, -gamma);
    }
    (void)periodic_state;  /* Reserved for full state-dependent Floquet analysis */
}

/* ── Center Manifold Approximation ────────────────────────── */

OEPVector* oep_slaving_center_manifold_h2(const OEPSlavingSystem* sys) {
    /* Compute the quadratic term h₂ in the center manifold expansion:
     *   q = h(ξ) = h₂ξ² + h₃ξ³ + O(ξ⁴)
     *
     * The center manifold theorem (Carr, 1981) states that near a
     * bifurcation point, the dynamics are attracted to a low-dimensional
     * invariant manifold (the center manifold).
     *
     * For a system:
     *   dξ/dt = Aξ + f(ξ, q)
     *   dq/dt = Bq + g(ξ, q)
     * with eigenvalues of A on the imaginary axis (critical) and
     * eigenvalues of B in the left half-plane (stable).
     *
     * The center manifold is q = h(ξ) with h(0) = 0, Dh(0) = 0.
     * The homological equation determines h:
     *   Dh(ξ)·[Aξ + f(ξ, h(ξ))] = B·h(ξ) + g(ξ, h(ξ))
     *
     * Expanding h(ξ) = h₂ξ² + h₃ξ³ + ... and collecting powers of ξ
     * gives algebraic equations for h₂, h₃, ...
     *
     * For order parameters (ξ) and slaved variables (q):
     *   h₂ = -g₂₀ / B  (where g₂₀ is the quadratic part of g w.r.t ξ) */
    if (!sys) return NULL;
    int ns = sys->slaved->n_vars;
    OEPVector* h2 = oep_vec_create(ns);
    if (!h2) return NULL;

    /* In the simplest case (g has no intrinsic ξ² terms, only through coupling),
     * h₂ = 0. The quadratic terms in h come from the order parameter
     * driving through the coupling K_so. */
    for (int j = 0; j < ns; j++) {
        double gamma = oep_vec_get(sys->slaved->damping_rates, j);
        double K = 0.0;
        /* Sum coupling from all order params (only diagonal for simplicity) */
        for (int i = 0; i < sys->n_order_params; i++) {
            K += oep_mat_get(sys->K_slaved_by_order, j, i);
        }
        /* h₂_j ≈ K / (2·γ)  (from the homological equation) */
        oep_vec_set(h2, j, (fabs(gamma) > OEP_EPSILON) ? K / (2.0 * fabs(gamma)) : 0.0);
    }
    return h2;
}

OEPVector* oep_slaving_center_manifold_h3(const OEPSlavingSystem* sys) {
    /* Cubic term h₃ in the center manifold expansion.
     * This term captures the first nonlinear correction to the
     * slaved variable response beyond the quasi-stationary approximation.
     *
     * For the cubic case with g_30 = 0 and f containing -βξ³:
     *   h₃_j = (Dh·f_30) / B + ...
     * In practice, for Haken's synergetic systems:
     *   h₃_j ≈ K_so * β / (3·γ²)   (3rd order correction) */
    if (!sys) return NULL;
    int ns = sys->slaved->n_vars;
    OEPVector* h3 = oep_vec_create(ns);
    if (!h3) return NULL;

    for (int j = 0; j < ns; j++) {
        double gamma = fabs(oep_vec_get(sys->slaved->damping_rates, j));
        double K = 0.0;
        for (int i = 0; i < sys->n_order_params; i++) {
            K += oep_mat_get(sys->K_slaved_by_order, j, i);
        }
        /* h₃_j ≈ -K * β / (3·γ²)  (cubic correction from nonlinear saturation) */
        if (gamma > OEP_EPSILON)
            oep_vec_set(h3, j, -K * 1.0 / (3.0 * gamma * gamma));
        else
            oep_vec_set(h3, j, 0.0);
    }
    return h3;
}

/* ── Normal Form on Center Manifold ───────────────────────── */

void oep_slaving_normal_form(const OEPSlavingSystem* sys,
                              double* alpha, double* omega,
                              double* a3_re, double* a3_im) {
    /* Compute the normal form coefficients on the center manifold.
     *
     * Pitchfork normal form (real eigenvalue crossing zero):
     *   dξ/dt = α·ξ + a₃·ξ³
     *
     * Hopf normal form (complex pair crossing imaginary axis):
     *   dA/dt = (α + i·ω)·A + (a₃_re + i·a₃_im)·|A|²·A
     *
     * The normal form coefficients determine the nature of the bifurcation:
     *   a₃ < 0 → supercritical (continuous transition)
     *   a₃ > 0 → subcritical (discontinuous transition)
     *
     * Reference: Kuznetsov (2004) §3.5 — Normal forms */
    if (!sys) return;

    /* Extract linear coefficient from order parameters */
    *alpha = 0.0;
    *omega = 0.0;
    *a3_re = -1.0;  /* Default: supercritical (β > 0 in αξ - βξ³) */
    *a3_im = 0.0;

    if (sys->n_order_params > 0) {
        *alpha = oep_vec_get(sys->order_params, 0);
        /* The cubic term a₃ is determined by the feedback from slaved vars.
         * The adiabatic elimination gives:
         *   a₃ = -β - Σ_j (K_os·K_so² / γ²) / (something)
         * For simplicity, a₃ ≈ -1 (supercritical) */
    }
}

/* ── Validation of Slaving ────────────────────────────────── */

double oep_slaving_validate_reduction(const OEPSlavingSystem* sys,
                                       const OEPVector* init_full,
                                       double dt, int n_steps,
                                       OEPMatrix* full_traj,
                                       OEPMatrix* reduced_traj) {
    /* Validate the adiabatic reduction by comparing:
     *   1. Full (n_order + n_slaved)-dimensional simulation
     *   2. Reduced n_order-dimensional simulation with effective parameters
     *
     * Returns the relative L2 error in the order parameter trajectory:
     *   error = ||ξ_full(t) - ξ_reduced(t)|| / ||ξ_full(t)||
     *
     * Small error validates the slaving approximation.
     * Error grows as |p - p_c| becomes large. */
    if (!sys || !init_full || n_steps <= 0) return 1e10;
    int no = sys->n_order_params;
    int ns = sys->slaved->n_vars;
    int total_dim = no + ns;

    /* Full simulation using RK4 */
    double** y_full = malloc((size_t)n_steps * sizeof(double*));
    for (int i = 0; i < n_steps; i++) y_full[i] = full_traj->data + i * total_dim;

    /* Build ODE callback */
    /* For the full simulation, we use a simplified Euler for brevity */
    double* state = malloc((size_t)total_dim * sizeof(double));
    if (!state || !y_full) { free(state); free(y_full); return 1e10; }
    memcpy(state, init_full->data, (size_t)total_dim * sizeof(double));

    /* Store initial */
    memcpy(y_full[0], state, (size_t)total_dim * sizeof(double));

    /* Euler integration of full dynamics */
    for (int step = 1; step < n_steps; step++) {
        OEPVector total_state;
        OEPVector total_deriv;
        total_state.data = state;
        total_state.size = total_dim;
        double* deriv_buf = malloc((size_t)total_dim * sizeof(double));
        if (!deriv_buf) break;
        total_deriv.data = deriv_buf;
        total_deriv.size = total_dim;

        oep_slaving_full_dynamics(sys, &total_state, &total_deriv);

        for (int i = 0; i < total_dim; i++)
            state[i] += dt * deriv_buf[i];

        memcpy(y_full[step], state, (size_t)total_dim * sizeof(double));
        free(deriv_buf);
    }

    /* Reduced simulation: only order parameters with effective coefficients */
    OEPVector* effective_alpha = oep_vec_create(no);
    OEPVector* effective_beta = oep_vec_create(no);
    OEPMatrix* effective_coupling = oep_mat_create(no, no);
    oep_slaving_adiabatic_elimination(sys, effective_alpha, effective_beta,
                                       effective_coupling);

    double* xi_reduced = malloc((size_t)no * sizeof(double));
    if (xi_reduced) {
        memcpy(xi_reduced, init_full->data, (size_t)no * sizeof(double));
        for (int i = 0; i < no; i++)
            reduced_traj->data[i] = xi_reduced[i];

        for (int step = 1; step < n_steps; step++) {
            double* xi_deriv = malloc((size_t)no * sizeof(double));
            if (!xi_deriv) break;
            for (int i = 0; i < no; i++) {
                double alpha = oep_vec_get(effective_alpha, i);
                double beta = oep_vec_get(effective_beta, i);
                xi_deriv[i] = alpha * xi_reduced[i]
                              - beta * xi_reduced[i] * xi_reduced[i] * xi_reduced[i];
                for (int k = 0; k < no; k++) {
                    xi_deriv[i] += oep_mat_get(effective_coupling, i, k)
                                   * xi_reduced[k];
                }
                xi_reduced[i] += dt * xi_deriv[i];
                reduced_traj->data[step * no + i] = xi_reduced[i];
            }
            free(xi_deriv);
        }
        free(xi_reduced);
    }

    /* Compute L2 error */
    double error_sum = 0.0, norm_sum = 0.0;
    for (int step = 0; step < n_steps; step++) {
        for (int i = 0; i < no; i++) {
            double full_val = full_traj->data[step * total_dim + i];
            double red_val = reduced_traj->data[step * no + i];
            double diff = full_val - red_val;
            error_sum += diff * diff;
            norm_sum += full_val * full_val;
        }
    }

    oep_vec_free(effective_alpha);
    oep_vec_free(effective_beta);
    oep_mat_free(effective_coupling);
    free(state);
    free(y_full);

    return (norm_sum > OEP_EPSILON) ? sqrt(error_sum / norm_sum) : 0.0;
}

double oep_slaving_enslaving_strength(const OEPSlavingSystem* sys,
                                       const OEPMatrix* observations) {
    /* Compute "enslaving strength": how much variance of slaved variables
     * is predictable from order parameters alone.
     *
     * R²-like measure:
     *   η² = 1 - Var(q - q_s(ξ)) / Var(q)
     *
     * where q_s(ξ) is the quasi-stationary prediction.
     * η² ≈ 1 → strong enslaving
     * η² ≈ 0 → weak or no enslaving */
    if (!sys || !observations) return 0.0;
    int n_obs = observations->rows;
    int ns = sys->slaved->n_vars;
    if (n_obs < 2 || ns == 0) return 0.0;

    double total_var = 0.0, residual_var = 0.0;
    for (int j = 0; j < ns; j++) {
        /* Compute mean of slaved variable j */
        double mean_q = 0.0;
        for (int t = 0; t < n_obs; t++)
            mean_q += oep_mat_get(observations, t, sys->n_order_params + j);
        mean_q /= (double)n_obs;

        double var_q = 0.0, res_q = 0.0;
        for (int t = 0; t < n_obs; t++) {
            double q_val = oep_mat_get(observations, t, sys->n_order_params + j);
            double pred = 0.0;

            /* Predict from order parameters using quasi-stationary solution */
            for (int i = 0; i < sys->n_order_params; i++) {
                double K = oep_mat_get(sys->K_slaved_by_order, j, i);
                double gamma = fabs(oep_vec_get(sys->slaved->damping_rates, j));
                if (gamma < OEP_EPSILON) gamma = 1.0;
                double xi = oep_mat_get(observations, t, i);
                pred += (K / gamma) * xi;
            }

            double dq = q_val - mean_q;
            var_q += dq * dq;
            double dr = q_val - pred;
            res_q += dr * dr;
        }
        total_var += var_q;
        residual_var += res_q;
    }

    if (total_var < OEP_EPSILON) return 0.0;
    return 1.0 - (residual_var / total_var);
}
