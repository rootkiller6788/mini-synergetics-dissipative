#include "oep_collective_dynamics.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==============================================================
 * Collective Dynamics & Self-Organization Implementation
 *
 * This module implements classic models of collective behavior
 * where macroscopic order parameters emerge from microscopic
 * interaction rules. Each model illustrates a different aspect
 * of the order parameter concept:
 *
 *   Kuramoto model → synchronization order parameter r
 *   Vicsek model    → flocking order parameter φ
 *   Ising model     → magnetization order parameter m
 *   Wilson-Cowan    → neural oscillation order parameter (amplitude, frequency)
 *
 * References:
 *   Kuramoto (1984) — Chemical Oscillations, Waves, and Turbulence
 *   Vicsek et al. (1995) — Phys. Rev. Lett. 75:1226
 *   Onsager (1944) — Crystal Statistics (Ising model)
 *   Wilson & Cowan (1972) — Biophys. J. 12:1
 * ============================================================== */

/* ── Kuramoto Model ────────────────────────────────────────── */

OEPKuramotoSystem* oep_kura_create(int n) {
    if (n <= 0 || n > OEP_MAX_DIM) return NULL;
    OEPKuramotoSystem* sys = calloc(1, sizeof(OEPKuramotoSystem));
    if (!sys) return NULL;
    sys->n_oscillators = n;
    sys->oscillators = calloc((size_t)n, sizeof(OEPOscillator));
    if (!sys->oscillators) { free(sys); return NULL; }
    sys->coupling_K = 1.0;
    return sys;
}

void oep_kura_free(OEPKuramotoSystem* sys) {
    if (!sys) return;
    free(sys->oscillators);
    free(sys);
}

void oep_kura_set_frequencies_distribution(OEPKuramotoSystem* sys,
                                            double mean, double std) {
    if (!sys) return;
    /* Draw natural frequencies from normal distribution.
     * The width γ of the frequency distribution determines K_c.
     * For Lorentzian: K_c = 2γ/π. For Gaussian, approximate K_c ≈ 2σ/√π. */
    for (int i = 0; i < sys->n_oscillators; i++) {
        sys->oscillators[i].natural_frequency = oep_rng_normal(mean, std);
        sys->oscillators[i].phase = oep_rng_uniform() * 2.0 * OEP_PI;
        sys->oscillators[i].coupling_strength = 1.0;
    }
}

void oep_kura_set_coupling(OEPKuramotoSystem* sys, double K) {
    if (sys) sys->coupling_K = K;
}

void oep_kura_compute_order_parameter(OEPKuramotoSystem* sys) {
    /* Compute the Kuramoto order parameter:
     *   r·exp(i·ψ) = (1/N) Σⱼ exp(i·θⱼ)
     *
     * r: coherence (0 = incoherent, 1 = fully synchronized)
     * ψ: mean phase
     *
     * This is the fundamental emergent quantity in the Kuramoto model. */
    if (!sys) return;
    int N = sys->n_oscillators;
    double sum_re = 0.0, sum_im = 0.0;
    for (int i = 0; i < N; i++) {
        sum_re += cos(sys->oscillators[i].phase);
        sum_im += sin(sys->oscillators[i].phase);
    }
    double r = sqrt(sum_re * sum_re + sum_im * sum_im) / (double)N;
    sys->mean_field_coherence = r;
    sys->mean_field_phase = atan2(sum_im, sum_re);
}

void oep_kura_step(OEPKuramotoSystem* sys, double dt) {
    /* Single time step of the Kuramoto model:
     *   dθᵢ/dt = ωᵢ + (K/N) Σⱼ sin(θⱼ - θᵢ)
     *
     * Using the mean-field: Σⱼ sin(θⱼ - θᵢ) = r·N·sin(ψ - θᵢ)
     *   dθᵢ/dt = ωᵢ + K·r·sin(ψ - θᵢ)
     *
     * We use Euler integration for simplicity. */
    if (!sys) return;
    int N = sys->n_oscillators;
    oep_kura_compute_order_parameter(sys);

    double Kr = sys->coupling_K * sys->mean_field_coherence;
    double psi = sys->mean_field_phase;

    for (int i = 0; i < N; i++) {
        double dtheta = sys->oscillators[i].natural_frequency
                        + Kr * sin(psi - sys->oscillators[i].phase);
        sys->oscillators[i].phase += dt * dtheta;
        /* Wrap phase */
        while (sys->oscillators[i].phase > 2.0 * OEP_PI)
            sys->oscillators[i].phase -= 2.0 * OEP_PI;
        while (sys->oscillators[i].phase < 0.0)
            sys->oscillators[i].phase += 2.0 * OEP_PI;
    }
}

void oep_kura_evolve(OEPKuramotoSystem* sys, double dt, int n_steps,
                      double* order_param_trajectory) {
    if (!sys || !order_param_trajectory || n_steps <= 0) return;
    for (int i = 0; i < n_steps; i++) {
        order_param_trajectory[i] = sys->mean_field_coherence;
        oep_kura_step(sys, dt);
    }
}

double oep_kura_critical_coupling(double freq_spread) {
    /* Critical coupling for the Kuramoto transition.
     * For Lorentzian g(ω) = γ/[π(ω² + γ²)]:
     *   K_c = 2 / [π·g(0)] = 2γ
     *
     * For Gaussian g(ω) = exp(-ω²/(2σ²)) / (√(2π)σ):
     *   K_c = (8/π)^{1/2}·σ ≈ 1.596·σ
     *
     * Here, freq_spread is treated as the standard deviation σ. */
    return sqrt(8.0 / OEP_PI) * freq_spread;
}

void oep_kura_sync_transition(double K_min, double K_max, int n_K,
                               double freq_spread, int n_oscillators,
                               double* K_values, double* r_values) {
    /* Compute the synchronization transition curve r(K).
     *
     * For K < K_c: r ≈ 0 (incoherent state)
     * For K > K_c: r ∝ √(K - K_c) (continuous phase transition)
     * At K = K_c: supercritical Hopf bifurcation in the mean-field equations.
     *
     * Reference: Strogatz (2000) Physica D 143:1 */
    if (n_oscillators <= 0) return;

    for (int ik = 0; ik < n_K; ik++) {
        double K = K_min + (double)ik * (K_max - K_min) / (double)(n_K - 1);

        OEPKuramotoSystem* sys = oep_kura_create(n_oscillators);
        if (!sys) continue;
        oep_kura_set_frequencies_distribution(sys, 0.0, freq_spread);
        oep_kura_set_coupling(sys, K);

        /* Evolve to steady state */
        for (int t = 0; t < 500; t++) oep_kura_step(sys, 0.01);
        /* Time average */
        double r_avg = 0.0;
        for (int t = 0; t < 100; t++) {
            oep_kura_step(sys, 0.01);
            r_avg += sys->mean_field_coherence;
        }
        r_avg /= 100.0;

        K_values[ik] = K;
        r_values[ik] = r_avg;

        oep_kura_free(sys);
    }
}

/* ── Vicsek Model ──────────────────────────────────────────── */

OEPVicsekSystem* oep_vicsek_create(int n_agents, double box_size) {
    if (n_agents <= 0 || box_size <= 0.0) return NULL;
    OEPVicsekSystem* sys = calloc(1, sizeof(OEPVicsekSystem));
    if (!sys) return NULL;
    sys->n_agents = n_agents;
    sys->agents = calloc((size_t)n_agents, sizeof(OEPVicsekAgent));
    if (!sys->agents) { free(sys); return NULL; }
    sys->box_size = box_size;
    sys->interaction_radius = box_size * 0.1;
    sys->speed_v0 = 0.5;
    sys->noise_eta = 0.1;
    sys->density = (double)n_agents / (box_size * box_size);
    return sys;
}

void oep_vicsek_free(OEPVicsekSystem* sys) {
    if (!sys) return;
    free(sys->agents);
    free(sys);
}

void oep_vicsek_random_init(OEPVicsekSystem* sys, double speed) {
    if (!sys) return;
    sys->speed_v0 = speed;
    for (int i = 0; i < sys->n_agents; i++) {
        sys->agents[i].x = oep_rng_uniform() * sys->box_size;
        sys->agents[i].y = oep_rng_uniform() * sys->box_size;
        double theta = oep_rng_uniform() * 2.0 * OEP_PI;
        sys->agents[i].heading = theta;
        sys->agents[i].vx = speed * cos(theta);
        sys->agents[i].vy = speed * sin(theta);
        sys->agents[i].speed = speed;
    }
}

void oep_vicsek_step(OEPVicsekSystem* sys, double dt) {
    if (!sys) return;
    int N = sys->n_agents;
    double R = sys->interaction_radius;
    double L = sys->box_size;

    /* Compute new headings (alignment with neighbors within radius R) */
    double* new_headings = malloc((size_t)N * sizeof(double));
    if (!new_headings) return;

    for (int i = 0; i < N; i++) {
        double sin_sum = 0.0, cos_sum = 0.0;
        int n_neighbors = 0;

        for (int j = 0; j < N; j++) {
            if (j == i) continue;

            /* Minimum image distance (periodic BC) */
            double dx = sys->agents[j].x - sys->agents[i].x;
            double dy = sys->agents[j].y - sys->agents[i].y;
            if (dx > L / 2.0) dx -= L;
            if (dx < -L / 2.0) dx += L;
            if (dy > L / 2.0) dy -= L;
            if (dy < -L / 2.0) dy += L;

            double dist = sqrt(dx * dx + dy * dy);
            if (dist < R) {
                sin_sum += sin(sys->agents[j].heading);
                cos_sum += cos(sys->agents[j].heading);
                n_neighbors++;
            }
        }

        double avg_heading;
        if (n_neighbors > 0) {
            avg_heading = atan2(sin_sum, cos_sum);
        } else {
            avg_heading = sys->agents[i].heading;  /* Keep own heading if isolated */
        }

        /* Add noise */
        avg_heading += sys->noise_eta * (oep_rng_uniform() - 0.5) * 2.0;
        new_headings[i] = avg_heading;
    }

    /* Update positions with new headings */
    double v0 = sys->speed_v0;
    for (int i = 0; i < N; i++) {
        sys->agents[i].heading = new_headings[i];
        sys->agents[i].vx = v0 * cos(new_headings[i]);
        sys->agents[i].vy = v0 * sin(new_headings[i]);

        sys->agents[i].x += dt * sys->agents[i].vx;
        sys->agents[i].y += dt * sys->agents[i].vy;

        /* Periodic boundary conditions */
        while (sys->agents[i].x < 0.0) sys->agents[i].x += L;
        while (sys->agents[i].x >= L) sys->agents[i].x -= L;
        while (sys->agents[i].y < 0.0) sys->agents[i].y += L;
        while (sys->agents[i].y >= L) sys->agents[i].y -= L;
    }

    /* Update order parameter */
    sys->order_parameter = oep_vicsek_order_parameter(sys);
    free(new_headings);
}

double oep_vicsek_order_parameter(const OEPVicsekSystem* sys) {
    /* φ = |Σ vᵢ| / (N·v₀) ∈ [0, 1]
     * φ ≈ 1: collective motion (flocking)
     * φ ≈ 0: disordered motion */
    if (!sys || sys->n_agents == 0) return 0.0;
    double sum_vx = 0.0, sum_vy = 0.0;
    for (int i = 0; i < sys->n_agents; i++) {
        sum_vx += sys->agents[i].vx;
        sum_vy += sys->agents[i].vy;
    }
    double norm = sqrt(sum_vx * sum_vx + sum_vy * sum_vy);
    return norm / ((double)sys->n_agents * sys->speed_v0);
}

void oep_vicsek_evolve(OEPVicsekSystem* sys, double dt, int n_steps,
                        double* order_param_trajectory) {
    if (!sys || !order_param_trajectory || n_steps <= 0) return;
    for (int i = 0; i < n_steps; i++) {
        order_param_trajectory[i] = sys->order_parameter;
        oep_vicsek_step(sys, dt);
    }
}

void oep_vicsek_phase_diagram(double noise_min, double noise_max,
                               int n_noise, double density_min,
                               double density_max, int n_density,
                               int n_agents, double box_size,
                               OEPMatrix* phi_matrix) {
    if (!phi_matrix) return;
    for (int ni = 0; ni < n_noise; ni++) {
        double eta = noise_min + (double)ni * (noise_max - noise_min) / (double)(n_noise - 1);
        for (int nd = 0; nd < n_density; nd++) {
            double rho = density_min + (double)nd * (density_max - density_min) / (double)(n_density - 1);

            OEPVicsekSystem* sys = oep_vicsek_create(n_agents, box_size);
            if (!sys) continue;
            oep_vicsek_random_init(sys, 0.5);
            sys->noise_eta = eta;

            /* Evolve and time-average φ */
            for (int t = 0; t < 200; t++) oep_vicsek_step(sys, 1.0);
            double phi_avg = 0.0;
            for (int t = 0; t < 50; t++) {
                oep_vicsek_step(sys, 1.0);
                phi_avg += sys->order_parameter;
            }
            oep_mat_set(phi_matrix, ni, nd, phi_avg / 50.0);
            oep_vicsek_free(sys);
        }
    }
}

/* ── Ising Model ───────────────────────────────────────────── */

OEPIsingSystem* oep_ising_create(int width, int height) {
    if (width <= 0 || height <= 0) return NULL;
    OEPIsingSystem* sys = calloc(1, sizeof(OEPIsingSystem));
    if (!sys) return NULL;
    sys->grid_width = width;
    sys->grid_height = height;
    sys->n_spins = width * height;
    sys->spins = calloc((size_t)sys->n_spins, sizeof(int));
    if (!sys->spins) { free(sys); return NULL; }
    sys->coupling_J = 1.0;
    sys->temperature = 2.0;  /* Below T_c ≈ 2.269 for 2D Ising */
    return sys;
}

void oep_ising_free(OEPIsingSystem* sys) {
    if (!sys) return;
    free(sys->spins);
    free(sys);
}

void oep_ising_random_init(OEPIsingSystem* sys) {
    if (!sys) return;
    for (int i = 0; i < sys->n_spins; i++)
        sys->spins[i] = (oep_rng_uniform() > 0.5) ? 1 : -1;
}

bool oep_ising_metropolis_step(OEPIsingSystem* sys) {
    /* Single spin-flip Metropolis algorithm.
     *
     * 1. Choose random spin s_i
     * 2. Compute energy change ΔE if spin is flipped
     * 3. Accept flip with probability min(1, exp(-ΔE/kT))
     *
     * For 2D square lattice with nearest-neighbor coupling J:
     *   ΔE = 2·J·s_i·Σ_{j∈nn(i)} s_j
     *
     * Returns true if flip was attempted (and accepted or rejected). */
    if (!sys) return false;

    int w = sys->grid_width, h = sys->grid_height;
    int i = (int)(oep_rng_uniform() * (double)sys->n_spins);
    int ix = i % w, iy = i / w;

    /* Nearest neighbors (periodic BC) */
    int left  = ((ix - 1 + w) % w) + iy * w;
    int right = ((ix + 1) % w) + iy * w;
    int up    = ix + ((iy - 1 + h) % h) * w;
    int down  = ix + ((iy + 1) % h) * w;

    double sum_nn = (double)(sys->spins[left] + sys->spins[right]
                            + sys->spins[up] + sys->spins[down]);
    double dE = 2.0 * sys->coupling_J * (double)sys->spins[i] * sum_nn
                + 2.0 * sys->field_h * (double)sys->spins[i];

    /* Metropolis acceptance */
    if (dE <= 0.0 || oep_rng_uniform() < exp(-dE / sys->temperature)) {
        sys->spins[i] = -sys->spins[i];
        return true;
    }
    return false;
}

void oep_ising_monte_carlo_sweep(OEPIsingSystem* sys) {
    /* One Monte Carlo sweep = N attempted spin flips (N = total spins).
     * After many sweeps, the system equilibrates at temperature T. */
    if (!sys) return;
    for (int i = 0; i < sys->n_spins; i++)
        oep_ising_metropolis_step(sys);
}

void oep_ising_compute_order_params(OEPIsingSystem* sys) {
    if (!sys) return;
    int N = sys->n_spins;
    double sum_s = 0.0;
    double E = 0.0;
    int w = sys->grid_width, h = sys->grid_height;

    for (int iy = 0; iy < h; iy++) {
        for (int ix = 0; ix < w; ix++) {
            int i = ix + iy * w;
            int s = sys->spins[i];
            sum_s += (double)s;

            /* Interaction energy (count each bond once) */
            int right = ((ix + 1) % w) + iy * w;
            int down = ix + ((iy + 1) % h) * w;
            E -= sys->coupling_J * (double)s
                 * (double)(sys->spins[right] + sys->spins[down]);
        }
    }
    sys->magnetization = sum_s / (double)N;
    sys->energy = E - sys->field_h * sum_s;  /* E = -J Σ s_i s_j - h Σ s_i */
}

double oep_ising_binder_cumulant(const double* magnetization_series,
                                  int n_samples) {
    /* Binder cumulant:
     *   U_4 = 1 - <m⁴> / (3 <m²>²)
     *
     * At the critical point T_c, U_4 is independent of system size
     * (crossing point for different L). This provides a precise
     * method for locating T_c without knowing critical exponents.
     *
     * For the 2D Ising model: U_4(T_c) ≈ 0.611 (Onsager) */
    if (!magnetization_series || n_samples < 4) return 0.0;

    double m2_avg = 0.0, m4_avg = 0.0;
    for (int i = 0; i < n_samples; i++) {
        double m = magnetization_series[i];
        double m2 = m * m;
        m2_avg += m2;
        m4_avg += m2 * m2;
    }
    m2_avg /= (double)n_samples;
    m4_avg /= (double)n_samples;

    if (m2_avg < OEP_EPSILON) return 1e10;  /* Avoid division by zero */
    return 1.0 - m4_avg / (3.0 * m2_avg * m2_avg);
}

void oep_ising_magnetization_curve(double T_min, double T_max,
                                    int n_T, int width, int height,
                                    int mc_sweeps,
                                    double* T_values, double* m_values) {
    /* Compute magnetization curve m(T) for the 2D Ising model.
     *
     * For T < T_c: spontaneous magnetization m ≠ 0 (ferromagnetic)
     * For T > T_c: m = 0 (paramagnetic)
     * Near T_c: m ∝ (T_c - T)^β with β = 1/8 (exact 2D Ising) */
    for (int it = 0; it < n_T; it++) {
        double T = T_min + (double)it * (T_max - T_min) / (double)(n_T - 1);
        OEPIsingSystem* sys = oep_ising_create(width, height);
        if (!sys) continue;
        sys->temperature = T;
        oep_ising_random_init(sys);

        /* Equilibrate */
        int eq_sweeps = mc_sweeps / 2;
        for (int s = 0; s < eq_sweeps; s++) oep_ising_monte_carlo_sweep(sys);

        /* Measure */
        double m_sum = 0.0;
        for (int s = 0; s < eq_sweeps; s++) {
            oep_ising_monte_carlo_sweep(sys);
            oep_ising_compute_order_params(sys);
            m_sum += sys->magnetization;
        }
        T_values[it] = T;
        m_values[it] = fabs(m_sum / (double)eq_sweeps);
        oep_ising_free(sys);
    }
}

/* ── Wilson-Cowan Model ────────────────────────────────────── */

OEPWilsonCowan* oep_wc_create(void) {
    OEPWilsonCowan* model = calloc(1, sizeof(OEPWilsonCowan));
    return model;
}

void oep_wc_free(OEPWilsonCowan* model) { free(model); }

void oep_wc_set_default_params(OEPWilsonCowan* model) {
    if (!model) return;
    model->tau_E = 10.0;
    model->tau_I = 20.0;
    model->w_EE = 12.0;
    model->w_EI = 4.0;
    model->w_IE = 13.0;
    model->w_II = 11.0;
    model->theta_E = 2.5;
    model->theta_I = 3.5;
    model->external_E = 1.0;
    model->external_I = 0.0;
    model->gain = 1.0;
    model->E = 0.1;
    model->I = 0.1;
}

void oep_wc_rhs(const OEPWilsonCowan* model, double* dE, double* dI) {
    /* Wilson-Cowan equations:
     *
     *   τ_E·dE/dt = -E + S(w_EE·E - w_EI·I + external_E - θ_E)
     *   τ_I·dI/dt = -I + S(w_IE·E - w_II·I + external_I - θ_I)
     *
     * where S(x) = 1 / (1 + exp(-gain·x)) is the sigmoid response function.
     *
     * This model exhibits:
     *   - Fixed points (steady activity)
     *   - Limit cycles (oscillations → brain rhythms)
     *   - Bistability (two stable fixed points → working memory)
     *   - Hopf bifurcation (transition to oscillations as drive changes) */
    if (!model || !dE || !dI) return;

    double sE = 1.0 / (1.0 + exp(-model->gain * (
        model->w_EE * model->E - model->w_EI * model->I
        + model->external_E - model->theta_E)));

    double sI = 1.0 / (1.0 + exp(-model->gain * (
        model->w_IE * model->E - model->w_II * model->I
        + model->external_I - model->theta_I)));

    *dE = (-model->E + sE) / model->tau_E;
    *dI = (-model->I + sI) / model->tau_I;
}

void oep_wc_evolve(OEPWilsonCowan* model, double dt, int n_steps,
                    double* E_traj, double* I_traj) {
    if (!model || !E_traj || !I_traj || n_steps <= 0) return;

    E_traj[0] = model->E;
    I_traj[0] = model->I;

    for (int i = 1; i < n_steps; i++) {
        double dE, dI;
        oep_wc_rhs(model, &dE, &dI);
        model->E += dt * dE;
        model->I += dt * dI;
        /* Clamp to [0, 1] */
        if (model->E < 0.0) model->E = 0.0;
        if (model->E > 1.0) model->E = 1.0;
        if (model->I < 0.0) model->I = 0.0;
        if (model->I > 1.0) model->I = 1.0;
        E_traj[i] = model->E;
        I_traj[i] = model->I;
    }
}

bool oep_wc_detect_oscillation(const OEPWilsonCowan* model) {
    /* Detect if the Wilson-Cowan model is in an oscillatory regime.
     *
     * Criterion: compute eigenvalues of the Jacobian at the fixed point.
     * If eigenvalues are complex with positive real part → oscillations.
     *
     * Jacobian: J = [-1/τ_E + ∂S_E/∂E·w_EE/τ_E,  ∂S_E/∂I·(-w_EI)/τ_E;
     *                ∂S_I/∂E·w_IE/τ_I,           -1/τ_I + ∂S_I/∂I·(-w_II)/τ_I]
     *
     * Simplified check: simulate a few steps and check for sign changes
     * in the derivative (indicates oscillation). */
    if (!model) return false;

    /* Quick linear stability check at current state */
    double sE = 1.0 / (1.0 + exp(-model->gain * (
        model->w_EE * model->E - model->w_EI * model->I
        + model->external_E - model->theta_E)));
    double sI = 1.0 / (1.0 + exp(-model->gain * (
        model->w_IE * model->E - model->w_II * model->I
        + model->external_I - model->theta_I)));

    double dsE = model->gain * sE * (1.0 - sE);
    double dsI = model->gain * sI * (1.0 - sI);

    /* Jacobian entries */
    double J_EE = (-1.0 + dsE * model->w_EE) / model->tau_E;
    double J_EI = (-dsE * model->w_EI) / model->tau_E;
    double J_IE = (dsI * model->w_IE) / model->tau_I;
    double J_II = (-1.0 - dsI * model->w_II) / model->tau_I;

    /* Trace and determinant */
    double trace = J_EE + J_II;
    double det = J_EE * J_II - J_EI * J_IE;
    double disc = trace * trace - 4.0 * det;

    /* Complex eigenvalues with positive real part */
    return (disc < 0.0) && (trace > 0.0);
}

double oep_wc_dominant_frequency(const double* E_traj, int n_steps,
                                  double dt) {
    /* Estimate dominant frequency from E-trajectory using zero-crossings.
     *
     * For a limit cycle, count upward zero crossings of the centered signal.
     * Frequency f = (n_crossings_unique) / (2·T_total) */
    if (!E_traj || n_steps < 10 || dt <= 0.0) return 0.0;

    /* Center the signal */
    double mu = oep_ts_mean(E_traj, n_steps);

    int crossings = 0;
    for (int i = 1; i < n_steps; i++) {
        if ((E_traj[i-1] - mu) * (E_traj[i] - mu) < 0.0)
            crossings++;
    }

    /* Each cycle has 2 zero crossings */
    if (crossings < 2) return 0.0;
    double period = 2.0 * (double)(n_steps - 1) * dt / (double)crossings;
    return 1.0 / period;
}

/* ── Self-Organization Metrics ─────────────────────────────── */

double oep_self_organization_index(const OEPMatrix* micro_states,
                                    const OEPVector* macro_order,
                                    int n_obs) {
    /* Self-organization index (SOI):
     *   SOI = 1 - H(micro | macro) / H(micro)
     *
     * where H(micro|macro) is the conditional entropy of micro states
     * given the macro order parameter.
     *
     * SOI ≈ 1: strong self-organization (macro reduces micro uncertainty)
     * SOI ≈ 0: weak/no self-organization */
    if (!micro_states || !macro_order || n_obs < 2) return 0.0;

    /* Simplified: use variance reduction as proxy for entropy reduction.
     * Full entropy computation requires binning and is computationally
     * more expensive. */
    int n_micro = micro_states->cols;
    double var_unexplained = 0.0, var_total = 0.0;

    for (int j = 0; j < n_micro; j++) {
        /* Total variance of micro variable j */
        double* col_data = malloc((size_t)n_obs * sizeof(double));
        if (!col_data) continue;
        for (int t = 0; t < n_obs; t++) col_data[t] = oep_mat_get(micro_states, t, j);
        double var_total_j = oep_ts_variance(col_data, n_obs);

        /* Variance explained by macro: use linear regression residual */
        double mu_micro = oep_ts_mean(col_data, n_obs);
        double mu_macro = oep_ts_mean(macro_order->data, n_obs);
        double cov = 0.0, var_macro = 0.0;
        for (int t = 0; t < n_obs; t++) {
            cov += (col_data[t] - mu_micro) * (macro_order->data[t] - mu_macro);
            var_macro += (macro_order->data[t] - mu_macro)
                         * (macro_order->data[t] - mu_macro);
        }
        double slope = (var_macro > OEP_EPSILON) ? (cov / var_macro) : 0.0;
        double intercept = mu_micro - slope * mu_macro;

        double residual_var = 0.0;
        for (int t = 0; t < n_obs; t++) {
            double pred = slope * macro_order->data[t] + intercept;
            double resid = col_data[t] - pred;
            residual_var += resid * resid;
        }
        residual_var /= (double)(n_obs - 2);

        var_total += var_total_j;
        var_unexplained += residual_var;
        free(col_data);
    }

    if (var_total < OEP_EPSILON) return 0.0;
    return 1.0 - (var_unexplained / var_total);
}

void oep_spatial_correlation(const double* field, int width, int height,
                              int max_r, double* C_r) {
    /* Spatial correlation function:
     *   C(r) = <ξ(x)ξ(x+r)> - <ξ>²
     *
     * For an Ising-like field near T_c: C(r) ∝ r^{-η} exp(-r/ξ)
     * where ξ is the correlation length. */
    if (!field || !C_r || max_r <= 0) return;

    double mu = oep_ts_mean(field, width * height);

    for (int r = 0; r < max_r && r < width / 2; r++) {
        double sum = 0.0;
        int count = 0;
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                int idx = y * width + x;
                int idx_r = y * width + ((x + r) % width);
                sum += (field[idx] - mu) * (field[idx_r] - mu);
                count++;
            }
        }
        C_r[r] = (count > 0) ? (sum / (double)count) : 0.0;
    }
}

double oep_correlation_length_fit(const double* C_r, int max_r) {
    /* Estimate correlation length ξ from exponential fit:
     *   C(r) ≈ A·exp(-r/ξ)
     * → log C(r) ≈ log A - r/ξ
     *
     * Linear regression of log(C(r)) vs r gives ξ = -1/slope. */
    if (!C_r || max_r < 3) return 0.0;

    /* Find first zero crossing or significantly negative value */
    int valid_r = max_r;
    for (int r = 0; r < max_r; r++) {
        if (C_r[r] < OEP_EPSILON) { valid_r = r; break; }
    }
    if (valid_r < 2) return 0.0;

    double s_r = 0.0, s_lnC = 0.0, s_r2 = 0.0, s_r_lnC = 0.0;
    int n = 0;
    for (int r = 0; r < valid_r; r++) {
        double lnC = log(C_r[r]);
        s_r += (double)r;
        s_lnC += lnC;
        s_r2 += (double)(r * r);
        s_r_lnC += (double)r * lnC;
        n++;
    }
    if (n < 2) return 0.0;

    double slope = ((double)n * s_r_lnC - s_r * s_lnC)
                   / ((double)n * s_r2 - s_r * s_r);
    if (slope >= 0.0) return 1e10;  /* No decay → infinite correlation */
    return -1.0 / slope;
}

void oep_structure_factor(const double* field, int width, int height,
                           int max_k, double* S_k) {
    /* Structure factor (static):
     *   S(k) = |Σ_x ξ(x) exp(ik·x)|² / N
     *
     * For a periodic pattern with wavelength λ: S(k) peaks at k = 2π/λ.
     *
     * Near a phase transition: S(k) ∝ k^{-2+η} for large k
     * (Ornstein-Zernike form: S(k) = k_B T / (a + c k²)) */
    if (!field || !S_k || max_k <= 0) return;

    int N = width * height;
    double mu = oep_ts_mean(field, N);

    for (int k = 0; k < max_k && k < width / 2; k++) {
        double kx = 2.0 * OEP_PI * (double)k / (double)width;
        double sum_re = 0.0, sum_im = 0.0;

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                double val = field[y * width + x] - mu;
                double theta = kx * (double)x;
                sum_re += val * cos(theta);
                sum_im += val * sin(theta);
            }
        }
        S_k[k] = (sum_re * sum_re + sum_im * sum_im) / (double)N;
    }
}
