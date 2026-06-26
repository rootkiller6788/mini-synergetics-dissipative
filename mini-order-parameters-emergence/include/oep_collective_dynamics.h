#ifndef OEP_COLLECTIVE_DYNAMICS_H
#define OEP_COLLECTIVE_DYNAMICS_H
#include "oep_core.h"

/* ==============================================================
 * Collective Dynamics & Self-Organization
 * Based on: Vicsek et al. (1995), Kuramoto (1984), Toner & Tu (1998)
 *
 * Collective dynamics studies how macroscopic patterns emerge from
 * microscopic interaction rules in multi-agent / multi-variable
 * complex systems.
 *
 * Key models:
 *   - Kuramoto model: coupled phase oscillators → synchronization
 *   - Vicsek model: self-propelled particles → flocking
 *   - Ising-like: spin models → ferromagnetic ordering
 *   - Reaction-diffusion: Turing pattern formation
 *   - Neural mass models: Wilson-Cowan → brain rhythms
 * ============================================================== */

/* ── Oscillator (for Kuramoto model) ───────────────────── */

typedef struct {
    double phase;             /* θᵢ: current phase in [0, 2π) */
    double natural_frequency; /* ωᵢ: intrinsic frequency */
    double coupling_strength; /* K: coupling to mean field */
} OEPOscillator;

typedef struct {
    OEPOscillator* oscillators;
    int            n_oscillators;
    double         mean_field_phase;    /* Ψ: mean field phase */
    double         mean_field_coherence;/* r: Kuramoto order parameter (0..1) */
    double         coupling_K;          /* Global coupling strength */
} OEPKuramotoSystem;

/* ── Vicsek Model Agents ───────────────────────────────── */

typedef struct {
    double x;          /* Position x */
    double y;          /* Position y */
    double vx;         /* Velocity x component */
    double vy;         /* Velocity y component */
    double speed;      /* Fixed speed magnitude v0 */
    double heading;    /* Orientation angle θ */
} OEPVicsekAgent;

typedef struct {
    OEPVicsekAgent* agents;
    int             n_agents;
    double          interaction_radius;  /* R: neighbor distance */
    double          speed_v0;            /* Fixed speed */
    double          noise_eta;           /* Angular noise amplitude */
    double          box_size;            /* Periodic boundary box L */
    double          order_parameter;     /* φ = |Σ vᵢ| / (N v0) */
    double          density;             /* ρ = N / L² */
} OEPVicsekSystem;

/* ── Ising-type Spin System ────────────────────────────── */

typedef struct {
    int*     spins;          /* sᵢ = ±1 */
    int      n_spins;
    int      grid_width;     /* Linear dimension for 2D lattice */
    int      grid_height;
    double   coupling_J;     /* Nearest-neighbor coupling */
    double   field_h;        /* External field */
    double   temperature;    /* T (in units where k_B = 1) */
    double   magnetization;  /* m = <sᵢ> */
    double   energy;         /* E = -J Σ sᵢsⱼ - h Σ sᵢ */
} OEPIsingSystem;

/* ── Wilson-Cowan Neural Mass Model ────────────────────── */

typedef struct {
    double E;           /* Excitatory population activity */
    double I;           /* Inhibitory population activity */
    double tau_E;       /* Excitatory time constant */
    double tau_I;       /* Inhibitory time constant */
    double w_EE;        /* E→E weight */
    double w_EI;        /* I→E weight */
    double w_IE;        /* E→I weight */
    double w_II;        /* I→I weight */
    double theta_E;     /* E threshold */
    double theta_I;     /* I threshold */
    double external_E;  /* External drive to E */
    double external_I;  /* External drive to I */
    double gain;        /* Sigmoid steepness */
} OEPWilsonCowan;

/* ── Kuramoto Model API ────────────────────────────────── */

OEPKuramotoSystem* oep_kura_create(int n);
void               oep_kura_free(OEPKuramotoSystem* sys);
void               oep_kura_set_frequencies_distribution(OEPKuramotoSystem* sys,
                                                          double mean, double std);
void               oep_kura_set_coupling(OEPKuramotoSystem* sys, double K);

/* Compute order parameter r*exp(i*ψ) = (1/N) Σⱼ exp(iθⱼ) */
void               oep_kura_compute_order_parameter(OEPKuramotoSystem* sys);

/* Single time step: dθᵢ/dt = ωᵢ + (K/N) Σⱼ sin(θⱼ - θᵢ) */
void               oep_kura_step(OEPKuramotoSystem* sys, double dt);

/* Evolve for n_steps with step size dt */
void               oep_kura_evolve(OEPKuramotoSystem* sys, double dt, int n_steps,
                                   double* order_param_trajectory);

/* Compute critical coupling K_c for given frequency distribution.
 * For Lorentzian g(ω) = γ/[π(ω²+γ²)], K_c = 2γ/[π g(0)] = 2γ */
double             oep_kura_critical_coupling(double freq_spread);

/* Compute the synchronization transition curve r(K). */
void               oep_kura_sync_transition(double K_min, double K_max,
                                             int n_K, double freq_spread,
                                             int n_oscillators,
                                             double* K_values, double* r_values);

/* ── Vicsek Model API ──────────────────────────────────── */

OEPVicsekSystem*   oep_vicsek_create(int n_agents, double box_size);
void               oep_vicsek_free(OEPVicsekSystem* sys);
void               oep_vicsek_random_init(OEPVicsekSystem* sys, double speed);

/* Single time step:
 *   1. Find neighbors within radius R
 *   2. Average heading of neighbors + noise
 *   3. Update positions */
void               oep_vicsek_step(OEPVicsekSystem* sys, double dt);

/* Compute order parameter φ */
double             oep_vicsek_order_parameter(const OEPVicsekSystem* sys);

/* Evolve system */
void               oep_vicsek_evolve(OEPVicsekSystem* sys, double dt, int n_steps,
                                      double* order_param_trajectory);

/* Phase diagram analysis: φ(η, ρ) */
void               oep_vicsek_phase_diagram(double noise_min, double noise_max,
                                             int n_noise, double density_min,
                                             double density_max, int n_density,
                                             int n_agents, double box_size,
                                             OEPMatrix* phi_matrix);

/* ── Ising Model API ───────────────────────────────────── */

OEPIsingSystem*    oep_ising_create(int width, int height);
void               oep_ising_free(OEPIsingSystem* sys);
void               oep_ising_random_init(OEPIsingSystem* sys);

/* Single spin flip Metropolis algorithm */
bool               oep_ising_metropolis_step(OEPIsingSystem* sys);

/* Multiple sweeps (one sweep = width*height attempted flips) */
void               oep_ising_monte_carlo_sweep(OEPIsingSystem* sys);

/* Compute current magnetization and energy */
void               oep_ising_compute_order_params(OEPIsingSystem* sys);

/* Find critical temperature via Binder cumulant:
 *   U_4 = 1 - <m⁴>/(3<m²>²) */
double             oep_ising_binder_cumulant(const double* magnetization_series,
                                              int n_samples);

/* Phase transition analysis: m(T) near T_c ≈ 2.269 for 2D Ising */
void               oep_ising_magnetization_curve(double T_min, double T_max,
                                                  int n_T, int width, int height,
                                                  int mc_sweeps,
                                                  double* T_values, double* m_values);

/* ── Wilson-Cowan Model API ────────────────────────────── */

OEPWilsonCowan*    oep_wc_create(void);
void               oep_wc_free(OEPWilsonCowan* model);
void               oep_wc_set_default_params(OEPWilsonCowan* model);

/* RHS of Wilson-Cowan equations:
 *   τ_E dE/dt = -E + S(w_EE*E - w_EI*I + external_E - theta_E)
 *   τ_I dI/dt = -I + S(w_IE*E - w_II*I + external_I - theta_I)
 * where S(x) = 1/(1 + exp(-gain*x)) */
void               oep_wc_rhs(const OEPWilsonCowan* model, double* dE, double* dI);

/* Evolve the neural mass model */
void               oep_wc_evolve(OEPWilsonCowan* model, double dt, int n_steps,
                                  double* E_traj, double* I_traj);

/* Bifurcation analysis over external drive: detect Hopf bifurcation */
bool               oep_wc_detect_oscillation(const OEPWilsonCowan* model);

/* Compute dominant frequency via FFT of activity time series */
double             oep_wc_dominant_frequency(const double* E_traj, int n_steps,
                                              double dt);

/* ── Generic Self-Organization Metrics ──────────────────── */

/* Measure the degree of self-organization: how much the global order
 * exceeds what random interactions would produce. */
double             oep_self_organization_index(const OEPMatrix* micro_states,
                                                const OEPVector* macro_order,
                                                int n_obs);

/* Spatial correlation function C(r) = <s(x)·s(x+r)> - <s>² */
void               oep_spatial_correlation(const double* field, int width, int height,
                                           int max_r, double* C_r);

/* Correlation length ξ from exponential fit C(r) ∝ exp(-r/ξ) */
double             oep_correlation_length_fit(const double* C_r, int max_r);

/* Structure factor S(k) = |Σ s(x) exp(ikx)|² / N */
void               oep_structure_factor(const double* field, int width, int height,
                                        int max_k, double* S_k);

#endif
