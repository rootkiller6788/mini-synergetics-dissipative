#ifndef NEPT_STOCHASTIC_DYNAMICS_H
#define NEPT_STOCHASTIC_DYNAMICS_H

#include <stdbool.h>
#include "phase_transition.h"

/* ============================================================================
 * Stochastic Dynamics of Non-Equilibrium Phase Transitions
 *
 * Near critical points, fluctuations become dominant. The deterministic
 * Landau theory must be supplemented by stochastic dynamics. Three
 * equivalent frameworks exist:
 *
 * 1. Master Equation: dP(n,t)/dt = sum_{m}[W(n|m)P(m) - W(m|n)P(n)]
 *    Describes probability evolution on a discrete state space.
 *
 * 2. Fokker-Planck Equation: dP(x,t)/dt = -d/dx[A(x)P] + (1/2) d^2/dx^2[B(x)P]
 *    Continuous approximation to the master equation (Kramers-Moyal expansion).
 *
 * 3. Langevin Equation: dx/dt = A(x) + sqrt(B(x))*xi(t)
 *    Stochastic differential equation (SDE), <xi(t)>=0, <xi(t)xi(t')>=delta(t-t').
 *
 * Key phenomena: noise-induced transitions, stochastic resonance,
 *   critical fluctuations, nucleation via rare events.
 *
 * References:
 *   N.G. van Kampen — Stochastic Processes in Physics and Chemistry (2007)
 *   C.W. Gardiner — Handbook of Stochastic Methods (2004)
 *   H. Risken — The Fokker-Planck Equation (1996)
 *   W. Horsthemke & R. Lefever — Noise-Induced Transitions (1984)
 *   L. Gammaitoni et al. — Stochastic Resonance, Rev. Mod. Phys. 70, 223 (1998)
 * ============================================================================ */

/* --- Master Equation Types --- */
typedef struct {
    int num_states;                /* Number of discrete states N              */
    double *transition_matrix;     /* W[i*N + j] = transition rate j -> i     */
    double *stationary_dist;       /* P_eq: W*P_eq = 0                        */
    double *current_state;         /* Current probability distribution        */
} MasterEquation;

/* --- Fokker-Planck Equation --- */
typedef struct {
    int n_grid;                    /* Number of grid points for spatial discretization */
    double x_min, x_max;           /* Spatial domain [x_min, x_max]           */
    double dx;                     /* Grid spacing                             */
    double *drift;                 /* A(x) evaluated at grid points            */
    double *diffusion;             /* B(x) evaluated at grid points            */
    double *probability;           /* P(x, t) current distribution             */
    double *stationary_pdf;        /* P_eq(x): steady-state solution           */
    double potential_barrier;      /* Delta Phi between stable states          */
} FokkerPlanckEq;

/* --- Langevin Equation / SDE --- */
typedef struct {
    int dimension;                 /* State dimension                          */
    double *state;                 /* Current state vector x(t)                */
    double *drift_coeff;           /* Drift coefficient A_i(x)                 */
    double *diffusion_matrix;      /* Diffusion matrix B_ij(x) [dim x dim]     */
    double *cholesky;              /* Lower-triangular Cholesky decomp of B    */
    double noise_intensity;        /* Overall noise scale factor               */
    unsigned long seed;            /* Random number generator seed              */
} LangevinEq;

/* --- Kramers Escape Rate --- */
typedef struct {
    double barrier_height;         /* Delta Phi = potential barrier            */
    double noise_intensity;        /* D = noise strength (temperature analog)  */
    double attempt_frequency;      /* omega_0: well frequency (prefactor)      */
    double barrier_frequency;      /* omega_b: curvature at barrier            */
    double escape_rate;            /* r_K = (omega_0*omega_b/2*pi)*exp(-Delta_Phi/D) */
} KramersEscape;

/* --- Noise-Induced Transition Descriptor --- */
typedef struct {
    double critical_noise;         /* Noise intensity at which transition occurs */
    int num_stable_states_before;  /* Before noise-induced transition           */
    int num_stable_states_after;  /* After noise-induced transition            */
    bool is_noise_induced;         /* True if no deterministic transition exists */
    double shift_in_maximum;       /* How much the peak of P(x) shifted         */
} NoiseInducedTransition;

/* --- Stochastic Resonance --- */
typedef struct {
    double signal_amplitude;       /* A: weak periodic driving amplitude        */
    double signal_frequency;       /* omega: driving frequency                  */
    double noise_optimum;          /* D_opt: noise level maximizing SNR         */
    double snr;                    /* Signal-to-noise ratio                     */
    double kramers_rate_match;     /* r_K ~ omega/2 for SR matching condition   */
} StochasticResonance;

/* --- Core API --- */

/* Master Equation */

/* Allocate master equation for N discrete states */
MasterEquation *nept_master_create(int num_states);

/* Set transition rate from state j to state i */
void nept_master_set_rate(MasterEquation *me, int from_state, int to_state,
                           double rate);

/* Solve for stationary distribution: W * P_eq = 0 subject to sum(P) = 1 */
/* Uses iterative method (power iteration on the transition matrix) */
int nept_master_stationary(MasterEquation *me, int max_iter, double tol);

/* Evolve probability distribution for time dt using Euler method */
void nept_master_evolve(MasterEquation *me, double dt);

/* Free the master equation */
void nept_master_free(MasterEquation *me);

/* Fokker-Planck Equation */

/* Create Fokker-Planck solver on uniform grid */
FokkerPlanckEq *nept_fpe_create(int n_grid, double x_min, double x_max);

/* Set drift coefficient A(x) */
void nept_fpe_set_drift_func(FokkerPlanckEq *fpe,
                              double (*drift_func)(double x, void *params),
                              void *params);

/* Set diffusion coefficient B(x) (assumed > 0) */
void nept_fpe_set_diffusion_func(FokkerPlanckEq *fpe,
                                  double (*diff_func)(double x, void *params),
                                  void *params);

/* Solve for stationary PDF using Crank-Nicolson implicit method */
int nept_fpe_stationary(FokkerPlanckEq *fpe, double dt, int max_iter,
                         double tolerance);

/* Evolve P(x,t) forward by dt using Crank-Nicolson */
void nept_fpe_evolve(FokkerPlanckEq *fpe, double dt);

/* Compute the effective potential Phi(x) = -integral A(x')/B(x') dx' */
void nept_fpe_effective_potential(FokkerPlanckEq *fpe);

/* Find all local maxima of the stationary PDF (stable states) */
int nept_fpe_find_stable_states(const FokkerPlanckEq *fpe,
                                 double *maxima_x, int max_maxima);

/* Free FPE solver */
void nept_fpe_free(FokkerPlanckEq *fpe);

/* Langevin Equation */

/* Allocate Langevin equation for dim-dimensional SDE */
LangevinEq *nept_langevin_create(int dim, unsigned long seed);

/* Set drift and diffusion from user functions */
void nept_langevin_set_drift(LangevinEq *le,
                              void (*drift_fun)(const double *x, double *a, void *p),
                              void *params);

void nept_langevin_set_diffusion(LangevinEq *le,
                                  void (*diff_fun)(const double *x, double *b, void *p),
                                  void *params);

/* Perform Cholesky decomposition of B for noise generation */
int nept_langevin_cholesky(LangevinEq *le);

/* Integrate SDE for one step using Euler-Maruyama scheme:
   x(t+dt) = x(t) + A(x)*dt + sqrt(B)*sqrt(dt)*N(0,1) */
void nept_langevin_step_euler_maruyama(LangevinEq *le, double dt);

/* Integrate SDE for one step using Milstein scheme (strong order 1.0) */
void nept_langevin_step_milstein(LangevinEq *le, double dt,
                                  void (*diff_deriv)(const double *x, double *db, void *p),
                                  void *params);

/* Generate a standard normal random number (Box-Muller) */
double nept_randn(unsigned long *seed);

/* Kramers Escape Rate */

/* Compute Kramers escape rate for a double-well potential */
KramersEscape nept_kramers_rate(double A_well, double A_barrier,
                                 double B_noise);

/* Estimate barrier height from FPE stationary distribution */
double nept_estimate_barrier_from_pdf(const FokkerPlanckEq *fpe);

/* Noise-Induced Transition */

/* Detect noise-induced transition: scan noise levels and count stable states */
NoiseInducedTransition nept_detect_noise_transition(
    double (*drift)(double x, double mu),
    double (*diffusion)(double x),
    double mu_control,
    double noise_min, double noise_max, int n_noise_levels);

/* Stochastic Resonance */

/* Compute SNR for a bistable system with periodic driving */
StochasticResonance nept_stochastic_resonance_analyze(
    double barrier, double driving_amp, double driving_freq,
    double noise_min, double noise_max, int n_noise);

#endif /* NEPT_STOCHASTIC_DYNAMICS_H */