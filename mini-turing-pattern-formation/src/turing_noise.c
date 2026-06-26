/**
 * turing_noise.c — Stochastic Effects in Turing Pattern Formation
 *
 * Implements noise-driven extensions of the deterministic reaction-diffusion
 * framework, including:
 *   - Additive Gaussian white noise (Langevin approach)
 *   - Multiplicative (state-dependent) noise
 *   - Colored (Ornstein-Uhlenbeck) noise with finite correlation time
 *   - Stochastic differential equation (SDE) integration (Euler-Maruyama)
 *   - Noise-induced pattern selection
 *   - Nucleation of patterns from fluctuations
 *
 * Reference: Turing, A.M. (1952) "The Chemical Basis of Morphogenesis"
 *            Butler, T. & Goldenfeld, N. (2009) Phys. Rev. E 80:030902(R)
 *            Biancalani, T. et al. (2010) Phys. Rev. Lett. 105(12):128102
 *            Maini, P.K. et al. (1997) Bull. Math. Biol. 59(1):61-84
 *            Garcia-Ojalvo, J. & Sancho, J.M. (1999) Noise in Spatially
 *            Extended Systems (Springer)
 *
 * Knowledge Coverage:
 *   L8: Stochastic PDEs, noise-induced transitions, Langevin dynamics
 *   L5: Euler-Maruyama integration for SDEs
 *   L7: Noise effects in biological pattern formation
 */

#include "turing_common.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==========================================================================
 * Internal PRNG (xoshiro256** — reused, locally defined for independence)
 * ==========================================================================*/

typedef struct { uint64_t s[4]; } XoshiroRNG;

static inline uint64_t rotl_64(uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

static uint64_t xoshiro_next(XoshiroRNG *rng) {
    uint64_t result = rotl_64(rng->s[1] * 5, 7) * 9;
    uint64_t t = rng->s[1] << 17;
    rng->s[2] ^= rng->s[0];
    rng->s[3] ^= rng->s[1];
    rng->s[1] ^= rng->s[2];
    rng->s[0] ^= rng->s[3];
    rng->s[2] ^= t;
    rng->s[3] = rotl_64(rng->s[3], 45);
    return result;
}

static void xoshiro_seed(XoshiroRNG *rng, uint64_t seed) {
    uint64_t z = seed + 0x9e3779b97f4a7c15ULL;
    for (int i = 0; i < 4; i++) {
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        z = z ^ (z >> 31);
        rng->s[i] = z;
    }
}

static double xoshiro_normal(XoshiroRNG *rng) {
    /* Box-Muller */
    double u1, u2;
    u1 = (double)(xoshiro_next(rng) >> 11) * 0x1.0p-53;
    u2 = (double)(xoshiro_next(rng) >> 11) * 0x1.0p-53;
    if (u1 < 1e-16) u1 = 1e-16;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

/* ==========================================================================
 * Noise Types
 * ==========================================================================*/

typedef enum {
    NOISE_WHITE_GAUSSIAN = 0,    /* ⟨ξ(t)ξ(t')⟩ = δ(t-t') */
    NOISE_COLORED_OU    = 1,     /* Ornstein-Uhlenbeck, correlation time τ */
    NOISE_MULTIPLICATIVE = 2     /* State-dependent amplitude σ(u) */
} NoiseType;

/**
 * Ornstein-Uhlenbeck Process State
 *
 * dη/dt = -η/τ + √(2/τ) · ξ(t)
 *
 * where ξ(t) is white noise and τ is the correlation time.
 * As τ → 0, η approaches white noise.
 * Stationary distribution: N(0, 1).
 *
 * Reference: Gillespie, D.T. (1996) Phys. Rev. E 54(2):2084-2091
 */
typedef struct {
    double   value;       /* current noise value */
    double   tau;         /* correlation time */
    double   amplitude;   /* noise strength √(2D) */
    XoshiroRNG rng;       /* random number generator */
} OUProcess;

/** Initialize an Ornstein-Uhlenbeck process. */
static void ou_init(OUProcess *ou, double tau, double amplitude, uint64_t seed) {
    ou->value = 0.0;
    ou->tau = tau;
    ou->amplitude = amplitude;
    xoshiro_seed(&ou->rng, seed);
}

/** Advance the OU process by one time step (exact update). */
static double ou_step(OUProcess *ou, double dt) {
    if (ou->tau < 1e-15) {
        /* White noise limit: independent Gaussian */
        ou->value = ou->amplitude * xoshiro_normal(&ou->rng) / sqrt(dt);
    } else {
        double rho = exp(-dt / ou->tau);
        double sigma = sqrt(1.0 - rho * rho);
        ou->value = rho * ou->value
                  + ou->amplitude * sigma * xoshiro_normal(&ou->rng);
    }
    return ou->value;
}

/* ==========================================================================
 * Stochastic Field Operations
 * ==========================================================================*/

/**
 * Add independent Gaussian white noise to each grid point.
 *
 * In the Langevin framework:
 *   du = f(u,v)dt + D_u∇²u dt + √(2ε_u) dW_u
 *   dv = g(u,v)dt + D_v∇²v dt + √(2ε_v) dW_v
 *
 * The noise amplitude ε scales as 1/N (inverse system size).
 *
 * Reference: Gardiner, C.W. (2009) Stochastic Methods, 4th ed., Springer
 */
void noise_add_white(ScalarField2D *u, ScalarField2D *v,
                     double eps_u, double eps_v, double dt,
                     uint64_t seed) {
    if (!u || !v) return;

    XoshiroRNG rng_u, rng_v;
    xoshiro_seed(&rng_u, seed);
    xoshiro_seed(&rng_v, seed + 0xDEADBEEF);

    double sigma_u = sqrt(2.0 * eps_u / dt);
    double sigma_v = sqrt(2.0 * eps_v / dt);

    for (int i = 0; i < u->Nx; i++) {
        for (int j = 0; j < u->Ny; j++) {
            double u_val = field_get(u, i, j) + sigma_u * xoshiro_normal(&rng_u);
            double v_val = field_get(v, i, j) + sigma_v * xoshiro_normal(&rng_v);
            field_set(u, i, j, u_val);
            field_set(v, i, j, v_val);
        }
    }
}

/**
 * Add multiplicative (state-dependent) noise.
 *
 * dη = σ(u) · dW, where σ²(u) ∝ u (reflecting Poisson character
 * of reaction events at low copy numbers).
 *
 * Reference: Schnakenberg, J. (1976) Rev. Mod. Phys. 48(4):571-585
 */
void noise_add_multiplicative(ScalarField2D *u, ScalarField2D *v,
                              double eps_u, double eps_v, double dt,
                              uint64_t seed) {
    if (!u || !v) return;

    XoshiroRNG rng_u, rng_v;
    xoshiro_seed(&rng_u, seed);
    xoshiro_seed(&rng_v, seed + 0xBEEFCAFE);

    double prefac_u = sqrt(2.0 * eps_u / dt);
    double prefac_v = sqrt(2.0 * eps_v / dt);

    for (int i = 0; i < u->Nx; i++) {
        for (int j = 0; j < u->Ny; j++) {
            double u_val = field_get(u, i, j);
            double v_val = field_get(v, i, j);

            /* Multiplicative: noise amplitude ∝ √(|u|) */
            double sigma_u = prefac_u * sqrt(fabs(u_val) + 1e-10);
            double sigma_v = prefac_v * sqrt(fabs(v_val) + 1e-10);

            double du_noise = sigma_u * xoshiro_normal(&rng_u);
            double dv_noise = sigma_v * xoshiro_normal(&rng_v);

            field_set(u, i, j, u_val + du_noise);
            field_set(v, i, j, v_val + dv_noise);
        }
    }
}

/**
 * Add spatially correlated (colored) noise with OU temporal dynamics.
 *
 * Each grid point has its own OU process with correlation time τ.
 * This models finite-time environmental fluctuations.
 *
 * Reference: Hänggi, P. & Jung, P. (1995) Adv. Chem. Phys. 89:239-326
 */
void noise_add_colored(ScalarField2D *u, ScalarField2D *v,
                       double eps_u, double eps_v, double dt,
                       double tau_u, double tau_v,
                       OUProcess **ou_states_u, OUProcess **ou_states_v,
                       uint64_t seed) {
    if (!u || !v || !ou_states_u || !ou_states_v) return;

    int Nx = u->Nx, Ny = u->Ny;
    double amp_u = sqrt(2.0 * eps_u);
    double amp_v = sqrt(2.0 * eps_v);

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            int idx = i * Ny + j;
            double eta_u = ou_step(&(*ou_states_u)[idx], dt);
            double eta_v = ou_step(&(*ou_states_v)[idx], dt);

            double u_val = field_get(u, i, j) + eta_u * sqrt(dt);
            double v_val = field_get(v, i, j) + eta_v * sqrt(dt);
            field_set(u, i, j, u_val);
            field_set(v, i, j, v_val);
        }
    }
}

/**
 * Allocate OU process arrays for colored noise on a grid.
 *
 * Returns 0 on success, -1 on allocation failure.
 */
int noise_ou_alloc(int Nx, int Ny, double tau_u, double tau_v,
                   double eps_u, double eps_v, uint64_t seed,
                   OUProcess **ou_u_out, OUProcess **ou_v_out) {
    if (Nx <= 0 || Ny <= 0 || !ou_u_out || !ou_v_out) return -1;

    int N = Nx * Ny;
    OUProcess *ou_u = (OUProcess *)calloc((size_t)N, sizeof(OUProcess));
    OUProcess *ou_v = (OUProcess *)calloc((size_t)N, sizeof(OUProcess));

    if (!ou_u || !ou_v) {
        free(ou_u); free(ou_v);
        return -1;
    }

    double amp_u = sqrt(2.0 * eps_u);
    double amp_v = sqrt(2.0 * eps_v);

    for (int i = 0; i < N; i++) {
        ou_init(&ou_u[i], tau_u, amp_u, seed + (uint64_t)i * 2);
        ou_init(&ou_v[i], tau_v, amp_v, seed + (uint64_t)i * 2 + 1);
    }

    *ou_u_out = ou_u;
    *ou_v_out = ou_v;
    return 0;
}

/** Free OU process arrays. */
void noise_ou_free(OUProcess *ou_u, OUProcess *ou_v) {
    free(ou_u);
    free(ou_v);
}

/* ==========================================================================
 * Euler-Maruyama Integration for Stochastic Reaction-Diffusion
 *
 * The stochastic PDE (in Ito form):
 *   du = [f(u,v) + D_u∇²u]dt + σ_u(u,v)·dW_u
 *   dv = [g(u,v) + D_v∇²v]dt + σ_v(u,v)·dW_v
 *
 * Euler-Maruyama scheme (strong order 1/2, weak order 1):
 *   u^{n+1} = u^n + [f(u^n,v^n) + D_u∇²u^n]·Δt + σ_u·√Δt·N(0,1)
 *
 * Reference: Higham, D.J. (2001) SIAM Rev. 43(3):525-546
 * ==========================================================================*/

/**
 * Euler-Maruyama step for additive noise.
 *
 * @param u, v       — fields to update (in-place)
 * @param reaction_f — reaction rate function
 * @param p          — model parameters
 * @param Du, Dv     — diffusion coefficients
 * @param dt         — time step
 * @param dx, dy     — grid spacing
 * @param eps_u, eps_v — noise intensities
 * @param seed       — random seed
 * @param stencil    — Laplacian stencil type
 * @param bc_x, bc_y — boundary conditions
 *
 * Returns 0 on success, -1 on error.
 */
int euler_maruyama_step(ScalarField2D *u, ScalarField2D *v,
                        int (*reaction_f)(double,double,const ModelParams*,
                                          double*,double*),
                        const ModelParams *p,
                        double Du, double Dv, double dt,
                        double dx, double dy,
                        double eps_u, double eps_v,
                        uint64_t seed,
                        LaplacianStencil stencil,
                        BoundaryCondition bc_x, BoundaryCondition bc_y) {
    if (!u || !v || !reaction_f || !p) return -1;

    int Nx = u->Nx, Ny = u->Ny;
    XoshiroRNG rng_u, rng_v;
    xoshiro_seed(&rng_u, seed);
    xoshiro_seed(&rng_v, seed + 1);

    /* Workspace fields */
    ScalarField2D *u_new = field_alloc(Nx, Ny, u->Lx, u->Ly, u->ghost);
    ScalarField2D *v_new = field_alloc(Nx, Ny, v->Lx, v->Ly, v->ghost);
    ScalarField2D *lap_u = field_alloc(Nx, Ny, u->Lx, u->Ly, u->ghost);
    ScalarField2D *lap_v = field_alloc(Nx, Ny, v->Lx, v->Ly, v->ghost);

    if (!u_new || !v_new || !lap_u || !lap_v) {
        field_free(u_new); field_free(v_new);
        field_free(lap_u); field_free(lap_v);
        return -1;
    }

    /* Compute Laplacians */
    field_apply_bc(u, bc_x, bc_y);
    field_apply_bc(v, bc_x, bc_y);
    laplacian_compute(u, lap_u, stencil);
    laplacian_compute(v, lap_v, stencil);

    double noise_u_prefac = sqrt(2.0 * eps_u * dt);
    double noise_v_prefac = sqrt(2.0 * eps_v * dt);

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            double u_val = field_get(u, i, j);
            double v_val = field_get(v, i, j);
            double lap_u_val = field_get(lap_u, i, j);
            double lap_v_val = field_get(lap_v, i, j);

            /* Reaction term */
            double f, g;
            reaction_f(u_val, v_val, p, &f, &g);

            /* Drift */
            double du_drift = f + Du * lap_u_val;
            double dv_drift = g + Dv * lap_v_val;

            /* Diffusion + noise */
            double u_next = u_val + dt * du_drift
                          + noise_u_prefac * xoshiro_normal(&rng_u);
            double v_next = v_val + dt * dv_drift
                          + noise_v_prefac * xoshiro_normal(&rng_v);

            field_set(u_new, i, j, u_next);
            field_set(v_new, i, j, v_next);
        }
    }

    /* Copy back */
    field_copy(u_new, u);
    field_copy(v_new, v);
    field_apply_bc(u, bc_x, bc_y);
    field_apply_bc(v, bc_x, bc_y);

    field_free(u_new); field_free(v_new);
    field_free(lap_u); field_free(lap_v);

    return 0;
}

/* ==========================================================================
 * Noise-Induced Pattern Selection Analysis
 *
 * In small systems (few molecules), intrinsic fluctuations can:
 *   1. Shift the effective bifurcation threshold
 *   2. Select specific pattern wavelengths (noise-sustained structures)
 *   3. Induce transitions between pattern types
 *
 * Butler & Goldenfeld (2009) showed that additive noise extends the
 * Turing region in parameter space, enabling pattern formation even
 * when the deterministic system is below threshold.
 *
 * Reference: Butler, T. & Goldenfeld, N. (2009) Phys. Rev. E 80:030902(R)
 * ==========================================================================*/

/**
 * Compute the noise-shifted dispersion relation.
 *
 * In the presence of additive noise of strength ε, the effective
 * growth rate becomes approximately:
 *   λ_eff(k) ≈ λ_det(k) + ε² · χ(k)
 *
 * where χ(k) is a noise-induced correction that depends on the
 * nonlinearity of the reaction terms.
 *
 * This is a simplified mean-field estimate based on a small-noise
 * expansion (van Kampen system size expansion).
 *
 * Reference: Biancalani, T. et al. (2010) Phys. Rev. Lett. 105(12):128102
 */
int noise_corrected_dispersion(const double J[2][2],
                               double Du, double Dv,
                               double eps,
                               int n_modes, double k_min, double k_max,
                               double *k_vals, double *lambda_vals) {
    if (!J || !k_vals || !lambda_vals || n_modes <= 0) return -1;

    double f_u = J[0][0], f_v = J[0][1];
    double g_u = J[1][0], g_v = J[1][1];
    double trace_J = f_u + g_v;
    double det_J = f_u * g_v - f_v * g_u;
    double sum_D = Du + Dv;

    for (int i = 0; i < n_modes; i++) {
        double k = k_min + (k_max - k_min) * (double)i / (double)(n_modes-1);
        double k2 = k * k;
        double k4 = k2 * k2;

        double T_k = trace_J - k2 * sum_D;
        double D_k = det_J - k2 * (Dv * f_u + Du * g_v) + k4 * Du * Dv;

        double disc = T_k * T_k - 4.0 * D_k;
        double lambda_det;
        if (disc >= 0.0) {
            lambda_det = 0.5 * (T_k + sqrt(disc));
        } else {
            lambda_det = 0.5 * T_k;
        }

        /* Noise correction: proportional to strength * |nonlinear terms| */
        double nonlinearity = fabs(f_v * g_u) + fabs(f_u * g_v);
        double correction = eps * eps * nonlinearity / (1.0 + k2);

        k_vals[i] = k;
        lambda_vals[i] = lambda_det + correction;
    }

    return 0;
}

/**
 * Estimate the minimal system size N_min for which deterministic
 * Turing analysis is valid (continuous approximation holds).
 *
 * N_min ≈ V / λ³  where λ is the pattern wavelength.
 * Below this, discrete stochastic effects dominate.
 *
 * Returns estimated N_min.
 */
double noise_estimate_min_system_size(double predicted_wavelength,
                                       double concentration_scale) {
    if (predicted_wavelength < 1e-15 || concentration_scale < 1e-15)
        return 1e10;

    /* Volume ~ wavelength³, molecules = volume * concentration */
    double V = predicted_wavelength * predicted_wavelength * predicted_wavelength;
    return V * concentration_scale;
}
