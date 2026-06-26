#include "oep_order_parameter.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==============================================================
 * Order Parameter Dynamics Implementation
 *
 * Implements the canonical order parameter equation:
 *   dξ/dt = αξ - βξ³ + √Q Γ(t)
 *
 * derived from the Ginzburg-Landau free energy:
 *   V(ξ) = -αξ²/2 + βξ⁴/4
 *
 * with the deterministic dynamics dξ/dt = -∂V/∂ξ.
 *
 * Near criticality (α → 0), the relaxation time diverges:
 *   τ = 1/α  (critical slowing down)
 *
 * References:
 *   Haken (1977) §5.1-5.3 — Derivation of order parameter equation
 *   Cross & Hohenberg (1993) §IV — Normal form theory
 *   Risken (1984) §5 — Fokker-Planck approach
 * ============================================================== */

/* ── Scalar Order Parameter ───────────────────────────── */

OEPScalarOP* oep_op_scalar_create(double alpha, double beta, double noise) {
    OEPScalarOP* op = calloc(1, sizeof(OEPScalarOP));
    if (!op) return NULL;
    op->value = 0.0;
    op->critical_value = 0.0;
    op->alpha = alpha;
    op->beta = beta > 0.0 ? beta : 1.0;  /* Ensure β > 0 for stability */
    op->noise_intensity = noise;
    if (fabs(alpha) > OEP_EPSILON)
        op->relaxation_time = 1.0 / fabs(alpha);
    else
        op->relaxation_time = 1e10;  /* Near critical → very slow */
    op->is_active = (alpha > 0.0);
    return op;
}

void oep_op_scalar_free(OEPScalarOP* op) { free(op); }

void oep_op_scalar_set_critical(OEPScalarOP* op, double pc, double p) {
    if (!op) return;
    op->alpha = (pc - p) / pc;  /* α ∝ (p_c - p) near criticality */
    op->critical_value = pc;
    double a = fabs(op->alpha);
    op->relaxation_time = (a > OEP_EPSILON) ? 1.0 / a : 1e10;
    op->is_active = (op->alpha > 0.0);
}

void oep_op_scalar_set_bifurcation(OEPScalarOP* op, double control_param,
                                    double critical_point, double sensitivity) {
    if (!op) return;
    op->alpha = sensitivity * (control_param - critical_point);
    op->critical_value = control_param;
    double a = fabs(op->alpha);
    op->relaxation_time = (a > OEP_EPSILON) ? 1.0 / a : 1e10;
    op->is_active = (op->alpha > 0.0);
}

double oep_op_scalar_rhs(double xi, double alpha, double beta, double noise) {
    /* dξ/dt = αξ - βξ³ + noise
     * This is the normal form of a supercritical pitchfork bifurcation. */
    return alpha * xi - beta * xi * xi * xi + noise;
}

double oep_op_scalar_drift(double xi, double alpha, double beta) {
    return alpha * xi - beta * pow(xi, 3.0);
}

double oep_op_scalar_potential(double xi, double alpha, double beta) {
    /* V(ξ) = -αξ²/2 + βξ⁴/4
     * The minima of V correspond to stable equilibria.
     * For α > 0: minima at ξ = ±√(α/β), maximum at ξ = 0
     * For α ≤ 0: single minimum at ξ = 0 */
    return -0.5 * alpha * xi * xi + 0.25 * beta * pow(xi, 4.0);
}

void oep_op_scalar_equilibria(double alpha, double beta,
                               double* eq0, double* eq1, double* eq2,
                               int* n_equilibria) {
    /* Solve αξ - βξ³ = 0 → ξ(α - βξ²) = 0
     * Solutions: ξ = 0, and if α/β > 0: ξ = ±√(α/β) */
    *eq0 = 0.0;
    if (alpha > 0.0 && beta > OEP_EPSILON) {
        double amp = sqrt(alpha / beta);
        *eq1 = amp;
        *eq2 = -amp;
        *n_equilibria = 3;
    } else if (alpha < 0.0) {
        *eq1 = 0.0;
        *eq2 = 0.0;
        *n_equilibria = 1;  /* Only ξ=0 is real (others are imaginary) */
    } else {
        *eq1 = 0.0;
        *eq2 = 0.0;
        *n_equilibria = 1;  /* At critical point, triple root at 0 */
    }
}

double oep_op_scalar_stability(double xi_eq, double alpha, double beta) {
    /* Linear stability: λ = ∂(dξ/dt)/∂ξ|_{\xi_eq}
     *   = α - 3βξ²_eq
     * For ξ_eq = 0: λ = α  → stable if α < 0
     * For ξ_eq = ±√(α/β): λ = α - 3β(α/β) = α - 3α = -2α  → stable if α > 0 */
    return alpha - 3.0 * beta * xi_eq * xi_eq;
}

bool oep_op_scalar_is_stable(double xi_eq, double alpha, double beta) {
    return oep_op_scalar_stability(xi_eq, alpha, beta) < -OEP_EPSILON;
}

double oep_op_scalar_relaxation_time(double alpha, double xi_eq) {
    /* τ = 1/|λ| where λ = α - 3βξ²_eq */
    double lam = alpha - 3.0 * 1.0 * xi_eq * xi_eq;  /* β=1 assumed for scaling */
    if (fabs(lam) < OEP_EPSILON) return 1e10;
    return 1.0 / fabs(lam);
}

double oep_op_scalar_susceptibility(double alpha, double beta, double xi_eq) {
    /* Susceptibility χ = ∂ξ/∂f|_ext
     * From dξ/dt = αξ - βξ³ + f_ext
     * At equilibrium: αξ_eq - βξ_eq³ + f_ext = 0
     * Differentiate: (α - 3βξ_eq²) dξ/df + 1 = 0
     * → dξ/df = -1/(α - 3βξ_eq²) = χ
     * At ξ_eq = 0: χ = -1/α, diverges as α → 0 */
    double denom = alpha - 3.0 * beta * xi_eq * xi_eq;
    if (fabs(denom) < OEP_EPSILON) return 1e10;
    return -1.0 / denom;
}

double oep_op_scalar_critical_exponent_susceptibility(double alpha) {
    /* Near critical point (α → 0 from above), χ ∝ 1/α
     * This gives γ = 1 (mean-field susceptibility exponent). */
    if (fabs(alpha) < OEP_EPSILON) return 1e10;
    return 1.0 / fabs(alpha);
}

double oep_op_scalar_correlation_length(double xi, double xi_sq_mean) {
    /* Correlation length proxy: ξ_corr ∝ 1/√(<ξ²> - <ξ>²)
     * The variance diverges at criticality. */
    double var = xi_sq_mean - xi * xi;
    if (var < OEP_EPSILON) return 1e10;
    return 1.0 / sqrt(var);
}

/* ── ODE Callback for Scalar OP ────────────────────────── */
typedef struct {
    double alpha;
    double beta;
} OEPScalarOPCtx;

static void oep_op_scalar_ode_f(double t, const double* y, double* dydt, void* ctx) {
    (void)t;
    OEPScalarOPCtx* c = (OEPScalarOPCtx*)ctx;
    dydt[0] = c->alpha * y[0] - c->beta * y[0] * y[0] * y[0];
}

void oep_op_scalar_evolve(OEPScalarOP* op, double dt, int n_steps,
                           double* trajectory) {
    if (!op || !trajectory || n_steps <= 0) return;
    OEPScalarOPCtx ctx;
    ctx.alpha = op->alpha;
    ctx.beta = op->beta;

    double y[1];
    y[0] = op->value;
    trajectory[0] = y[0];
    for (int i = 1; i < n_steps; i++) {
        double y_next[1];
        oep_ode_rk4_step((double)(i - 1) * dt, y, 1, dt,
                          oep_op_scalar_ode_f, &ctx, y_next);
        y[0] = y_next[0];
        trajectory[i] = y[0];
    }
}

void oep_op_scalar_evolve_stochastic(OEPScalarOP* op, double dt, int n_steps,
                                      double* trajectory) {
    /* Euler-Maruyama scheme for:
     *   dξ = (αξ - βξ³) dt + √(Q) dW
     * Discretization:
     *   ξ_{n+1} = ξ_n + (αξ_n - βξ_n³) dt + √(Q·dt) * N(0,1)
     *
     * Convergence: strong order 1/2, weak order 1.
     * Reference: Kloeden & Platen (1992) §10.2 */
    if (!op || !trajectory || n_steps <= 0) return;

    double xi = op->value;
    trajectory[0] = xi;
    double sqrt_Qdt = sqrt(op->noise_intensity * dt);

    for (int i = 1; i < n_steps; i++) {
        double drift = op->alpha * xi - op->beta * xi * xi * xi;
        double noise = sqrt_Qdt * oep_rng_normal(0.0, 1.0);
        xi = xi + drift * dt + noise;
        trajectory[i] = xi;
    }
    op->value = xi;
}

double oep_op_scalar_stationary_pdf(double xi, double alpha, double beta, double Q) {
    /* Stationary distribution of the Fokker-Planck equation:
     *   P_s(ξ) = (1/Z) exp(-2V(ξ)/Q)
     * with V(ξ) = -αξ²/2 + βξ⁴/4
     *
     * This is accurate for additive white noise with intensity Q.
     * Normalization Z is not computed here; returns unnormalized density. */
    double V = -0.5 * alpha * xi * xi + 0.25 * beta * xi * xi * xi * xi;
    return exp(-2.0 * V / Q);
}

double oep_op_scalar_most_probable(double alpha, double beta, double Q) {
    (void)Q;  /* Additive noise does not shift the mode */
    /* Mode of P_s(ξ) = argmin_ξ V(ξ)
     * Extrema: dV/dξ = -αξ + βξ³ = 0
     *   ξ = 0, and if α > 0: ξ = ±√(α/β)
     * For α > 0, V has minima at nonzero ξ (broken symmetry).
     * For α ≤ 0, minimum at ξ = 0.
     *
     * Noise Q shifts the effective potential but the mode is unchanged
     * for additive noise. */
    if (alpha > OEP_EPSILON && beta > OEP_EPSILON)
        return sqrt(alpha / beta);
    return 0.0;
}

/* ── Complex Order Parameter ───────────────────────────── */

OEPComplexOP* oep_op_complex_create(double growth_rate, double frequency,
                                     double nonlinear_coeff, double noise_real,
                                     double noise_imag) {
    OEPComplexOP* op = calloc(1, sizeof(OEPComplexOP));
    if (!op) return NULL;
    op->amplitude = 0.01;  /* Small initial perturbation */
    op->phase = 0.0;
    op->frequency = frequency;
    op->growth_rate = growth_rate;
    op->nonlinear_coeff = nonlinear_coeff;
    op->noise_real = noise_real;
    op->noise_imag = noise_imag;
    return op;
}

void oep_op_complex_free(OEPComplexOP* op) { free(op); }

void oep_op_complex_rhs(double amplitude, double phase, double sigma_re,
                         double sigma_im, double g_re, double g_im,
                         double* dA_re, double* dA_imag) {
    /* Complex Ginzburg-Landau (CGL) equation:
     *   ∂A/∂t = σA - g|A|²A
     * with A = r·exp(iθ), σ = σ_re + i·σ_im, g = g_re + i·g_im
     *
     * In polar coordinates:
     *   dr/dt = σ_re·r - g_re·r³
     *   dθ/dt = σ_im - g_im·r²
     *
     * Then transform back to Cartesian for the complex RHS. */
    double r = amplitude;
    double drdt = sigma_re * r - g_re * r * r * r;
    double dthetadt = sigma_im - g_im * r * r;

    /* Convert polar derivatives to Cartesian */
    double cos_th = cos(phase);
    double sin_th = sin(phase);

    *dA_re = drdt * cos_th - r * dthetadt * sin_th;
    *dA_imag = drdt * sin_th + r * dthetadt * cos_th;
}

void oep_op_complex_evolve(OEPComplexOP* op, double dt, int n_steps,
                            double* re_traj, double* im_traj) {
    /* Evolve the complex order parameter using semi-analytical integration.
     * The amplitude and phase equations are decoupled at leading order,
     * allowing efficient integration. */
    if (!op || !re_traj || !im_traj || n_steps <= 0) return;

    double r = op->amplitude;
    double theta = op->phase;
    double sigma_re = op->growth_rate;
    double g_re = op->nonlinear_coeff;

    /* For simplicity, assume g_im = 0 and σ_im = ω (pure frequency) */
    double sigma_im = op->frequency;
    double g_im = 0.0;

    re_traj[0] = r * cos(theta);
    im_traj[0] = r * sin(theta);

    for (int i = 1; i < n_steps; i++) {
        /* Amplitude: dr/dt = σ_re·r - g_re·r³ (RK4) */
        double k1r = dt * (sigma_re * r - g_re * r * r * r);
        double r_mid1 = r + 0.5 * k1r;
        double k2r = dt * (sigma_re * r_mid1 - g_re * r_mid1 * r_mid1 * r_mid1);
        double r_mid2 = r + 0.5 * k2r;
        double k3r = dt * (sigma_re * r_mid2 - g_re * r_mid2 * r_mid2 * r_mid2);
        double r_end = r + k3r;
        double k4r = dt * (sigma_re * r_end - g_re * r_end * r_end * r_end);
        r = r + (k1r + 2.0 * k2r + 2.0 * k3r + k4r) / 6.0;
        if (r < 0.0) r = 0.0;

        /* Phase: dθ/dt = σ_im - g_im·r² */
        theta += dt * (sigma_im - g_im * r * r);

        /* Wrap phase */
        while (theta > 2.0 * OEP_PI) theta -= 2.0 * OEP_PI;
        while (theta < 0.0) theta += 2.0 * OEP_PI;

        re_traj[i] = r * cos(theta);
        im_traj[i] = r * sin(theta);
    }
    op->amplitude = r;
    op->phase = theta;
}

double oep_op_complex_equilibrium_amplitude(double sigma_re, double g_re) {
    /* Equilibrium: dr/dt = σ_re·r - g_re·r³ = 0
     *   r = 0 (trivial) or r = √(σ_re / g_re) (for σ_re > 0, g_re > 0)
     * This is the Hopf bifurcation normal form. */
    if (sigma_re > OEP_EPSILON && g_re > OEP_EPSILON)
        return sqrt(sigma_re / g_re);
    return 0.0;
}

double oep_op_complex_phase_velocity(double sigma_im, double g_im, double amplitude) {
    /* Phase velocity on the limit cycle:
     *   Ω = σ_im - g_im·r²
     * For r_eq = √(σ_re / g_re):
     *   Ω_eq = σ_im - g_im·(σ_re / g_re) */
    return sigma_im - g_im * amplitude * amplitude;
}

/* ── Multi-component Order Parameter ────────────────────── */

OEPMultiOP* oep_op_multi_create(int n_components) {
    if (n_components <= 0 || n_components > OEP_MAX_DIM) return NULL;
    OEPMultiOP* op = calloc(1, sizeof(OEPMultiOP));
    if (!op) return NULL;
    op->n_components = n_components;
    op->values = oep_vec_create(n_components);
    op->coupling = oep_mat_create(n_components, n_components);
    op->linear_coeffs = oep_vec_create(n_components);
    op->cubic_coeffs = oep_vec_create(n_components);
    if (!op->values || !op->coupling || !op->linear_coeffs || !op->cubic_coeffs) {
        oep_op_multi_free(op);
        return NULL;
    }
    return op;
}

void oep_op_multi_free(OEPMultiOP* op) {
    if (!op) return;
    oep_vec_free(op->values);
    oep_mat_free(op->coupling);
    oep_vec_free(op->linear_coeffs);
    oep_vec_free(op->cubic_coeffs);
    free(op);
}

void oep_op_multi_set_linear(OEPMultiOP* op, int i, double alpha) {
    if (op && i >= 0 && i < op->n_components)
        oep_vec_set(op->linear_coeffs, i, alpha);
}

void oep_op_multi_set_cubic(OEPMultiOP* op, int i, double beta) {
    if (op && i >= 0 && i < op->n_components)
        oep_vec_set(op->cubic_coeffs, i, beta);
}

void oep_op_multi_set_coupling(OEPMultiOP* op, int i, int j, double val) {
    if (op && i >= 0 && i < op->n_components && j >= 0 && j < op->n_components)
        oep_mat_set(op->coupling, i, j, val);
}

void oep_op_multi_rhs(const OEPMultiOP* op, const OEPVector* state,
                       OEPVector* deriv) {
    /* dξᵢ/dt = αᵢξᵢ - βᵢξᵢ³ + Σⱼ C_{ij}ξⱼ
     * The coupling term allows interaction between order parameter components,
     * which can produce mode competition, synchronization, or sequential ordering. */
    if (!op || !state || !deriv) return;
    int n = op->n_components;

    /* Compute C * state */
    OEPVector* coupled = oep_mat_vec_mul(op->coupling, state);
    if (!coupled) return;

    for (int i = 0; i < n; i++) {
        double xi = oep_vec_get(state, i);
        double alpha = oep_vec_get(op->linear_coeffs, i);
        double beta = oep_vec_get(op->cubic_coeffs, i);
        double coupling_term = oep_vec_get(coupled, i);
        deriv->data[i] = alpha * xi - beta * xi * xi * xi + coupling_term;
    }
    oep_vec_free(coupled);
}

/* Multi-component ODE callback */
typedef struct {
    OEPMultiOP* op;
} OEPMultiOPCtx;

static void oep_op_multi_ode_f(double t, const double* y, double* dydt, void* ctx) {
    (void)t;
    OEPMultiOPCtx* c = (OEPMultiOPCtx*)ctx;
    int n = c->op->n_components;
    /* Build temporary vectors */
    OEPVector state_vec;
    OEPVector deriv_vec;
    state_vec.data = (double*)y;  /* Const-correctness hack for internal use */
    state_vec.size = n;
    deriv_vec.data = dydt;
    deriv_vec.size = n;
    oep_op_multi_rhs(c->op, &state_vec, &deriv_vec);
}

void oep_op_multi_evolve(OEPMultiOP* op, double dt, int n_steps,
                          OEPMatrix* trajectory) {
    /* Each row of trajectory holds the state at a time step.
     * trajectory must be pre-allocated with n_steps rows and n_components cols. */
    if (!op || !trajectory || n_steps <= 0) return;
    int n = op->n_components;
    OEPMultiOPCtx ctx;
    ctx.op = op;

    double** y_traj = malloc((size_t)n_steps * sizeof(double*));
    if (!y_traj) return;
    for (int i = 0; i < n_steps; i++)
        y_traj[i] = trajectory->data + i * n;

    double* y0 = malloc((size_t)n * sizeof(double));
    if (!y0) { free(y_traj); return; }
    for (int i = 0; i < n; i++) y0[i] = oep_vec_get(op->values, i);

    oep_ode_rk4_multi(0.0, y0, n, n_steps, dt, oep_op_multi_ode_f, &ctx, y_traj);

    free(y0);
    free(y_traj);
}

void oep_op_multi_normal_modes(const OEPMultiOP* op, OEPVector* eigenvalues,
                                 OEPMatrix* eigenvectors) {
    /* Compute normal modes of the linearized system:
     *   dξ/dt = Λ·ξ + nonlinear terms
     * where Λ_ij = α_i·δ_ij + C_ij
     *
     * Diagonalize Λ using QR algorithm (simplified for small matrices).
     * Reference: Golub & Van Loan (2013) §8.3. */
    if (!op || !eigenvalues || !eigenvectors) return;
    int n = op->n_components;

    /* Build the linearized matrix Λ */
    OEPMatrix* Lambda = oep_mat_create(n, n);
    if (!Lambda) return;

    /* Λ = diag(α) + C */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double val = oep_mat_get(op->coupling, i, j);
            if (i == j) val += oep_vec_get(op->linear_coeffs, i);
            oep_mat_set(Lambda, i, j, val);
        }
    }

    /* Compute eigenvalues via direct QR iteration for small matrices.
     * For symmetric Λ: eigenvalues are real; for asymmetric: may be complex.
     * We use a simplified approach: tridiagonalize via Householder, then QR. */
    if (n == 1) {
        oep_vec_set(eigenvalues, 0, oep_mat_get(Lambda, 0, 0));
        oep_mat_set(eigenvectors, 0, 0, 1.0);
        oep_mat_free(Lambda);
        return;
    }

    /* For n ≤ OEP_MAX_DIM: use power iteration + deflation.
     * This gives the dominant eigenvalues and vectors. */
    double* A = Lambda->data;
    double* v_work = malloc((size_t)(3 * n) * sizeof(double));
    if (!v_work) { oep_mat_free(Lambda); return; }

    for (int k = 0; k < n; k++) {
        /* Power iteration for eigenvalue k */
        double* v = v_work;
        /* Initialize v with random values */
        for (int i = 0; i < n; i++) v[i] = oep_rng_normal(0.0, 1.0);
        /* Normalize */
        double norm = 0.0;
        for (int i = 0; i < n; i++) norm += v[i] * v[i];
        norm = sqrt(norm);
        if (norm > 0.0) for (int i = 0; i < n; i++) v[i] /= norm;

        double lambda_old = 0.0, lambda_new = 0.0;
        for (int iter = 0; iter < 100; iter++) {
            /* w = A * v */
            double* w = v_work + n;
            for (int i = 0; i < n; i++) {
                w[i] = 0.0;
                for (int j = 0; j < n; j++) w[i] += A[i * n + j] * v[j];
            }
            /* λ = v^T * w (Rayleigh quotient) */
            lambda_new = 0.0;
            norm = 0.0;
            for (int i = 0; i < n; i++) {
                lambda_new += v[i] * w[i];
                norm += w[i] * w[i];
            }
            norm = sqrt(norm);
            if (norm < OEP_EPSILON) break;
            for (int i = 0; i < n; i++) v[i] = w[i] / norm;

            if (fabs(lambda_new - lambda_old) < OEP_EPSILON) break;
            lambda_old = lambda_new;
        }

        oep_vec_set(eigenvalues, k, lambda_new);
        for (int i = 0; i < n; i++)
            oep_mat_set(eigenvectors, i, k, v[i]);

        /* Deflate A: A ← A - λ·v·v^T */
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                A[i * n + j] -= lambda_new * v[i] * v[j];
    }

    free(v_work);
    oep_mat_free(Lambda);
}

int oep_op_multi_adiabatic_reduce(OEPMultiOP* op, double timescale_threshold,
                                   OEPVector* remaining_indices) {
    /* Adiabatic elimination: remove fast-relaxing modes (damping > threshold).
     * Modes with large negative eigenvalues relax quickly and can be
     * adiabatically eliminated (slaved to the slow modes).
     *
     * Returns number of remaining (slow) modes. */
    if (!op || !remaining_indices) return 0;
    int n = op->n_components;

    OEPVector* eigenvalues = oep_vec_create(n);
    OEPMatrix* eigenvectors = oep_mat_create(n, n);
    if (!eigenvalues || !eigenvectors) {
        oep_vec_free(eigenvalues);
        oep_mat_free(eigenvectors);
        return 0;
    }

    oep_op_multi_normal_modes(op, eigenvalues, eigenvectors);

    int count = 0;
    for (int i = 0; i < n; i++) {
        double lam = oep_vec_get(eigenvalues, i);
        /* Keep modes with |λ| < threshold (slow modes)
         * or with λ > 0 (unstable/active modes — order parameter candidates) */
        if (fabs(lam) < timescale_threshold || lam > -OEP_EPSILON) {
            oep_vec_set(remaining_indices, count, (double)i);
            count++;
        }
    }

    oep_vec_free(eigenvalues);
    oep_mat_free(eigenvectors);
    return count;
}

/* ── Order Parameter Identification from Data ───────────── */

OEPIdentifiedOP** oep_op_identify_from_data(const OEPMatrix* data,
                                             int* n_identified) {
    /* Identify candidate order parameters from multivariate time series.
     *
     * Method: Principal Component Analysis (PCA) on the data matrix.
     * Dimensions with large variance and slow autocorrelation decay
     * are candidates for order parameters.
     *
     * Reference: Haken (1983) §8 — "Information compression by order parameters" */
    if (!data || !n_identified || data->rows < 3 || data->cols < 1) {
        if (n_identified) *n_identified = 0;
        return NULL;
    }
    int n_obs = data->rows;
    int n_vars = data->cols;

    /* Compute PCA */
    OEPVector* eigenvalues = oep_vec_create(n_vars);
    OEPMatrix* loadings = oep_mat_create(n_vars, n_vars);
    OEPVector* explained_var = oep_vec_create(n_vars);
    oep_op_pca(data, eigenvalues, loadings, explained_var);

    /* Allocate result array (max: all components could be OPs) */
    OEPIdentifiedOP** results = calloc((size_t)n_vars, sizeof(OEPIdentifiedOP*));
    if (!results) {
        oep_vec_free(eigenvalues);
        oep_mat_free(loadings);
        oep_vec_free(explained_var);
        *n_identified = 0;
        return NULL;
    }

    int count = 0;
    double total_var = 0.0;
    for (int i = 0; i < n_vars; i++) total_var += oep_vec_get(eigenvalues, i);

    /* Threshold: components explaining >5% of variance */
    for (int i = 0; i < n_vars; i++) {
        double eig = oep_vec_get(eigenvalues, i);
        double frac = (total_var > OEP_EPSILON) ? eig / total_var : 0.0;
        if (frac < 0.05) continue;

        OEPIdentifiedOP* op = calloc(1, sizeof(OEPIdentifiedOP));
        if (!op) continue;
        op->significance = frac;
        op->eigenvalue = eig;
        op->direction = oep_vec_create(n_vars);
        if (!op->direction) { free(op); continue; }
        for (int j = 0; j < n_vars; j++)
            oep_vec_set(op->direction, j, oep_mat_get(loadings, j, i));

        /* Estimate relaxation time from autocorrelation of the component */
        double* component = malloc((size_t)n_obs * sizeof(double));
        if (component) {
            for (int t = 0; t < n_obs; t++) {
                double score = 0.0;
                for (int j = 0; j < n_vars; j++)
                    score += oep_mat_get(data, t, j) * oep_mat_get(loadings, j, i);
                component[t] = score;
            }
            double ac1 = oep_ts_autocorrelation(component, n_obs, 1);
            /* τ = -1/log(ac1) for AR(1) process */
            if (ac1 > 0.0 && ac1 < 1.0 - OEP_EPSILON)
                op->relaxation_time = -1.0 / log(ac1);
            else
                op->relaxation_time = 1.0;
            free(component);
        }

        op->nonlinear_coeff = 1.0;  /* Default, to be refined */
        op->is_critical = (frac > 0.3);  /* High variance fraction → near critical */

        results[count] = op;
        count++;
    }

    oep_vec_free(eigenvalues);
    oep_mat_free(loadings);
    oep_vec_free(explained_var);
    *n_identified = count;
    return results;
}

void oep_op_id_free(OEPIdentifiedOP* op) {
    if (!op) return;
    oep_vec_free(op->direction);
    free(op);
}

void oep_op_pca(const OEPMatrix* data, OEPVector* eigenvalues,
                OEPMatrix* loadings, OEPVector* explained_variance) {
    /* Principal Component Analysis via covariance eigendecomposition.
     *
     * Steps:
     *   1. Center the data
     *   2. Compute covariance matrix C = X^T·X / (n-1)
     *   3. Eigendecomposition via power iteration with deflation
     *
     * Complexity: O(n_obs·n_vars² + n_vars³)
     * Reference: Jolliffe (2002) "Principal Component Analysis" */
    if (!data || !eigenvalues || !loadings) return;
    int n_obs = data->rows;
    int n_vars = data->cols;
    if (n_obs < 2 || n_vars < 2) return;

    /* Center the data */
    double* X = malloc((size_t)(n_obs * n_vars) * sizeof(double));
    if (!X) return;
    memcpy(X, data->data, (size_t)(n_obs * n_vars) * sizeof(double));
    for (int j = 0; j < n_vars; j++) {
        double mean = 0.0;
        for (int i = 0; i < n_obs; i++) mean += X[i * n_vars + j];
        mean /= (double)n_obs;
        for (int i = 0; i < n_obs; i++) X[i * n_vars + j] -= mean;
    }

    /* Covariance matrix: C = X^T·X / (n-1), size: n_vars × n_vars */
    double* C = calloc((size_t)(n_vars * n_vars), sizeof(double));
    if (!C) { free(X); return; }
    for (int i = 0; i < n_obs; i++) {
        for (int p = 0; p < n_vars; p++) {
            double xip = X[i * n_vars + p];
            for (int q = p; q < n_vars; q++) {
                C[p * n_vars + q] += xip * X[i * n_vars + q];
            }
        }
    }
    double scale = 1.0 / (double)(n_obs - 1);
    for (int p = 0; p < n_vars; p++) {
        for (int q = p; q < n_vars; q++) {
            C[p * n_vars + q] *= scale;
            C[q * n_vars + p] = C[p * n_vars + q];  /* Symmetrize */
        }
    }

    /* Power iteration with deflation for eigenvalues/vectors */
    double* A = C;
    double* v_work = malloc((size_t)(2 * n_vars) * sizeof(double));
    if (!v_work) { free(C); free(X); return; }
    double total_eig = 0.0;

    for (int k = 0; k < n_vars; k++) {
        double* v = v_work;
        /* Initialize with small random */
        for (int i = 0; i < n_vars; i++) v[i] = oep_rng_normal(0.0, 1.0);
        double norm = 0.0;
        for (int i = 0; i < n_vars; i++) norm += v[i] * v[i];
        norm = sqrt(norm);
        if (norm > 0.0) for (int i = 0; i < n_vars; i++) v[i] /= norm;

        double lam_old = 0.0, lam_new = 0.0;
        for (int iter = 0; iter < 100; iter++) {
            double* w = v_work + n_vars;
            for (int i = 0; i < n_vars; i++) {
                w[i] = 0.0;
                for (int j = 0; j < n_vars; j++)
                    w[i] += A[j * n_vars + i] * v[j];  /* A^T·v for symmetric A */
            }
            lam_new = 0.0;
            norm = 0.0;
            for (int i = 0; i < n_vars; i++) {
                lam_new += v[i] * w[i];
                norm += w[i] * w[i];
            }
            norm = sqrt(norm);
            if (norm < OEP_EPSILON) break;
            for (int i = 0; i < n_vars; i++) v[i] = w[i] / norm;
            if (fabs(lam_new - lam_old) < OEP_EPSILON) break;
            lam_old = lam_new;
        }
        if (lam_new < 0.0) lam_new = 0.0;  /* Covariance eigenvalues ≥ 0 */
        total_eig += lam_new;
        oep_vec_set(eigenvalues, k, lam_new);
        for (int i = 0; i < n_vars; i++)
            oep_mat_set(loadings, i, k, v[i]);

        /* Deflate */
        for (int i = 0; i < n_vars; i++)
            for (int j = 0; j < n_vars; j++)
                A[j * n_vars + i] -= lam_new * v[j] * v[i];
    }

    /* Explained variance fraction */
    if (total_eig > OEP_EPSILON) {
        for (int k = 0; k < n_vars; k++)
            oep_vec_set(explained_variance, k,
                        oep_vec_get(eigenvalues, k) / total_eig);
    }

    free(v_work);
    free(C);
    free(X);
}

bool oep_op_is_valid_order_parameter(const OEPIdentifiedOP* op,
                                      double var_threshold,
                                      double relaxation_ratio) {
    (void)relaxation_ratio;  /* Reserved for multi-mode comparison */
    /* Validation criteria (Haken, 1983):
     *   1. Large variance fraction (significance > var_threshold)
     *   2. Slow relaxation compared to other modes
     *   3. Symmetry-breaking near criticality */
    if (!op) return false;
    return (op->significance >= var_threshold);
}

double oep_op_circular_causality(const OEPMatrix* data,
                                  const OEPVector* op_direction,
                                  const OEPMatrix* microscopic) {
    /* Circular causality measure: quantifies the bidirectional influence
     * between the order parameter (macro) and microscopic variables.
     *
     * Reference: Haken (1983) §7.5 — "The concept of circular causality"
     *
     * Returns: correlation between OP state and micro-state reconstruction.
     *
     * The idea: if the OP both emerges from and constrains the micro variables,
     * then the micro variables should be predictable from the OP alone.
     * The strength of this predictability is the circular causality measure. */
    if (!data || !op_direction) return 0.0;
    int n_obs = data->rows;
    int n_vars = data->cols;
    (void)microscopic;  /* Reserved for future use */

    if (n_obs < 3 || n_vars < 2) return 0.0;

    /* Project data onto OP direction */
    double* op_ts = malloc((size_t)n_obs * sizeof(double));
    if (!op_ts) return 0.0;
    for (int t = 0; t < n_obs; t++) {
        op_ts[t] = 0.0;
        for (int j = 0; j < n_vars; j++)
            op_ts[t] += oep_mat_get(data, t, j) * oep_vec_get(op_direction, j);
    }

    /* Measure mutual information between OP and first micro variable */
    double* micro_ts = malloc((size_t)n_obs * sizeof(double));
    if (!micro_ts) { free(op_ts); return 0.0; }
    for (int t = 0; t < n_obs; t++)
        micro_ts[t] = oep_mat_get(data, t, 0);

    /* Use correlation as proxy for circular causality */
    double m_op = oep_ts_mean(op_ts, n_obs);
    double m_mic = oep_ts_mean(micro_ts, n_obs);
    double cov = 0.0, var_op = 0.0, var_mic = 0.0;
    for (int t = 0; t < n_obs; t++) {
        double d_op = op_ts[t] - m_op;
        double d_mic = micro_ts[t] - m_mic;
        cov += d_op * d_mic;
        var_op += d_op * d_op;
        var_mic += d_mic * d_mic;
    }
    double corr = (var_op > OEP_EPSILON && var_mic > OEP_EPSILON)
                  ? cov / sqrt(var_op * var_mic) : 0.0;

    free(op_ts);
    free(micro_ts);
    return fabs(corr);
}
