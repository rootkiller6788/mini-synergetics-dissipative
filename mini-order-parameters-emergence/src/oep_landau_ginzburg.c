#include "oep_landau_ginzburg.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==============================================================
 * Landau-Ginzburg Theory Implementation
 *
 * Implements the spatially-extended theory of phase transitions
 * with gradient terms for pattern formation.
 *
 * Spatial TDGL equation (Model A):
 *   ∂ξ/∂t = -Γ [a(p)·ξ + b·ξ³ - c·∇²ξ] + η(r, t)
 *
 * Domain wall solutions:
 *   ξ(x) = ξ₀·tanh(x / (√2·ξ))
 *
 * Superconductivity applications (Ginzburg-Landau equations):
 *   αΨ + β|Ψ|²Ψ + (1/2m*)(ħ∇/i - e*A)²Ψ = 0
 *
 * References:
 *   Ginzburg & Landau (1950) — Superconductivity
 *   Hohenberg & Halperin (1977) — Dynamic critical phenomena
 *   Bray (1994) — Phase ordering dynamics
 * ============================================================== */

/* ── Spatial Field ─────────────────────────────────────────── */

OEPSpatialField* oep_field_create(int nx, int ny, double Lx, double Ly) {
    if (nx <= 0 || ny <= 0 || Lx <= 0.0 || Ly <= 0.0) return NULL;
    OEPSpatialField* field = calloc(1, sizeof(OEPSpatialField));
    if (!field) return NULL;
    field->nx = nx;
    field->ny = ny;
    field->dx = Lx / (double)nx;
    field->dy = Ly / (double)ny;
    field->Lx = Lx;
    field->Ly = Ly;
    field->data = calloc((size_t)(nx * ny), sizeof(double));
    if (!field->data) { free(field); return NULL; }
    return field;
}

void oep_field_free(OEPSpatialField* field) {
    if (!field) return;
    free(field->data);
    free(field);
}

void oep_field_set(OEPSpatialField* field, int i, int j, double val) {
    if (field && i >= 0 && i < field->nx && j >= 0 && j < field->ny)
        field->data[i + j * field->nx] = val;
}

double oep_field_get(const OEPSpatialField* field, int i, int j) {
    if (!field || i < 0 || i >= field->nx || j < 0 || j >= field->ny) return 0.0;
    return field->data[i + j * field->nx];
}

void oep_field_zero(OEPSpatialField* field) {
    if (field) memset(field->data, 0, (size_t)(field->nx * field->ny) * sizeof(double));
}

void oep_field_set_random(OEPSpatialField* field, double amplitude) {
    if (!field) return;
    for (int k = 0; k < field->nx * field->ny; k++)
        field->data[k] = amplitude * (oep_rng_uniform() - 0.5) * 2.0;
}

/* ── Spatial Derivatives (centered finite differences, periodic BC) ── */

double oep_field_laplacian(const OEPSpatialField* field, int i, int j) {
    /* 5-point stencil for ∇² on 2D periodic grid:
     *   ∇²ξ = (ξ_{i+1,j} + ξ_{i-1,j} + ξ_{i,j+1} + ξ_{i,j-1} - 4ξ_{i,j}) / dx² */
    if (!field) return 0.0;
    int nx = field->nx, ny = field->ny;
    double dx2 = field->dx * field->dx;
    double dy2 = field->dy * field->dy;
    if (dx2 < OEP_EPSILON || dy2 < OEP_EPSILON) return 0.0;

    int ip = (i + 1) % nx;
    int im = (i - 1 + nx) % nx;
    int jp = (j + 1) % ny;
    int jm = (j - 1 + ny) % ny;

    double c = field->data[i + j * nx];
    double d2x = (field->data[ip + j * nx] + field->data[im + j * nx] - 2.0 * c) / dx2;
    double d2y = (field->data[i + jp * nx] + field->data[i + jm * nx] - 2.0 * c) / dy2;

    return d2x + d2y;
}

double oep_field_grad_x(const OEPSpatialField* field, int i, int j) {
    /* Centered difference: ∂ξ/∂x = (ξ_{i+1,j} - ξ_{i-1,j}) / (2·dx) */
    if (!field) return 0.0;
    int nx = field->nx;
    int ip = (i + 1) % nx;
    int im = (i - 1 + nx) % nx;
    return (field->data[ip + j * nx] - field->data[im + j * nx]) / (2.0 * field->dx);
}

double oep_field_grad_y(const OEPSpatialField* field, int i, int j) {
    if (!field) return 0.0;
    int ny = field->ny;
    int jp = (j + 1) % ny;
    int jm = (j - 1 + ny) % ny;
    return (field->data[i + jp * field->nx] - field->data[i + jm * field->nx])
           / (2.0 * field->dy);
}

double oep_field_grad_sq(const OEPSpatialField* field, int i, int j) {
    double gx = oep_field_grad_x(field, i, j);
    double gy = oep_field_grad_y(field, i, j);
    return gx * gx + gy * gy;
}

double oep_field_spatial_mean(const OEPSpatialField* field) {
    if (!field) return 0.0;
    double s = 0.0;
    int N = field->nx * field->ny;
    for (int k = 0; k < N; k++) s += field->data[k];
    return s / (double)N;
}

double oep_field_spatial_variance(const OEPSpatialField* field) {
    if (!field) return 0.0;
    double mu = oep_field_spatial_mean(field);
    double s = 0.0;
    int N = field->nx * field->ny;
    for (int k = 0; k < N; k++) {
        double d = field->data[k] - mu;
        s += d * d;
    }
    return s / (double)N;
}

double oep_field_free_energy(const OEPSpatialField* field,
                              double a, double b, double c,
                              double dx, double dy) {
    /* Compute the discretized Landau-Ginzburg free energy:
     *   F[ξ] = Σ_{i,j} [a·ξ_{ij}²/2 + b·ξ_{ij}⁴/4 + c/2·|∇ξ_{ij}|²] · dx·dy
     *
     * The gradient term penalizes spatial variations, leading to
     * smooth domain structures and a well-defined correlation length. */
    if (!field) return 0.0;
    double F = 0.0;
    for (int j = 0; j < field->ny; j++) {
        for (int i = 0; i < field->nx; i++) {
            double xi = oep_field_get(field, i, j);
            double local = 0.5 * a * xi * xi + 0.25 * b * xi * xi * xi * xi;
            local += 0.5 * c * oep_field_grad_sq(field, i, j);
            F += local;
        }
    }
    return F * dx * dy;
}

/* ── Landau-Ginzburg System ────────────────────────────────── */

OEPLandauGinzburg* oep_lg_create(int nx, int ny, double Lx, double Ly,
                                  double a, double b, double c, double Gamma) {
    if (nx <= 0 || ny <= 0 || Lx <= 0.0 || Ly <= 0.0) return NULL;
    OEPLandauGinzburg* lg = calloc(1, sizeof(OEPLandauGinzburg));
    if (!lg) return NULL;
    lg->nx = nx;
    lg->ny = ny;
    lg->dx = Lx / (double)nx;
    lg->dy = Ly / (double)ny;
    lg->Lx = Lx;
    lg->Ly = Ly;
    lg->a = a;
    lg->b = b;
    lg->c = c;
    lg->Gamma = Gamma;
    lg->noise_amplitude = 0.0;
    return lg;
}

void oep_lg_free(OEPLandauGinzburg* lg) { free(lg); }

void oep_lg_set_noise(OEPLandauGinzburg* lg, double amplitude) {
    if (lg) lg->noise_amplitude = amplitude;
}

double oep_lg_local_rhs(const OEPLandauGinzburg* lg,
                         const OEPSpatialField* field,
                         int i, int j) {
    /* Time-dependent Ginzburg-Landau RHS at grid point (i, j):
     *   ∂ξ/∂t = -Γ·[a·ξ + b·ξ³ - c·∇²ξ] + η
     *
     * This is Model A in the Hohenberg-Halperin classification,
     * appropriate for non-conserved order parameters. */
    if (!lg || !field) return 0.0;
    double xi = oep_field_get(field, i, j);
    double lap = oep_field_laplacian(field, i, j);
    double deterministic = -lg->Gamma * (lg->a * xi + lg->b * xi * xi * xi - lg->c * lap);
    double noise = lg->noise_amplitude * oep_rng_normal(0.0, 1.0);
    return deterministic + noise;
}

void oep_lg_step(OEPLandauGinzburg* lg, OEPSpatialField* field, double dt) {
    /* One explicit Euler step of TDGL.
     *   ξ_{n+1} = ξ_n + dt·(∂ξ/∂t)
     *
     * Stability condition: dt < dx²/(2·Γ·c) (diffusive CFL).
     * For explicit Euler, this is the limiting factor. */
    if (!lg || !field) return;

    double* new_data = malloc((size_t)(field->nx * field->ny) * sizeof(double));
    if (!new_data) return;

    for (int j = 0; j < field->ny; j++) {
        for (int i = 0; i < field->nx; i++) {
            double rhs = oep_lg_local_rhs(lg, field, i, j);
            int idx = i + j * field->nx;
            new_data[idx] = field->data[idx] + dt * rhs;
        }
    }

    memcpy(field->data, new_data, (size_t)(field->nx * field->ny) * sizeof(double));
    free(new_data);
}

void oep_lg_evolve(OEPLandauGinzburg* lg, OEPSpatialField* field,
                    double dt, int n_steps,
                    double* mean_trajectory,
                    double* variance_trajectory) {
    if (!lg || !field || n_steps <= 0) return;

    for (int step = 0; step < n_steps; step++) {
        if (mean_trajectory)
            mean_trajectory[step] = oep_field_spatial_mean(field);
        if (variance_trajectory)
            variance_trajectory[step] = oep_field_spatial_variance(field);
        oep_lg_step(lg, field, dt);
    }
}

/* ── Domain Wall Structure ─────────────────────────────────── */

void oep_lg_domain_wall_1d(OEPSpatialField* field,
                            double x_center, double a, double b,
                            double c, int nx) {
    /* 1D domain wall (kink) profile:
     *   ξ(x) = ξ₀·tanh((x - x_center) / (√2·ξ))
     * where ξ₀ = √(-a/b) is the equilibrium amplitude (a < 0),
     * and ξ = √(c/|a|) is the correlation length.
     *
     * This profile connects ξ = -ξ₀ (left) to ξ = +ξ₀ (right),
     * with a transition region of width ~ ξ. */
    if (!field || b <= 0.0 || nx <= 0) return;
    double a_abs = fabs(a);
    double xi0 = (a < 0.0) ? sqrt(a_abs / b) : 0.1;
    double xi_corr = (a_abs > OEP_EPSILON) ? sqrt(c / a_abs) : 1.0;
    double width = sqrt(2.0) * xi_corr;

    for (int i = 0; i < nx && i < field->nx; i++) {
        double x = (double)i * field->dx;
        double arg = (x - x_center) / width;
        for (int j = 0; j < field->ny; j++)
            oep_field_set(field, i, j, xi0 * tanh(arg));
    }
}

double oep_lg_domain_wall_energy(double a, double b, double c) {
    /* Domain wall surface tension (energy per unit area):
     *   σ = (2√2/3) · |a|^{3/2} · √c / b
     *
     * The domain wall energy controls the coarsening dynamics:
     * smaller σ → finer patterns, larger σ → coarser patterns. */
    if (b < OEP_EPSILON) return 0.0;
    double a_abs = fabs(a);
    return (2.0 * sqrt(2.0) / 3.0) * pow(a_abs, 1.5) * sqrt(c) / b;
}

double oep_lg_domain_wall_width(double a, double c) {
    /* Domain wall width:
     *   δ = ξ = √(c/|a|)
     *
     * This is the characteristic length scale over which the
     * order parameter changes from one phase to another.
     * It diverges as a → 0 (at the critical point). */
    double a_abs = fabs(a);
    if (a_abs < OEP_EPSILON) return 1e10;
    return sqrt(c / a_abs);
}

/* ── Correlation Function ──────────────────────────────────── */

double oep_lg_correlation_length(double a, double c) {
    /* Correlation length in the Gaussian (b = 0) theory:
     *   ξ = √(c/|a|)
     *
     * This is the length scale over which fluctuations decay.
     * At the critical point (a → 0): ξ → ∞ (critical opalescence). */
    double a_abs = fabs(a);
    if (a_abs < OEP_EPSILON) return 1e10;
    return sqrt(c / a_abs);
}

double oep_lg_ornstein_zernike_ck(double k, double a, double c,
                                   double temperature) {
    /* Ornstein-Zernike form of the static structure factor:
     *   C(k) = k_B·T / (a + c·k²)
     *
     * This is the Fourier transform of the Ornstein-Zernike
     * correlation function: C(r) ∝ (1/r) exp(-r/ξ) in 3D.
     *
     * At criticality (a = 0): C(k) ∝ 1/k² (divergence at k = 0).
     * This power-law behavior is the signature of scale invariance. */
    double denom = a + c * k * k;
    if (denom < OEP_EPSILON) return 1e10;
    return temperature / denom;
}

/* ── Superconductivity ─────────────────────────────────────── */

double oep_lg_superconductor_gap(double T, double Tc, double alpha0, double beta) {
    /* Superconducting order parameter magnitude near T_c:
     *   |Ψ|² = 0                     for T ≥ T_c
     *   |Ψ|² = (α₀/β)·(T_c - T)/T_c for T < T_c
     *
     * This is the mean-field Bardeen-Cooper-Schrieffer (BCS) result
     * expressed in Landau-Ginzburg form.
     *
     * The gap Δ(T) ∝ |Ψ| determines the superconducting transition. */
    if (T >= Tc || beta < OEP_EPSILON) return 0.0;
    return sqrt(alpha0 * (Tc - T) / (beta * Tc));
}

double oep_lg_coherence_length(double T, double Tc, double xi0) {
    /* Coherence length:
     *   ξ(T) = ξ₀·|T_c/(T_c - T)|^{1/2}
     *
     * ξ is the characteristic size of a Cooper pair
     * and the shortest length over which the order parameter
     * can vary appreciably. It diverges at T_c. */
    double dist = Tc - T;
    if (fabs(dist) < OEP_EPSILON) return 1e10;
    return xi0 * sqrt(fabs(Tc / dist));
}

double oep_lg_penetration_depth(double T, double Tc, double lambda0) {
    /* London penetration depth:
     *   λ_L(T) = λ_L(0)·|T_c/(T_c - T)|^{1/2}
     *
     * λ_L is the length scale over which an external magnetic field
     * decays inside a superconductor (Meissner effect). */
    double dist = Tc - T;
    if (fabs(dist) < OEP_EPSILON) return 1e10;
    return lambda0 * sqrt(fabs(Tc / dist));
}

double oep_lg_ginzburg_parameter(double lambda_L, double xi) {
    /* Ginzburg-Landau parameter:
     *   κ = λ_L / ξ
     *
     * This dimensionless ratio determines the superconductor type:
     *   κ < 1/√2 ≈ 0.71: Type I (complete flux expulsion)
     *   κ > 1/√2: Type II (vortex lattice in mixed state) */
    if (xi < OEP_EPSILON) return 1e10;
    return lambda_L / xi;
}

bool oep_lg_is_type_II(double kappa) {
    return kappa > (1.0 / sqrt(2.0));  /* ≈ 0.7071 */
}

double oep_lg_hc1(double H_c, double kappa) {
    /* Lower critical field:
     *   H_{c1} = H_c·(ln κ + 0.5) / (√2·κ)
     *
     * Below H_{c1}: perfect Meissner state (complete flux expulsion).
     * Above H_{c1}: flux begins to penetrate as quantized vortices. */
    if (kappa < OEP_EPSILON) return 0.0;
    return H_c * (log(kappa) + 0.5) / (sqrt(2.0) * kappa);
}

double oep_lg_hc2(double H_c, double kappa) {
    /* Upper critical field:
     *   H_{c2} = √2·κ·H_c
     *
     * Above H_{c2}: superconductivity is destroyed.
     * Between H_{c1} and H_{c2}: mixed state (Abrikosov vortex lattice). */
    return sqrt(2.0) * kappa * H_c;
}

/* ── Structural Phase Transitions ──────────────────────────── */

void oep_lg_strain_equilibria(double a, double b, double c,
                               double* eq_values, int* n_equilibria) {
    /* Find equilibrium strains for Landau-Devonshire free energy:
     *   F(e) = a·e²/2 + b·e³/3 + c·e⁴/4
     * → dF/de = a·e + b·e² + c·e³ = 0
     *
     * Solutions of ce³ + be² + ae = 0:
     *   e = 0 (always), and e = [-b ± √(b² - 4ac)]/(2c) (if discriminant ≥ 0)
     *
     * This describes martensitic transformations with:
     *   b = 0: symmetric (2nd order) transition at a = 0
     *   b ≠ 0: asymmetric (1st order) transition */
    if (!eq_values || !n_equilibria) return;
    eq_values[0] = 0.0;
    int n = 1;

    if (fabs(c) > OEP_EPSILON) {
        double disc = b * b - 4.0 * a * c;
        if (disc >= 0.0) {
            double sqrt_disc = sqrt(disc);
            eq_values[1] = (-b + sqrt_disc) / (2.0 * c);
            eq_values[2] = (-b - sqrt_disc) / (2.0 * c);
            n = 3;
        }
    } else if (fabs(b) > OEP_EPSILON) {
        eq_values[1] = -a / b;
        n = 2;
    }
    *n_equilibria = n;
}

double oep_lg_elastic_coupling_shift(double g, double e, double b) {
    /* Strain coupling to order parameter:
     *   F = aξ²/2 + bξ⁴/4 + g·e·ξ²
     *
     * Eliminating strain (minimize w.r.t e with elastic energy c·e²/2):
     *   ∂F/∂e = c·e + g·ξ² = 0 → e = -(g/c)·ξ²
     *
     * Effective free energy:
     *   F_eff(ξ) = aξ²/2 + (b - 2g²/c)·ξ⁴/4
     *
     * The coupling renormalizes the quartic coefficient:
     *   b_eff = b - 2g²/c
     *
     * If b_eff < 0: first-order transition (even if b > 0!).
     * This is the mechanism of the "strain-induced first-order transition." */
    return b - 2.0 * g * g * e;  /* e here represents 1/c (elastic compliance) */
}

double oep_lg_soft_mode_frequency(double T, double Tc, double omega0_sq) {
    /* Soft mode frequency near structural phase transition:
     *   ω²(T) = ω₀²·(T - T_c)/T_c   for T > T_c (soft mode)
     *   ω²(T) = ω₀²·2|T - T_c|/T_c   for T < T_c (hardening)
     *
     * The soft mode is the collective excitation whose frequency
     * goes to zero at the phase transition. It is the dynamical
     * signature of the incipient instability.
     *
     * Reference: Cochran (1960) — Soft mode theory of ferroelectrics */
    double t = (T - Tc) / Tc;
    if (t > 0.0)
        return omega0_sq * t;
    else
        return omega0_sq * 2.0 * fabs(t);
}
