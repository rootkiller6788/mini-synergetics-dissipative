#ifndef OEP_LANDAU_GINZBURG_H
#define OEP_LANDAU_GINZBURG_H
#include "oep_core.h"

/* ==============================================================
 * Landau-Ginzburg Theory of Phase Transitions
 * Based on: Landau & Lifshitz (1980), Ginzburg & Landau (1950)
 *
 * The Landau-Ginzburg formalism extends Landau's mean-field theory
 * to include spatial inhomogeneities through a gradient term:
 *
 *   F[ξ(r)] = ∫ dᵈr [ a(p)ξ²/2 + bξ⁴/4 + c/2 |∇ξ|² ]
 *
 * The dynamics follow the time-dependent Ginzburg-Landau (TDGL):
 *
 *   ∂ξ/∂t = -Γ δF/δξ + η(r,t)
 *         = -Γ [a(p)ξ + bξ³ - c∇²ξ] + η(r,t)
 *
 * where η(r,t) is spatiotemporal white noise (Model A dynamics).
 * ============================================================== */

/* ── Spatial Grid ──────────────────────────────────────── */

typedef struct {
    double* data;
    int     nx, ny;    /* Grid dimensions */
    double  dx, dy;    /* Grid spacing */
    double  Lx, Ly;    /* Physical domain size */
} OEPSpatialField;

typedef struct {
    double  a;          /* a(p) = (p - p_c) / p_c */
    double  b;          /* Quartic coefficient (> 0 for stability) */
    double  c;          /* Gradient stiffness coefficient */
    double  Gamma;      /* Kinetic coefficient */
    double  noise_amplitude; /* η intensity */
    int     nx, ny;
    double  dx, dy;
    double  Lx, Ly;
} OEPLandauGinzburg;

/* ── Spatial Field API ─────────────────────────────────── */

OEPSpatialField* oep_field_create(int nx, int ny, double Lx, double Ly);
void             oep_field_free(OEPSpatialField* field);
void             oep_field_set(OEPSpatialField* field, int i, int j, double val);
double           oep_field_get(const OEPSpatialField* field, int i, int j);
void             oep_field_zero(OEPSpatialField* field);
void             oep_field_set_random(OEPSpatialField* field, double amplitude);

/* Spatial derivatives (centered finite differences, periodic BC) */
double           oep_field_laplacian(const OEPSpatialField* field, int i, int j);
double           oep_field_grad_x(const OEPSpatialField* field, int i, int j);
double           oep_field_grad_y(const OEPSpatialField* field, int i, int j);
double           oep_field_grad_sq(const OEPSpatialField* field, int i, int j);

/* Compute spatial average <ξ> */
double           oep_field_spatial_mean(const OEPSpatialField* field);

/* Compute spatial variance <ξ²> - <ξ>² */
double           oep_field_spatial_variance(const OEPSpatialField* field);

/* Compute total free energy F[ξ] = ∫ [aξ²/2 + bξ⁴/4 + c/2|∇ξ|²] dᵈr */
double           oep_field_free_energy(const OEPSpatialField* field,
                                       double a, double b, double c,
                                       double dx, double dy);

/* ── Landau-Ginzburg Dynamics API ──────────────────────── */

OEPLandauGinzburg* oep_lg_create(int nx, int ny, double Lx, double Ly,
                                  double a, double b, double c, double Gamma);
void               oep_lg_free(OEPLandauGinzburg* lg);
void               oep_lg_set_noise(OEPLandauGinzburg* lg, double amplitude);

/* Time-dependent Ginzburg-Landau (TDGL) RHS at point (i,j):
 *   ∂ξ/∂t = -Γ [aξ + bξ³ - c∇²ξ] + η */
double             oep_lg_local_rhs(const OEPLandauGinzburg* lg,
                                     const OEPSpatialField* field,
                                     int i, int j);

/* One full time step of TDGL (explicit Euler):
 *   ξ_{n+1} = ξ_n + dt * ∂ξ/∂t */
void               oep_lg_step(OEPLandauGinzburg* lg, OEPSpatialField* field,
                                double dt);

/* Evolve for n_steps */
void               oep_lg_evolve(OEPLandauGinzburg* lg, OEPSpatialField* field,
                                  double dt, int n_steps,
                                  double* mean_trajectory,
                                  double* variance_trajectory);

/* ── Domain Wall Structure ─────────────────────────────── */

/* 1D domain wall profile:
 *   ξ(x) = ξ₀ tanh(x / (√2 ξ₀))  where ξ₀ = √(-a/b)  (a < 0) */
void               oep_lg_domain_wall_1d(OEPSpatialField* field,
                                          double x_center, double a, double b,
                                          double c, int nx);

/* Domain wall energy (surface tension):
 *   σ = (2√2/3) |a|^{3/2} √c / b */
double             oep_lg_domain_wall_energy(double a, double b, double c);

/* Domain wall width:
 *   δ = ξ₀  where ξ₀ = √(c/|a|)  (the correlation length) */
double             oep_lg_domain_wall_width(double a, double c);

/* ── Correlation Function & Structure Factor ───────────── */

/* Spatial correlation function:
 *   C(r) = <ξ(0)ξ(r)> - <ξ>²
 * For the Gaussian (b=0) theory:
 *   C(r) ∝ r^{-(d-2)} exp(-r/ξ)  with ξ = √(c/|a|) */
double             oep_lg_correlation_length(double a, double c);

/* Ornstein-Zernike form of the correlation function:
 *   C(k) = k_B T / (a + c k²)
 * in Fourier space. Compute C(k) for a given wavevector. */
double             oep_lg_ornstein_zernike_ck(double k, double a, double c,
                                               double temperature);

/* ── Superconductivity Application ─────────────────────── */

/* Ginzburg-Landau equations for superconductivity:
 *   αΨ + β|Ψ|²Ψ + (1/(2m*))(ħ∇/i - e*A)²Ψ = 0
 *   J = (eħ/m*) |Ψ|² (∇φ - (e/ħ)A)
 * where Ψ = order parameter (superconducting gap), A = vector potential */

/* Order parameter magnitude near T_c in zero field:
 *   |Ψ|² = -α/β = (α₀/β) (T_c - T)/T_c */
double             oep_lg_superconductor_gap(double T, double Tc, double alpha0, double beta);

/* Coherence length: ξ(T) = ξ₀ |T_c/(T_c - T)|^{1/2} */
double             oep_lg_coherence_length(double T, double Tc, double xi0);

/* London penetration depth: λ_L(T) = λ_L(0) |T_c/(T_c - T)|^{1/2} */
double             oep_lg_penetration_depth(double T, double Tc, double lambda0);

/* Ginzburg-Landau parameter: κ = λ_L/ξ */
double             oep_lg_ginzburg_parameter(double lambda_L, double xi);

/* Type-I vs Type-II superconductor:
 *   κ < 1/√2: Type I
 *   κ > 1/√2: Type II */
bool               oep_lg_is_type_II(double kappa);

/* Lower critical field H_{c1} */
double             oep_lg_hc1(double H_c, double kappa);

/* Upper critical field H_{c2} = √2 κ H_c */
double             oep_lg_hc2(double H_c, double kappa);

/* ── Structural Phase Transitions ──────────────────────── */

/* Strain order parameter (e.g. for martensitic transformation):
 *   F(e) = a e²/2 + b e³/3 + c e⁴/4
 * where e is the strain. This is the Landau-Devonshire form. */
void               oep_lg_strain_equilibria(double a, double b, double c,
                                             double* eq_values, int* n_equilibria);

/* Elastic free energy coupling to electromagnetic order parameter:
 *   F = aξ²/2 + bξ⁴/4 + g e ξ² */
double             oep_lg_elastic_coupling_shift(double g, double e, double b);

/* Soft mode frequency near structural transition:
 *   ω²(T) = ω₀² (T - T_c)/T_c  (soft mode) */
double             oep_lg_soft_mode_frequency(double T, double Tc, double omega0_sq);

#endif
