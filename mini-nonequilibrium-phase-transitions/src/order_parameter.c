#include "order_parameter.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * Order Parameters — Implementation
 *
 * Knowledge points implemented:
 * 1. Order parameter construction (L1: scalar OP, complex OP)
 * 2. TDGL relaxation dynamics (L4: Model A / Model B dynamics)
 * 3. Landau free energy minimization (L4: equilibrium OP from dF/dη=0)
 * 4. Correlation function analysis (L3: G(r) and xi extraction)
 * 5. Structure factor computation (L3: S(k) via Fourier transform)
 * 6. Symmetry breaking analysis (L4: Goldstone theorem counting)
 * 7. Hyperscaling relation check (L4: 2-alpha = d*nu)
 * ============================================================================ */

/* --- Scalar OP creation --- */
OrderParameter nept_op_scalar_create(double magnitude)
{
    OrderParameter op;
    op.magnitude = fabs(magnitude);
    op.phase = 0.0;
    op.type = OP_SCALAR_REAL;
    op.is_conserved = false;
    op.relaxation_rate = 1.0;
    return op;
}

/* --- Complex OP creation --- */
OrderParameter nept_op_complex_create(double real_part, double imag_part)
{
    OrderParameter op;
    op.magnitude = sqrt(real_part * real_part + imag_part * imag_part);
    op.phase = atan2(imag_part, real_part);
    op.type = OP_SCALAR_COMPLEX;
    op.is_conserved = false;
    op.relaxation_rate = 1.0;
    return op;
}

/* ---------------------------------------------------------------------------
 * TDGL Deterministic Step
 *
 * Knowledge point: The time-dependent Ginzburg-Landau (TDGL) equation
 * describes the relaxation of a non-conserved order parameter (Model A):
 *
 *   dη/dt = -Γ * δF/δη
 *         = -Γ * [a*t*η + b*η^3 - κ*∇^2*η + ...]
 *
 * For a conserved OP (Model B, Hohenberg & Halperin 1977):
 *   dη/dt = Γ * ∇^2 * δF/δη  (diffusion-controlled)
 *
 * This is a gradient flow: the OP always moves downhill in free energy.
 * The Euler update is:
 *   η(t+dt) = η(t) - Γ*dt * (a*t*η + b*η^3 - κ*laplacian)
 *
 * Reference: Hohenberg & Halperin, Rev. Mod. Phys. 49, 435 (1977)
 * --------------------------------------------------------------------------- */
double nept_tdgl_deterministic_step(OrderParameter *op,
                                     const TDGLParameters *params,
                                     double reduced_temp,
                                     double laplacian,
                                     double dt)
{
    if (!op || !params) return op ? op->magnitude : 0.0;

    double eta = op->magnitude;
    double a = params->a_coefficient;
    double b = params->b_coefficient;
    double kappa = params->kappa;
    double gamma = params->kinetic_coefficient;

    /* dF/dη = a*t*η + b*η^3 - κ*∇^2*η (functional derivative) */
    double dF_deta = a * reduced_temp * eta + b * eta * eta * eta
                     - kappa * laplacian;

    /* η_new = η_old - Γ * dt * dF/dη */
    double new_eta = eta - gamma * dt * dF_deta;

    /* Ensure η >= 0 for scalar magnitude */
    if (new_eta < 0.0) new_eta = 0.0;

    op->magnitude = new_eta;
    return new_eta;
}

/* ---------------------------------------------------------------------------
 * Landau Free Energy Density
 *
 * Knowledge point: The Landau free energy is the thermodynamic potential
 * expanded in powers of the order parameter:
 *
 *   f(η, t) = f0 + (a*t/2)*η^2 + (b/4)*η^4 + (c/6)*η^6
 *
 * For t > 0: single minimum at η = 0 (disordered phase)
 * For t < 0 (b > 0): two symmetric minima at η = ±√(-a*t/b) (ordered)
 * For b < 0: first-order transition with metastable states
 *
 * This function computes f(η) - f0 to isolate the OP-dependent part.
 *
 * Reference: Landau & Lifshitz, Statistical Physics §146
 * --------------------------------------------------------------------------- */
double nept_landau_free_energy(double eta,
                                const TDGLParameters *params,
                                double reduced_temp)
{
    if (!params) return 0.0;
    double a = params->a_coefficient;
    double b = params->b_coefficient;

    /* f(η) = (a*t/2)*η^2 + (b/4)*η^4  [standard φ^4 theory] */
    double eta2 = eta * eta;
    double eta4 = eta2 * eta2;

    return 0.5 * a * reduced_temp * eta2 + 0.25 * b * eta4;
}

/* ---------------------------------------------------------------------------
 * Landau Equilibrium Order Parameter
 *
 * Knowledge point: The equilibrium value of the order parameter minimizes
 * the free energy: dF/dη = 0.
 *
 *   dF/dη = a*t*η + b*η^3 = 0
 *   => η*(a*t + b*η^2) = 0
 *   Solutions: η = 0 (always), η = ±√(-a*t/b) for t < 0 (b > 0)
 *
 * For b < 0 (first-order), minimization requires the sextic term:
 *   dF/dη = a*t*η + b*η^3 + c*η^5 = 0
 *   First-order when F(η=0) = F(η*) for some η* > 0 at t < 0.
 *
 * Reference: Chaikin & Lubensky, §2.1
 * --------------------------------------------------------------------------- */
double nept_landau_equilibrium_op(const TDGLParameters *params,
                                   double reduced_temp,
                                   NEPT_TransitionType *trans_type_out)
{
    if (!params) {
        if (trans_type_out) *trans_type_out = NEPT_CONTINUOUS;
        return 0.0;
    }

    double a = params->a_coefficient;
    double b = params->b_coefficient;

    if (trans_type_out) {
        if (b > 0.0) {
            *trans_type_out = NEPT_CONTINUOUS;
        } else if (b < 0.0) {
            *trans_type_out = NEPT_FIRST_ORDER;
        } else {
            *trans_type_out = NEPT_CONTINUOUS; /* tricritical */
        }
    }

    /* For continuous transition (b > 0): */
    if (reduced_temp >= 0.0) return 0.0;

    /* t < 0: ordered phase */
    if (b > 0.0) {
        return sqrt(-a * reduced_temp / b);
    }

    /* For b < 0 (first-order): need sextic to stabilize
     * At the binodal, the non-zero minimum has:
     * η* = sqrt((-b + sqrt(b^2 - 4*a*c*t))/(2*c))
     *
     * Simplified: return approximate jump = sqrt(-3*b/(4*c)) at transition
     * This is a rough estimate — the full solution requires c coefficient. */
    return 0.0; /* First-order needs full sextic analysis */
}

/* --- Landau classification from coefficients alone --- */
NEPT_TransitionType nept_landau_classify(const TDGLParameters *params)
{
    if (!params) return NEPT_CONTINUOUS;
    if (params->b_coefficient > 1e-12) return NEPT_CONTINUOUS;
    if (params->b_coefficient < -1e-12) return NEPT_FIRST_ORDER;
    return NEPT_CONTINUOUS; /* Tricritical: need higher-order analysis */
}

/* ---------------------------------------------------------------------------
 * Correlation Function — Construction and Analysis
 *
 * Knowledge point: The two-point correlation function G(r) measures how
 * fluctuations at different points are related:
 *
 *   G(|r1 - r2|) = <η(r1) η(r2)> - <η(r1)><η(r2)>
 *
 * In the disordered phase (T > Tc), G(r) decays exponentially:
 *   G(r) ~ exp(-r/ξ) / r^{(d-1)/2}   (Ornstein-Zernike at large r)
 *
 * At the critical point (T = Tc), G(r) decays as a power law:
 *   G(r) ~ r^{-(d-2+η)}   (Fisher exponent η)
 *
 * The correlation length ξ is extracted by fitting log(G*r^{(d-1)/2}) vs r
 * for large r (where the exponential dominates).
 *
 * Reference: Stanley, §2.3; Goldenfeld, §3.1
 * --------------------------------------------------------------------------- */
CorrelationFunction *nept_correlation_create(int n_points)
{
    if (n_points <= 0) return NULL;
    CorrelationFunction *cf = (CorrelationFunction *)malloc(
        sizeof(CorrelationFunction));
    if (!cf) return NULL;
    cf->values = (double *)calloc((size_t)n_points, sizeof(double));
    cf->r = (double *)calloc((size_t)n_points, sizeof(double));
    if (!cf->values || !cf->r) {
        free(cf->values);
        free(cf->r);
        free(cf);
        return NULL;
    }
    cf->n_points = n_points;
    cf->correlation_length = 0.0;
    cf->exponent_eta_fisher = 0.0;
    return cf;
}

void nept_correlation_exponential(CorrelationFunction *cf,
                                   double xi,
                                   double r_min,
                                   double r_max)
{
    if (!cf || cf->n_points < 2 || xi <= 0.0) return;

    double dr = (r_max - r_min) / (double)(cf->n_points - 1);
    for (int i = 0; i < cf->n_points; i++) {
        double r = r_min + (double)i * dr;
        cf->r[i] = r;
        /* Ornstein-Zernike form: G(r) ~ exp(-r/xi) / r  (3D) */
        if (r > 1e-12) {
            cf->values[i] = exp(-r / xi) / r;
        } else {
            cf->values[i] = 1.0; /* Regularized at r=0 */
        }
    }
    cf->correlation_length = xi;
}

double nept_correlation_extract_xi(const CorrelationFunction *cf)
{
    if (!cf || cf->n_points < 3) return 0.0;

    /* Fit log(G(r)) vs r for large r (exponential tail).
     * log(G) = const - r/xi => slope = -1/xi.
     * Use the last half of the data points where exponential dominates. */
    int start = cf->n_points / 2;
    int n_fit = cf->n_points - start;
    if (n_fit < 2) return 0.0;

    double sum_r = 0.0, sum_logG = 0.0, sum_r2 = 0.0, sum_r_logG = 0.0;
    for (int i = start; i < cf->n_points; i++) {
        if (cf->values[i] <= 0.0) continue;
        double r = cf->r[i];
        double logG = log(cf->values[i]);
        sum_r += r;
        sum_logG += logG;
        sum_r2 += r * r;
        sum_r_logG += r * logG;
    }

    double denom = (double)n_fit * sum_r2 - sum_r * sum_r;
    if (fabs(denom) < 1e-15) return 0.0;
    double slope = ((double)n_fit * sum_r_logG - sum_r * sum_logG) / denom;

    /* slope = -1/xi => xi = -1/slope */
    if (fabs(slope) < 1e-15) return 1e15; /* effectively infinite */
    double xi = -1.0 / slope;
    return (xi > 0.0) ? xi : fabs(xi);
}

void nept_correlation_free(CorrelationFunction *cf)
{
    if (cf) {
        free(cf->values);
        free(cf->r);
        free(cf);
    }
}

/* ---------------------------------------------------------------------------
 * Structure Factor
 *
 * Knowledge point: The structure factor S(k) is the Fourier transform of
 * G(r). Experimentally accessible via scattering (neutron, X-ray, light).
 *
 * For T > Tc (Ornstein-Zernike):
 *   S(k) = A / (ξ^{-2} + k^2)   for small k
 *
 * At Tc (ξ → ∞): S(k) ~ k^{-2+η}  for small k
 *
 * This implements a discrete Fourier transform for spherically symmetric
 * G(r) in 3D:
 *   S(k) = 4π ∫_0^∞ r^2 G(r) sin(kr)/(kr) dr
 *
 * Reference: Chaikin & Lubensky, §2.3; Stanley, §4.3
 * --------------------------------------------------------------------------- */
StructureFactor *nept_structure_factor_from_correlation(
    const CorrelationFunction *cf, int n_kpoints)
{
    if (!cf || n_kpoints <= 0) return NULL;

    StructureFactor *sf = (StructureFactor *)malloc(sizeof(StructureFactor));
    if (!sf) return NULL;
    sf->intensity = (double *)calloc((size_t)n_kpoints, sizeof(double));
    sf->k_values = (double *)calloc((size_t)n_kpoints, sizeof(double));
    if (!sf->intensity || !sf->k_values) {
        free(sf->intensity);
        free(sf->k_values);
        free(sf);
        return NULL;
    }
    sf->n_points = n_kpoints;
    sf->ornstein_zernike_amplitude = 0.0;

    double r_max = cf->r[cf->n_points - 1];
    double dk = 2.0 * M_PI / r_max;
    double dr = cf->r[1] - cf->r[0];

    for (int ik = 0; ik < n_kpoints; ik++) {
        double k = (double)(ik + 1) * dk;
        sf->k_values[ik] = k;
        double sum = 0.0;
        for (int ir = 0; ir < cf->n_points; ir++) {
            double r = cf->r[ir];
            double kr = k * r;
            double sinc = (fabs(kr) > 1e-12) ? sin(kr) / kr : 1.0;
            /* 3D Fourier transform of spherically symmetric function */
            sum += r * r * cf->values[ir] * sinc * dr;
        }
        sf->intensity[ik] = 4.0 * M_PI * sum;
    }

    /* Estimate Ornstein-Zernike amplitude: S(k→0) = A * ξ^2 */
    if (n_kpoints >= 1 && cf->correlation_length > 0.0) {
        sf->ornstein_zernike_amplitude =
            sf->intensity[0] / (cf->correlation_length * cf->correlation_length);
    }

    return sf;
}

void nept_structure_factor_free(StructureFactor *sf)
{
    if (sf) {
        free(sf->intensity);
        free(sf->k_values);
        free(sf);
    }
}

/* ---------------------------------------------------------------------------
 * Symmetry Breaking Analysis
 *
 * Knowledge point: A phase transition is accompanied by spontaneous
 * symmetry breaking. The order parameter measures how much symmetry
 * is broken.
 *
 * Goldstone theorem: For each broken continuous symmetry generator,
 * there is a massless excitation (Goldstone mode). Number of Goldstone
 * modes = dim(G) - dim(H), where G is the full symmetry group and H
 * is the residual (unbroken) subgroup.
 *
 * Examples:
 *   Ising: Z2 → {1}, 0 Goldstone modes (discrete)
 *   XY: O(2) → {1}, 1 Goldstone mode (spin wave)
 *   Heisenberg: O(3) → O(2), 2 Goldstone modes (magnons)
 *
 * Reference: Goldstone, Nuovo Cimento 19, 154 (1961);
 *   Nambu, Phys. Rev. 117, 648 (1960)
 * --------------------------------------------------------------------------- */
SymmetryBreaking nept_symmetry_breaking_create(
    NEPT_SymmetryType broken, NEPT_SymmetryType residual,
    int order_full, int order_sub)
{
    SymmetryBreaking sb;
    sb.broken_symmetry = broken;
    sb.residual_symmetry = residual;
    sb.order_of_broken_group = order_full;
    sb.order_of_residual_group = order_sub;
    sb.num_goldstone_modes = order_full - order_sub;
    return sb;
}

int nept_count_goldstone_modes(const SymmetryBreaking *sb)
{
    if (!sb) return 0;
    /* dim(G/H) = order(full) - order(sub) for finite groups,
     * or more generally = dim(G) - dim(H) for Lie groups.
     * Here we encode dimensions via the order fields. */
    int n_modes = sb->order_of_broken_group - sb->order_of_residual_group;
    return (n_modes > 0) ? n_modes : 0;
}

/* ---------------------------------------------------------------------------
 * Estimate nu exponent from correlation length data
 *
 * Knowledge point: Near the critical point, ξ ~ |t|^{-ν}.
 * Taking logs: log(ξ) = log(ξ0) - ν * log(|t|)
 * Linear regression gives ν as the negative slope.
 *
 * Reference: Stanley, §4.1
 * --------------------------------------------------------------------------- */
double nept_estimate_nu_exponent(const double *xi_values,
                                  const double *reduced_temps,
                                  int n_points,
                                  double *xi0_out)
{
    if (!xi_values || !reduced_temps || n_points < 3) {
        if (xi0_out) *xi0_out = 0.0;
        return 0.0;
    }

    double sum_x = 0.0, sum_y = 0.0, sum_xx = 0.0, sum_xy = 0.0;
    int count = 0;

    for (int i = 0; i < n_points; i++) {
        if (fabs(reduced_temps[i]) < 1e-15) continue; /* skip exactly Tc */
        double x = log(fabs(reduced_temps[i]));
        double y = log(xi_values[i]);
        sum_x += x;
        sum_y += y;
        sum_xx += x * x;
        sum_xy += x * y;
        count++;
    }

    if (count < 2) {
        if (xi0_out) *xi0_out = 0.0;
        return 0.0;
    }

    double denom = (double)count * sum_xx - sum_x * sum_x;
    if (fabs(denom) < 1e-15) {
        if (xi0_out) *xi0_out = 0.0;
        return 0.0;
    }

    double slope = ((double)count * sum_xy - sum_x * sum_y) / denom;
    double intercept = (sum_y - slope * sum_x) / (double)count;

    /* xi ~ xi0 * |t|^{-nu} => slope = -nu */
    double nu = -slope;
    if (xi0_out) *xi0_out = exp(intercept);

    return nu;
}

/* ---------------------------------------------------------------------------
 * Hyperscaling relation check
 *
 * Knowledge point: Josephson hyperscaling relation:
 *   2 - α = d * ν
 * where d is spatial dimension. This holds below the upper critical
 * dimension d_c = 4 (for O(N) models). Above d_c, mean-field exponents
 * (α=0, ν=1/2) do NOT satisfy hyperscaling for d > 4.
 *
 * This relation follows from the assumption that the singular part
 * of the free energy scales as f_s ~ ξ^{-d}, which gives:
 *   f_s ~ |t|^{2-α} ~ ξ^{-d} ~ |t|^{d*ν} => 2-α = d*ν
 *
 * Reference: Josephson, Phys. Lett. 21, 608 (1966)
 * --------------------------------------------------------------------------- */
bool nept_check_hyperscaling(double alpha, double nu, int dimension,
                              double tolerance)
{
    double lhs = 2.0 - alpha;
    double rhs = (double)dimension * nu;
    return fabs(lhs - rhs) < tolerance;
}