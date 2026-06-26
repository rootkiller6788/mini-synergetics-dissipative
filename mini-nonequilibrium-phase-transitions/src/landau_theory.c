#include "landau_theory.h"
#include <stdlib.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * Landau-Ginzburg Theory — Implementation
 *
 * Knowledge points:
 * 1. Equilibrium OP from Landau free energy minimization (L4)
 * 2. Susceptibility in mean-field theory (L4: Curie-Weiss law)
 * 3. Specific heat jump at Tc (L4: Landau specific heat discontinuity)
 * 4. Binodal/spinodal for first-order transitions (L4: coexistence curves)
 * 5. Ginzburg criterion for mean-field validity (L8)
 * 6. Nucleation barrier for first-order transitions (L5: CNT)
 * 7. O(N) model equilibrium (L3: multi-component OP)
 * 8. Tricritical point analysis (L8)
 * 9. Hysteresis loop area (L5: energy dissipation proxy)
 * ============================================================================ */

/* ---------------------------------------------------------------------------
 * Landau Equilibrium — Full scalar phi^6 theory
 *
 * Free energy: F(η) = f0 + (a*t/2)*η^2 + (b/4)*η^4 + (c/6)*η^6
 * Stationarity: dF/dη = η*(a*t + b*η^2 + c*η^4) = 0
 *
 * Solutions:
 *   η = 0 (disordered)
 *   η^2 = [-b ± sqrt(b^2 - 4*a*c*t)] / (2*c)  for the ordered branch
 *
 * Selection rule: choose the solution with lowest F (global minimum).
 *
 * For b > 0 (continuous): at t < 0, non-zero η ~ |t|^{1/2}
 * For b < 0 (first-order): discontinuous jump at t_binodal > 0
 * --------------------------------------------------------------------------- */

double nept_landau_equilibrium(LandauCoefficients lc,
                                double reduced_temp,
                                MeanFieldResults *results)
{
    if (!results) return 0.0;

    double a = lc.a;
    double b = lc.b;
    double c = lc.c;
    double t = reduced_temp;

    /* Coefficient of η^2: A = a*t, coefficient of η^4: B = b, coefficient of η^6: C = c */
    double A = a * t;

    /* Compute free energies and find the global minimum */
    double eta_eq = 0.0;
    double F0 = lc.f0;                       /* F at η = 0 */
    double F_nonzero = F0;

    if (c > 1e-15) {
        /* Solve quadratic in η^2: c*η^4 + b*η^2 + a*t = 0
         * η^2 = [-b ± sqrt(b^2 - 4*a*c*t)] / (2*c) */
        double discriminant = b * b - 4.0 * a * c * t;
        if (discriminant >= 0.0) {
            double sqrt_d = sqrt(discriminant);
            double eta2_plus  = (-b + sqrt_d) / (2.0 * c);
            double eta2_minus = (-b - sqrt_d) / (2.0 * c);

            /* Evaluate F at each positive solution */
            if (eta2_plus > 1e-15) {
                double eta = sqrt(eta2_plus);
                double Fp = F0 + 0.5 * A * eta2_plus
                            + 0.25 * b * eta2_plus * eta2_plus
                            + (c / 6.0) * eta2_plus * eta2_plus * eta2_plus;
                if (Fp < F0) { eta_eq = eta; F_nonzero = Fp; }
            }
            if (eta2_minus > 1e-15 && eta2_minus != eta2_plus) {
                double eta = sqrt(eta2_minus);
                double Fm = F0 + 0.5 * A * eta2_minus
                            + 0.25 * b * eta2_minus * eta2_minus
                            + (c / 6.0) * eta2_minus * eta2_minus * eta2_minus;
                if (Fm < F_nonzero) { eta_eq = eta; F_nonzero = Fm; }
            }
        }
    } else {
        /* Standard phi^4: only η^2 and η^4 terms (c = 0) */
        if (b > 1e-15 && A < 0.0) {
            eta_eq = sqrt(-A / b);
            double eta2 = eta_eq * eta_eq;
            F_nonzero = F0 + 0.5 * A * eta2 + 0.25 * b * eta2 * eta2;
        }
    }

    /* Fill results */
    results->spontaneous_op = eta_eq;
    results->free_energy = (eta_eq > 1e-12) ? F_nonzero : F0;
    results->free_energy_barrier = 0.0;

    /* Susceptibility: χ = 1/(d^2F/dη^2)|_{η_eq}
     * d^2F/dη^2 = a*t + 3*b*η^2 + 5*c*η^4 */
    double d2F = A + 3.0 * b * eta_eq * eta_eq
                 + 5.0 * c * eta_eq * eta_eq * eta_eq * eta_eq;
    results->susceptibility = (fabs(d2F) > 1e-15) ? 1.0 / d2F : 1e15;

    /* Specific heat: C = -T * d^2F/dT^2
     * d^2F/dT^2 = (a/(Tc^2))*η_eq^2 + ... (simplified) */
    results->specific_heat = 0.0; /* needs Tc for full computation */
    results->specific_heat_jump = nept_landau_specific_heat_jump(lc, 1.0);

    /* Free energy barrier for first-order: difference between maxima and minima */
    if (b < -1e-12 && c > 1e-12) {
        /* Barrier = max(F) - min(F) between the two minima */
        double dF2_d2eta_at_0 = A; /* curvature at origin */
        if (dF2_d2eta_at_0 < 0.0 && eta_eq > 1e-12) {
            /* Barrier = F(0) - F(η_eq) if F(0) is a local max,
             * or F(local_max) - F(η_eq) more generally */
            results->free_energy_barrier = F0 - F_nonzero;
            if (results->free_energy_barrier < 0.0)
                results->free_energy_barrier = 0.0;
        }
    }

    return eta_eq;
}

/* --- Transition type from coefficients --- */
NEPT_TransitionType nept_landau_transition_type(LandauCoefficients lc)
{
    if (lc.b > 1e-12) return NEPT_CONTINUOUS;
    if (lc.b < -1e-12) return NEPT_FIRST_ORDER;
    /* b == 0: tricritical point (need higher-order analysis) */
    if (lc.c > 1e-12) return NEPT_CONTINUOUS; /* tricritical: still continuous at this point */
    return NEPT_CONTINUOUS;
}

/* ---------------------------------------------------------------------------
 * Mean-Field Susceptibility — Curie-Weiss Law
 *
 * Knowledge point: In Landau theory, the susceptibility diverges as
 * χ ~ |t|^{-γ} with γ = 1 (mean-field).
 *
 * For T > Tc (η = 0):
 *   χ = 1/(a*t)  =>  χ ~ t^{-1}  =>  Curie-Weiss law
 *
 * For T < Tc (η = √(-a*t/b)):
 *   χ = 1/(4*a*|t|)  =>  amplitude ratio χ_+/χ_- = b/(2*a)?
 *   More precisely: d^2F/dη^2 = a*t + 3*b*η_eq^2
 *                   = a*t + 3*b*(-a*t/b) = -2*a*t = 2*a*|t|
 *   So χ = 1/(2*a*|t|) for t < 0.
 *
 * Reference: Stanley, §5.2
 * --------------------------------------------------------------------------- */
double nept_landau_susceptibility(LandauCoefficients lc,
                                   double reduced_temp)
{
    double t = reduced_temp;
    double a = lc.a;
    double b = lc.b;

    if (t > 0.0) {
        /* Disordered phase: η = 0, d^2F/dη^2 = a*t */
        double denom = a * t;
        return (fabs(denom) > 1e-15) ? 1.0 / denom : 1e15;
    } else if (t < 0.0 && b > 1e-12) {
        /* Ordered phase: η_eq = sqrt(-a*t/b)
         * d^2F/dη^2 = a*t + 3*b*(-a*t/b) = -2*a*t = 2*a*|t| */
        double denom = 2.0 * a * fabs(t);
        return (fabs(denom) > 1e-15) ? 1.0 / denom : 1e15;
    } else {
        /* Exactly at Tc or first-order: divergent */
        return 1e15;
    }
}

/* ---------------------------------------------------------------------------
 * Specific Heat Jump — Landau 1937
 *
 * Knowledge point: In mean-field theory, the specific heat has a finite
 * discontinuity ΔC at Tc:
 *
 *   C(T > Tc) = C0(T)  (background)
 *   C(T < Tc) = C0(T) + a^2 * T / (2*b)
 *
 * The jump:
 *   ΔC = C(Tc^-) - C(Tc^+) = a^2 * Tc / (2*b)
 *
 * For first-order (b < 0), latent heat replaces specific heat jump:
 *   L = Tc * ΔS = Tc * a * Δ(η^2)/2
 *
 * Reference: Landau & Lifshitz, Statistical Physics §143
 * --------------------------------------------------------------------------- */
double nept_landau_specific_heat_jump(LandauCoefficients lc, double Tc)
{
    double a = lc.a;
    double b = lc.b;

    if (fabs(b) < 1e-15) return 1e15; /* Tricritical: C diverges (alpha = 1/2 in mean-field) */
    if (b < 0.0) return 0.0;           /* First-order: latent heat, not C jump */

    /* dC = a^2 * Tc / (2*b) */
    return a * a * Tc / (2.0 * b);
}

/* ---------------------------------------------------------------------------
 * Binodal Temperature — First-Order Transitions
 *
 * Knowledge point: The binodal is the coexistence curve where two phases
 * have equal free energy: F(η=0) = F(η=η*).
 *
 * For phi^6 theory with b < 0, c > 0:
 *   At the binodal: F(η*) = F(0) => (a*t/2)*η*^2 + (b/4)*η*^4 + (c/6)*η*^6 = 0
 *   Combined with dF/dη(η*) = 0 => a*t + b*η*^2 + c*η*^4 = 0
 *   Solving: η*^2 = -3*b/(4*c), t_binodal = 3*b^2/(16*a*c)
 *
 * Spinodal points: d^2F/dη^2 = 0 where metastability ends.
 *
 * Reference: Chaikin & Lubensky, §2.2.3
 * --------------------------------------------------------------------------- */
double nept_landau_binodal_temperature(LandauCoefficients lc,
                                        FirstOrderCoexistence *coex)
{
    double a = lc.a;
    double b = lc.b;
    double c = lc.c;

    /* Only first-order (b < 0, c > 0) has binodal */
    if (b >= -1e-12 || c < 1e-12) {
        if (coex) {
            coex->binodal_t = 0.0;
            coex->spinodal_upper = 0.0;
            coex->spinodal_lower = 0.0;
            coex->latent_heat = 0.0;
            coex->coexistence_width = 0.0;
        }
        return 0.0;
    }

    /* At binodal: a*t_b + (3*b^2)/(16*c) = 0 => t_b = -(3*b^2)/(16*a*c) */
    double t_binodal = -(3.0 * b * b) / (16.0 * a * c);
    double eta_star = sqrt(-3.0 * b / (4.0 * c));

    if (coex) {
        coex->binodal_t = t_binodal;
        coex->coexistence_width = eta_star;

        /* Spinodal: d^2F/dη^2 = 0 => a*t + 3*b*η^2 + 5*c*η^4 = 0
         * At η=0: d^2F/dη^2 = a*t = 0 => upper spinodal at t = 0
         * At η=η*: solve for t */
        coex->spinodal_upper = 0.0; /* trivial branch loses stability */

        /* Lower spinodal: non-trivial branch stability limit
         * d^2F/dη^2(η*) = 0 => a*t_s + 3*b*η*^2 + 5*c*η*^4 = 0
         * t_s = t_binodal - b^2/(4*a*c)
         * = -(3*b^2)/(16*a*c) - (4*b^2)/(16*a*c) = -(7*b^2)/(16*a*c) */
        double b2 = b * b;
        double ac = a * c;
        if (fabs(ac) > 1e-15) {
            coex->spinodal_lower = -(7.0 * b2) / (16.0 * ac);
        } else {
            coex->spinodal_lower = t_binodal;
        }

        /* Latent heat: L = T * ΔS, ΔS = -(∂F/∂T)_η * Δη
         * In Landau: ∂F/∂T = (a/(2*Tc))*η^2  (from t = (T-Tc)/Tc)
         * ΔS = (a/(2*Tc))*η*^2, L = T * (a/(2*Tc))*η*^2 ≈ (a/2)*η*^2 */
        coex->latent_heat = 0.5 * a * eta_star * eta_star;
    }

    return t_binodal;
}

/* --- Spinodal points (standalone) --- */
void nept_landau_spinodal_points(LandauCoefficients lc,
                                  double *spinodal_upper,
                                  double *spinodal_lower)
{
    FirstOrderCoexistence coex;
    nept_landau_binodal_temperature(lc, &coex);
    if (spinodal_upper) *spinodal_upper = coex.spinodal_upper;
    if (spinodal_lower) *spinodal_lower = coex.spinodal_lower;
}

/* ---------------------------------------------------------------------------
 * Ginzburg Criterion
 *
 * Knowledge point: Levanyuk-Ginzburg criterion determines when mean-field
 * theory breaks down due to fluctuations. Fluctuations of the order
 * parameter in a correlated volume ξ^d must be small compared to η^2:
 *
 *   <(Δη)^2>_ξ << η^2
 *
 * This gives:
 *   t_G ≈ (k_B / (ΔC * ξ0^3))^2   (Ginzburg number)
 *
 * For 3D: t_G ~ 10^{-2} to 10^{-4} in superconductors, ~ 0.1 in liquid crystals.
 * For d ≥ 4 (upper critical dimension): t_G → 0, mean-field is exact.
 *
 * Reference: Ginzburg, Sov. Phys. Solid State 2, 1824 (1960);
 *   Als-Nielsen & Birgeneau, Am. J. Phys. 45, 554 (1977)
 * --------------------------------------------------------------------------- */
GinzburgCriterion nept_ginzburg_criterion(double a0, double xi0,
                                           double delta_C, double Tc,
                                           int spatial_dimension)
{
    GinzburgCriterion gc;
    (void)a0; (void)Tc;
    gc.upper_critical_dimension = 4;

    if (spatial_dimension >= 4) {
        gc.ginzburg_reduced_temp = 0.0;
        gc.ginzburg_length = 1e15;
        gc.mean_field_valid = true;
        return gc;
    }

    /* t_G = (k_B * Tc / (a0 * xi0^3 * Delta_C))^2
     * Here we use a simplified form: t_G = (1/(Delta_C * xi0^3))^2
     * The exact prefactor involves kB and a0 which we normalize out. */
    double xi03 = xi0 * xi0 * xi0;
    if (fabs(delta_C * xi03) < 1e-15) {
        gc.ginzburg_reduced_temp = 0.0;
        gc.mean_field_valid = true;
    } else {
        double gi = 1.0 / (fabs(delta_C) * xi03);
        gc.ginzburg_reduced_temp = gi * gi;
        gc.mean_field_valid = false;
    }

    /* ξ_G = ξ0 * |t_G|^{-ν} with ν = 1/2 (mean-field) */
    if (gc.ginzburg_reduced_temp > 1e-15) {
        gc.ginzburg_length = xi0 / sqrt(gc.ginzburg_reduced_temp);
    } else {
        gc.ginzburg_length = 1e15;
    }

    return gc;
}

/* ---------------------------------------------------------------------------
 * Nucleation Barrier
 *
 * Knowledge point: In a first-order transition, the system can remain
 * in a metastable state (supercooled liquid, superheated solid).
 * Transition occurs via nucleation of the stable phase.
 *
 * Classical Nucleation Theory (CNT):
 *   ΔG(r) = -(4π/3)*r^3*Δg_v + 4π*r^2*σ
 *   Barrier: ΔG* = 16π*σ^3 / (3*Δg_v^2)
 *
 * In Landau theory, the nucleation barrier is the free energy of a
 * critical droplet, approximately equal to the barrier height between
 * the metastable minimum and the saddle point.
 *
 * This function computes the barrier as the difference between the
 * local maximum and local minimum of the Landau free energy.
 *
 * Reference: Langer, Ann. Phys. 54, 258 (1969)
 * --------------------------------------------------------------------------- */
double nept_nucleation_barrier(LandauCoefficients lc, double reduced_temp)
{
    double a = lc.a;
    double b = lc.b;
    double c = lc.c;
    double t = reduced_temp;

    if (b >= -1e-12) return 0.0; /* no barrier for continuous transitions */
    if (c < 1e-12) return 0.0;

    double A = a * t;

    /* F(η) = (a*t/2)*η^2 + (b/4)*η^4 + (c/6)*η^6
     * Find extrema: dF/dη = η*(A + b*η^2 + c*η^4) = 0
     *
     * Solutions (non-zero): η^2 = [-b ± sqrt(b^2 - 4*A*c)]/(2*c) */
    double disc = b * b - 4.0 * A * c;
    if (disc <= 0.0) return 0.0; /* single minimum, no barrier */

    double sqrt_d = sqrt(disc);
    double eta2_small = (-b - sqrt_d) / (2.0 * c); /* local maximum (unstable) */
    double eta2_large = (-b + sqrt_d) / (2.0 * c); /* local minimum (metastable) */

    /* Only valid if both are positive */
    if (eta2_small <= 1e-15 || eta2_large <= 1e-15) return 0.0;

    (void)(eta2_small); /* eta_small = sqrt(eta2_small) */
    (void)(eta2_large); /* eta_large = sqrt(eta2_large) */

    /* Compute free energy at each extremum:
       F(eta^2) = 0.5*A*e2 + 0.25*b*e2^2 + (c/6)*e2^3 */
    double F_small = 0.5*A*eta2_small + 0.25*b*eta2_small*eta2_small
                     + (c/6.0)*eta2_small*eta2_small*eta2_small;
    double F_large = 0.5*A*eta2_large + 0.25*b*eta2_large*eta2_large
                     + (c/6.0)*eta2_large*eta2_large*eta2_large;
    double barrier = F_small - F_large;

    return (barrier > 0.0) ? barrier : 0.0;
}

/* ---------------------------------------------------------------------------
 * O(N) Model
 *
 * Knowledge point: For O(N) symmetric models, the order parameter is an
 * N-component vector η = (η₁, ..., η_N). The free energy is:
 *
 *   F = (a*t/2)*|η|^2 + (b/4)*(|η|^2)^2 - h*η₁
 *
 * Ground state: η aligned with external field h (longitudinal).
 * For N ≥ 2 with h = 0: spontaneous symmetry breaking yields N-1
 * massless Goldstone modes (transverse fluctuations).
 *
 * Equilibrium: dF/dη_i = 0 => (a*t + b*|η|^2)*η_i = h*δ_{i,1}
 * For h = 0, t < 0: |η_eq| = √(-a*t/b)
 *
 * Reference: Polyakov, Gauge Fields and Strings (1987), Chapter 1
 * --------------------------------------------------------------------------- */
double nept_on_equilibrium_magnitude(const OVectorModel *model,
                                      double reduced_temp)
{
    if (!model) return 0.0;
    double a = model->a_coefficient;
    double b = model->b_coefficient;
    double h = model->external_field;
    double t = reduced_temp;

    if (fabs(h) < 1e-15 && b > 1e-12) {
        if (t >= 0.0) return 0.0;
        return sqrt(-a * t / b);
    }

    /* With external field: solve cubic for |η| */
    if (fabs(h) > 1e-15 && b > 1e-12) {
        /* For small h, η ≈ h/(a*t) in disordered phase,
         * η ≈ √(-a*t/b) in ordered phase */
        if (t >= 0.0) {
            return h / (a * t + b * h * h); /* approximate */
        } else {
            double eta0 = sqrt(-a * t / b);
            return eta0 + h / (2.0 * a * fabs(t)); /* small-h expansion */
        }
    }

    return 0.0;
}

double nept_on_longitudinal_susceptibility(const OVectorModel *model,
                                            double reduced_temp)
{
    if (!model) return 0.0;
    /* χ_∥ = 1/(a*t) for t > 0, 1/(2*a*|t|) for t < 0 (mean-field) */
    double t = reduced_temp;
    double a = model->a_coefficient;

    if (fabs(t) < 1e-15) return 1e15;
    if (t > 0.0) return 1.0 / (a * t);
    return 1.0 / (2.0 * a * fabs(t));
}

double nept_on_transverse_susceptibility(const OVectorModel *model,
                                          double external_field)
{
    if (!model || fabs(external_field) < 1e-15) {
        /* For N >= 2 with h = 0 and t < 0: χ_⊥ = ∞ (Goldstone modes) */
        return 1e15;
    }

    /* χ_⊥ = η_eq / h (transverse response to small perpendicular field)
     * This diverges as h → 0 below Tc: the Hohenberg-Mermin-Wagner theorem
     * (d ≤ 2) actually destroys long-range order for continuous symmetry. */
    double eta = model->external_field; /* using h as proxy for η_eq */
    if (fabs(external_field) < 1e-15) return 1e15;
    return eta / external_field;
}

/* --- Convenience constructors --- */
LandauCoefficients nept_landau_ising_create(double a_val, double b_val)
{
    LandauCoefficients lc;
    lc.a = a_val;
    lc.b = b_val;
    lc.c = 0.0;
    lc.f0 = 0.0;
    lc.has_cubic = false;
    lc.d_cubic = 0.0;
    return lc;
}

LandauCoefficients nept_landau_tricritical_create(double a_val, double c_val)
{
    LandauCoefficients lc;
    lc.a = a_val;
    lc.b = 0.0;
    lc.c = c_val;
    lc.f0 = 0.0;
    lc.has_cubic = false;
    lc.d_cubic = 0.0;
    return lc;
}

bool nept_landau_is_at_tricritical(LandauCoefficients lc, double tolerance)
{
    return (fabs(lc.b) < tolerance && lc.c > tolerance);
}

/* ---------------------------------------------------------------------------
 * Hysteresis Loop Area
 *
 * Knowledge point: First-order phase transitions exhibit hysteresis:
 * the transition point depends on the direction of parameter change.
 * The area enclosed by the hysteresis loop represents energy dissipated
 * as heat during one cycle. This maps to the area inside the
 * spinodal region of the phase diagram.
 *
 * Simplified model: sweep η_eq(t) for increasing and decreasing t.
 * The area A = ∮ η dt gives the dissipated work.
 *
 * Reference: Binder, Rep. Prog. Phys. 50, 783 (1987)
 * --------------------------------------------------------------------------- */
double nept_landau_hysteresis_area(LandauCoefficients lc,
                                    double temp_min, double temp_max,
                                    int n_steps)
{
    if (n_steps < 2) return 0.0;
    if (lc.b >= -1e-12) return 0.0; /* No hysteresis for continuous transitions */

    double dt = (temp_max - temp_min) / (double)n_steps;
    double area = 0.0;
    double eta_prev = 0.0;

    /* Forward sweep (increasing t) */
    for (int i = 0; i <= n_steps; i++) {
        double t = temp_min + (double)i * dt;
        MeanFieldResults res;
        nept_landau_equilibrium(lc, t, &res);
        if (i > 0) {
            area += 0.5 * (eta_prev + res.spontaneous_op) * dt;
        }
        eta_prev = res.spontaneous_op;
    }

    /* Backward sweep (decreasing t) */
    double area_forward = area;
    area = 0.0;
    eta_prev = 0.0;

    for (int i = 0; i <= n_steps; i++) {
        double t = temp_max - (double)i * dt;
        MeanFieldResults res;
        nept_landau_equilibrium(lc, t, &res);
        if (i > 0) {
            area += 0.5 * (eta_prev + res.spontaneous_op) * dt;
        }
        eta_prev = res.spontaneous_op;
    }

    /* Hysteresis area = |forward - backward| */
    return fabs(area_forward - area);
}