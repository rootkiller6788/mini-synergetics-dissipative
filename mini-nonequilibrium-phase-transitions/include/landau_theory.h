#ifndef NEPT_LANDAU_THEORY_H
#define NEPT_LANDAU_THEORY_H

#include <stdbool.h>
#include "phase_transition.h"
#include "order_parameter.h"

/* ============================================================================
 * Landau-Ginzburg Theory of Phase Transitions
 *
 * Landau theory is the mean-field theory of phase transitions. The core
 * idea is to expand the free energy F as a power series in the order
 * parameter eta, constrained by the symmetry of the disordered phase.
 *
 * For a scalar OP with Z2 symmetry (eta -> -eta):
 *   F = F0 + (a*t/2)*eta^2 + (b/4)*eta^4 + (c/6)*eta^6
 *   with t = (T - Tc)/Tc, a > 0, b can be +/- , c > 0 (stability).
 *
 * When b > 0: continuous transition at t = 0.
 * When b < 0: first-order transition when F(eta=0) = F(eta=eta*).
 *
 * Ginzburg criterion: Mean-field breaks down when fluctuations dominate.
 *   t_G = (kB * Tc) / (a0 * xi0^3 * DeltaC)^2
 * Below this |t|, critical exponents deviate from mean-field values.
 *
 * References:
 *   L.D. Landau — Zh. Eksp. Teor. Fiz. 7, 19 (1937)
 *   V.L. Ginzburg — Sov. Phys. Solid State 2, 1824 (1960)
 *   P.C. Hohenberg & B.I. Halperin — Rev. Mod. Phys. 49, 435 (1977)
 * ============================================================================ */

/* --- Landau Expansion Coefficients --- */
typedef struct {
    double a;       /* Coefficient of eta^2/2: a > 0, a*t changes sign at Tc */
    double b;       /* Coefficient of eta^4/4: b > 0 => continuous,
                       b < 0 => first-order transition                        */
    double c;       /* Coefficient of eta^6/6: c > 0 for stability            */
    double f0;      /* Background free energy F0(T), analytic at Tc           */
    bool has_cubic; /* Whether cubic term (eta^3) is allowed by symmetry      */
    double d_cubic; /* Coefficient of eta^3/3 (if symmetry allows)            */
} LandauCoefficients;

/* --- Mean-Field Thermodynamic Quantities --- */
typedef struct {
    double spontaneous_op;         /* Equilibrium order parameter eta_eq(t)     */
    double susceptibility;         /* chi = d(eta)/d(h)|_{h=0} ~ |t|^{-gamma} */
    double specific_heat;          /* C = -T * d^2F/dT^2                       */
    double specific_heat_jump;     /* Discontinuity Delta C at Tc (mean-field) */
    double free_energy;            /* F(eta_eq)                                */
    double free_energy_barrier;    /* Delta F between minima (first-order)     */
} MeanFieldResults;

/* --- Spinodal and Binodal (First-Order Transitions) --- */
typedef struct {
    double binodal_t;              /* Coexistence curve: F(min1) = F(min2)     */
    double spinodal_upper;         /* Limit of metastability (overheating)     */
    double spinodal_lower;         /* Limit of metastability (undercooling)    */
    double latent_heat;            /* T * Delta S across first-order transition */
    double coexistence_width;      /* |eta_ordered - eta_disordered| at binodal */
} FirstOrderCoexistence;

/* --- Ginzburg Criterion --- */
typedef struct {
    double ginzburg_reduced_temp;  /* t_G: temperature below which MF fails    */
    double ginzburg_length;        /* xi_G: correlation length at t_G          */
    int upper_critical_dimension;  /* d_c = 4 for Ising (scalar phi^4)         */
    bool mean_field_valid;         /* True for |t| > t_G or d >= d_c           */
} GinzburgCriterion;

/* --- Multi-Component Order Parameter Extension --- */
typedef struct {
    int n_components;              /* Number of OP components (O(N) model)     */
    double *eta_components;        /* Array of OP component values             */
    double a_coefficient;          /* (a*t/2)*|eta|^2                          */
    double b_coefficient;          /* (b/4)*(|eta|^2)^2                        */
    double c_coefficient;          /* (c/6)*(|eta|^2)^3                        */
    double external_field;         /* h couples to eta_1 (symmetry-breaking)   */
} OVectorModel;

/* --- Core API --- */

/* Compute mean-field solution for scalar Landau theory
 *   F = f0 + (a*t/2)*eta^2 + (b/4)*eta^4 + (c/6)*eta^6
 *   Minimize dF/d(eta) = 0 to find equilibrium OP.
 *   Returns the equilibrium order parameter eta_eq(t). */
double nept_landau_equilibrium(LandauCoefficients lc,
                                double reduced_temp,
                                MeanFieldResults *results);

/* Determine the transition type from Landau coefficients alone */
NEPT_TransitionType nept_landau_transition_type(LandauCoefficients lc);

/* Compute the susceptibility chi = d(eta_eq)/dh|_h=0 in mean-field */
double nept_landau_susceptibility(LandauCoefficients lc,
                                   double reduced_temp);

/* Compute the specific heat jump at Tc purely from Landau parameters */
double nept_landau_specific_heat_jump(LandauCoefficients lc, double Tc);

/* Find the first-order coexistence (binodal) temperature
 *   where F(eta=0) = F(eta=eta*)
 *   Returns reduced temperature t_binodal. */
double nept_landau_binodal_temperature(LandauCoefficients lc,
                                        FirstOrderCoexistence *coex);

/* Compute spinodal points: where d^2F/d(eta)^2 = 0 (limit of metastability) */
void nept_landau_spinodal_points(LandauCoefficients lc,
                                  double *spinodal_upper,
                                  double *spinodal_lower);

/* Compute the Ginzburg criterion for mean-field validity
 *   t_G = (kB * Tc / (a0 * xi0^3 * DeltaC))^2
 *   d_c = 2 * beta/gamma + 2 = 4 for Landau (beta=1/2, gamma=1) */
GinzburgCriterion nept_ginzburg_criterion(double a0, double xi0,
                                           double delta_C, double Tc,
                                           int spatial_dimension);

/* Free energy barrier for nucleation in first-order transitions */
double nept_nucleation_barrier(LandauCoefficients lc, double reduced_temp);

/* O(N) model: compute equilibrium magnitude of n-component order parameter */
double nept_on_equilibrium_magnitude(const OVectorModel *model,
                                      double reduced_temp);

/* O(N) model: compute susceptibility (longitudinal) */
double nept_on_longitudinal_susceptibility(const OVectorModel *model,
                                            double reduced_temp);

/* O(N) model: compute transverse susceptibility (Goldstone modes for N >= 2) */
double nept_on_transverse_susceptibility(const OVectorModel *model,
                                          double external_field);

/* Construct standard Landau coefficients for Ising universality class */
LandauCoefficients nept_landau_ising_create(double a_val, double b_val);

/* Construct standard phi^6 Landau coefficients for tricritical point */
LandauCoefficients nept_landau_tricritical_create(double a_val, double c_val);

/* Compute the tricritical point condition: a = 0, b = 0, c > 0 */
bool nept_landau_is_at_tricritical(LandauCoefficients lc, double tolerance);

/* Compute hysteresis loop area for first-order transition (proxy for dissipation) */
double nept_landau_hysteresis_area(LandauCoefficients lc,
                                    double temp_min, double temp_max,
                                    int n_steps);

#endif /* NEPT_LANDAU_THEORY_H */