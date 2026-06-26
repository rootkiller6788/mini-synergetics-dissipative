#ifndef NEPT_CRITICAL_PHENOMENA_H
#define NEPT_CRITICAL_PHENOMENA_H

#include <stdbool.h>
#include "phase_transition.h"
#include "order_parameter.h"

/* ============================================================================
 * Critical Phenomena, Scaling, and Universality
 *
 * Near a continuous phase transition, thermodynamic quantities exhibit
 * power-law singularities characterized by critical exponents:
 *
 *   Order parameter:    eta ~ |t|^beta           (t < 0)
 *   Susceptibility:     chi ~ |t|^{-gamma}       (t > 0 or t < 0)
 *   Specific heat:      C ~ |t|^{-alpha}         (t > 0 or t < 0)
 *   Correlation length: xi ~ |t|^{-nu}           (both sides)
 *   Critical isotherm:  eta ~ |h|^{1/delta}      (t = 0)
 *   Correlation function: G(r) ~ r^{-(d-2+eta)}  (t = 0, Fisher exponent)
 *
 * Scaling relations (exact, independent of microscopic details):
 *   Rushbrooke:  alpha + 2*beta + gamma = 2
 *   Widom:       gamma = beta*(delta - 1)
 *   Fisher:      gamma = nu*(2 - eta)
 *   Josephson:   nu*d = 2 - alpha   (hyperscaling, d <= d_c)
 *
 * Universality: critical exponents depend only on:
 *   - Spatial dimension d
 *   - Symmetry of the order parameter (O(N))
 *   - Range of interactions (short vs. long)
 *
 * References:
 *   H.E. Stanley — Introduction to Phase Transitions... (1971)
 *   K.G. Wilson & J. Kogut — The Renormalization Group... (1974)
 *   M.E. Fisher — Rep. Prog. Phys. 30, 615 (1967)
 *   J. Cardy — Scaling and Renormalization... (1996)
 *   N. Goldenfeld — Lectures on Phase Transitions... (1992)
 * ============================================================================ */

/* --- Critical Exponents --- */
typedef struct {
    double alpha;       /* Specific heat: C ~ |t|^{-alpha}                    */
    double beta;        /* Order parameter: eta ~ |t|^beta (t < 0)            */
    double gamma;       /* Susceptibility: chi ~ |t|^{-gamma}                 */
    double delta;       /* Critical isotherm: eta ~ |h|^{1/delta} (t=0)       */
    double nu;          /* Correlation length: xi ~ |t|^{-nu}                 */
    double eta_fisher;  /* Fisher exponent: G(r) ~ r^{-(d-2+eta)} at Tc       */
    double z;           /* Dynamic exponent: tau ~ xi^z (relaxation time)     */
    int dimension;      /* Spatial dimension d                                */
    int n_components;   /* Number of OP components (O(N) model)               */
} CriticalExponents;

/* --- Predefined universality classes --- */
typedef enum {
    UC_MEAN_FIELD = 0,       /* d >= 4, any N: alpha=0, beta=1/2, gamma=1     */
    UC_ISING_2D = 1,         /* 2D Ising (exact): beta=1/8, gamma=7/4, nu=1   */
    UC_ISING_3D = 2,         /* 3D Ising: beta~0.326, gamma~1.239, nu~0.630   */
    UC_XY_3D = 3,            /* 3D XY: beta~0.349, gamma~1.318, nu~0.672      */
    UC_HEISENBERG_3D = 4,    /* 3D Heisenberg: beta~0.369, gamma~1.396, nu~0.711 */
    UC_PERCOLATION_2D = 5,   /* 2D percolation: beta=5/36, gamma=43/18        */
    UC_DIRECTED_PERCOLATION = 6, /* Non-equilibrium universal class           */
    UC_NONEQUILIBRIUM_DP = 7 /* Absorbing-state transitions (DP conjecture)   */
} UniversalityClass;

/* --- Scaling Function --- */
typedef struct {
    double *scaled_x;          /* x = t/|h|^{1/(beta+gamma)} or t*L^{1/nu}   */
    double *scaled_y;          /* y = eta/|h|^{beta/(beta+gamma)} or eta*L^{beta/nu} */
    int n_points;
    double scaling_power;      /* Exponent used in scaling collapse           */
} ScalingFunction;

/* --- Finite-Size Scaling --- */
typedef struct {
    int *system_sizes;             /* L values                                  */
    double *critical_point_shift;  /* Tc(L) - Tc(inf) ~ L^{-1/nu}              */
    double *peak_susceptibility;   /* chi_max(L) ~ L^{gamma/nu}                */
    double *peak_specific_heat;    /* C_max(L) ~ (L^{alpha/nu} - 1)/alpha      */
    int n_sizes;
    double estimate_nu;            /* Finite-size estimate of nu               */
    double estimate_gamma;         /* Finite-size estimate of gamma            */
} FiniteSizeScaling;

/* --- Data Collapse Quality Metric --- */
typedef struct {
    double residual;               /* Sum of squared deviations from master curve */
    double best_exponent;          /* Exponent that minimizes residual           */
    double confidence_interval;    /* 95% CI for the exponent estimate           */
} DataCollapseResult;

/* --- Renormalization Group (Simplified) --- */
typedef struct {
    int n_iterations;              /* Number of RG transformations              */
    double *coupling_constants;    /* Flow of coupling constants u_n            */
    int n_couplings;
    double *fixed_points;          /* Found fixed points of RG transformation   */
    int n_fixed_points;
    double relevant_eigenvalue;    /* lambda = b^{1/nu} > 1 at fixed point     */
    double irrelevant_eigenvalue;  /* |lambda| < 1                             */
    double estimate_nu;            /* nu = ln(b)/ln(lambda_relevant)            */
} RenormalizationGroup;

/* --- Non-Equilibrium Critical Exponents --- */
typedef struct {
    CriticalExponents exponents;
    double beta_prime;             /* Order parameter exponent for absorbing state */
    double nu_parallel;            /* Correlation time exponent: tau_|| ~ xi^z   */
    double nu_perpendicular;       /* Correlation length exponent                 */
    double theta;                  /* Initial slip exponent                       */
    bool dp_conjecture_holds;      /* Directed percolation universality check     */
} NonEquilibriumExponents;

/* --- Core API --- */

/* Create critical exponents for a given universality class */
CriticalExponents nept_critical_exponents_mean_field(int dimension);
CriticalExponents nept_critical_exponents_ising_2d(void);
CriticalExponents nept_critical_exponents_ising_3d(void);
CriticalExponents nept_critical_exponents_xy_3d(void);
CriticalExponents nept_critical_exponents_heisenberg_3d(void);

/* Check Rushbrooke scaling relation: alpha + 2*beta + gamma = 2 */
bool nept_check_rushbrooke(const CriticalExponents *ce, double tolerance);

/* Check Widom scaling relation: gamma = beta*(delta - 1) */
bool nept_check_widom(const CriticalExponents *ce, double tolerance);

/* Check Fisher scaling relation: gamma = nu*(2 - eta) */
bool nept_check_fisher(const CriticalExponents *ce, double tolerance);

/* Check Josephson hyperscaling: nu*d = 2 - alpha */
bool nept_check_josephson(const CriticalExponents *ce, double tolerance);

/* Check all four scaling relations simultaneously */
int nept_check_all_scaling_relations(const CriticalExponents *ce,
                                      double tolerance,
                                      bool *results_out);

/* Estimate beta exponent from order parameter vs reduced temperature data */
double nept_estimate_beta(const double *eta, const double *t_reduced,
                           int n_points, double *amplitude_out);

/* Estimate gamma exponent from susceptibility data */
double nept_estimate_gamma(const double *chi, const double *t_reduced,
                            int n_points, double *amplitude_out);

/* Estimate nu exponent from correlation length data */
double nept_estimate_nu_exponent(const double *xi, const double *t_reduced,
                                  int n_points, double *xi0_out);

/* Perform data collapse for order parameter: y = eta*|t|^{-beta} vs x = h*|t|^{-(beta+gamma)} */
DataCollapseResult nept_data_collapse_order_param(
    const double *eta, const double *t, const double *h,
    int n_points, double beta_guess, double gamma_guess);

/* Finite-size scaling analysis: extract nu and gamma from L-dependent data */
FiniteSizeScaling *nept_finite_size_scaling_create(int n_sizes);
void nept_finite_size_scaling_free(FiniteSizeScaling *fss);
void nept_finite_size_add_size(FiniteSizeScaling *fss, int idx,
                                int L, double Tc_shift,
                                double chi_peak, double C_peak);
int nept_finite_size_extract_exponents(FiniteSizeScaling *fss);

/* Simple real-space RG transformation for 1D Ising */
RenormalizationGroup *nept_rg_create(int max_iter);
void nept_rg_free(RenormalizationGroup *rg);
int nept_rg_iterate_1d_ising(RenormalizationGroup *rg,
                              double initial_coupling, int decimation_factor);

/* Compute the dynamic exponent z from relaxation time vs correlation length */
double nept_estimate_dynamic_exponent(const double *tau, const double *xi,
                                       int n_points);

/* Check if measured exponents match a known universality class */
UniversalityClass nept_match_universality_class(const CriticalExponents *ce,
                                                  double tolerance);

/* Compute upper critical dimension from mean-field exponents:
   d_c = 2*beta/gamma + 2*beta*n/(gamma) ...  d_c = 4 for phi^4 */
int nept_upper_critical_dimension(const CriticalExponents *mf_exponents);

#endif /* NEPT_CRITICAL_PHENOMENA_H */