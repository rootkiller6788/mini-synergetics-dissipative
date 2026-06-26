#ifndef OEP_PHASE_TRANSITION_H
#define OEP_PHASE_TRANSITION_H
#include "oep_core.h"

/* ==============================================================
 * Nonequilibrium Phase Transitions
 * Based on: Prigogine & Nicolis (1977), Cross & Hohenberg (1993)
 *
 * Phase transitions far from thermodynamic equilibrium are
 * characterized by symmetry breaking, bifurcation, and the
 * emergence of spatio-temporal order.
 *
 * Landau theory of phase transitions:
 *   F(ξ, p) = F₀ + a(p)ξ²/2 + bξ⁴/4 + ...   (Free energy)
 *
 * Classification:
 *   - Type I (first-order): discontinuous ξ at critical point
 *   - Type II (second-order): continuous ξ, dξ/dp diverges
 *   - Supercritical bifurcation: stable branch beyond p_c
 *   - Subcritical bifurcation: unstable branch, jump transition
 * ============================================================== */

/* ── Phase Transition Types ────────────────────────────── */

typedef enum {
    OEP_PT_NONE = 0,           /* No transition */
    OEP_PT_SUPERCRITICAL = 1,  /* Continuous (2nd order / supercritical pitchfork) */
    OEP_PT_SUBCRITICAL = 2,    /* Discontinuous (1st order / subcritical pitchfork) */
    OEP_PT_HOPF = 3,           /* Hopf bifurcation (limit cycle emergence) */
    OEP_PT_SADDLE_NODE = 4,    /* Saddle-node (fold) bifurcation */
    OEP_PT_TRANSCRITICAL = 5,  /* Transcritical (exchange of stability) */
    OEP_PT_PERIOD_DOUBLING = 6 /* Period doubling (route to chaos) */
} OEPPhaseTransitionType;

/* ── Bifurcation Analysis Structures ───────────────────── */

/* Bifurcation diagram point */
typedef struct {
    double param;          /* Control parameter value */
    double equilibrium;    /* Equilibrium order parameter */
    double eigenvalue;     /* Linear stability eigenvalue */
    bool   is_stable;      /* Stability of this branch */
} OEPBifPoint;

/* Bifurcation diagram (multiple branches over parameter range) */
typedef struct {
    OEPBifPoint** branches;  /* Array of branches (each is array of points) */
    int*          branch_lengths; /* Number of points per branch */
    int           n_branches;
    double        param_min;
    double        param_max;
} OEPBifDiag;

/* Phase transition descriptor */
typedef struct {
    OEPPhaseTransitionType type;
    double  critical_param;       /* p_c: critical value of control parameter */
    double  critical_exponent;    /* β: ξ ∝ |p-p_c|^β (order parameter exponent) */
    double  susceptibility_exponent; /* γ: χ ∝ |p-p_c|^{-γ} */
    double  correlation_exponent; /* ν: ξ ∝ |p-p_c|^{-ν} */
    double  heat_capacity_exponent; /* α̅: C ∝ |p-p_c|^{-α̅} */
    double  dynamical_exponent;   /* z: τ ∝ ξ^z */
    bool    is_continuous;        /* 2nd order (true) or 1st order (false) */
    double  jump_magnitude;       /* For 1st order: size of discontinuous jump */
    double  hysteresis_width;     /* For 1st order: width of bistable region */
} OEPPhaseTransition;

/* ── Bifurcation Diagram API ───────────────────────────── */

OEPBifDiag* oep_bifdiag_create(int n_branches);
void        oep_bifdiag_free(OEPBifDiag* bd);
void        oep_bifdiag_set_point(OEPBifDiag* bd, int branch, int idx,
                                   double param, double eq, double eig, bool stable);

/* Generate pitchfork bifurcation diagram:
 *   dξ/dt = αξ - βξ³,  α(p) = sensitivity * (p - p_c) */
OEPBifDiag* oep_bifdiag_pitchfork(double p_min, double p_max, int n_points,
                                   double p_c, double sensitivity, double beta);

/* Generate transcritical bifurcation diagram:
 *   dξ/dt = αξ - βξ² */
OEPBifDiag* oep_bifdiag_transcritical(double p_min, double p_max, int n_points,
                                       double p_c, double sensitivity, double beta);

/* Generate saddle-node bifurcation diagram:
 *   dξ/dt = α + βξ²  (fold at α = 0) */
OEPBifDiag* oep_bifdiag_saddle_node(double p_min, double p_max, int n_points,
                                     double p_c, double sensitivity, double beta);

/* Generate subcritical pitchfork:
 *   dξ/dt = αξ + βξ³ - γξ⁵  (β < 0, γ > 0 gives subcritical) */
OEPBifDiag* oep_bifdiag_subcritical(double p_min, double p_max, int n_points,
                                     double p_c, double sensitivity,
                                     double beta, double gamma);

/* Detect bifurcation type from numerical continuation */
OEPPhaseTransitionType oep_bifdiag_classify(const OEPBifDiag* bd);

/* ── Phase Transition Analysis API ─────────────────────── */

/* Detect phase transition from order parameter vs control parameter data */
OEPPhaseTransition* oep_pt_detect(const double* control_params,
                                   const double* order_params,
                                   int n_points);

void                oep_pt_free(OEPPhaseTransition* pt);

/* Compute critical exponents from data near the transition */
void                oep_pt_fit_critical_exponents(const double* control_params,
                                                   const double* order_params,
                                                   int n_points,
                                                   double p_c,
                                                   double* beta, double* gamma,
                                                   double* delta);

/* Landau free energy at given order parameter and control parameter */
double              oep_pt_landau_free_energy(double xi, double p, double p_c,
                                               double a_coeff, double b_coeff);

/* Metastability region bounds (for 1st order transitions) */
void                oep_pt_metastability_bounds(const OEPPhaseTransition* pt,
                                                 double* lower, double* upper);

/* Maxwell construction: equal-area rule for 1st order transitions */
double              oep_pt_maxwell_construction(const double* control_params,
                                                 const double* order_params,
                                                 int n_points);

/* ── Critical Phenomena ────────────────────────────────── */

/* Scaling relations between critical exponents:
 *   Rushbrooke: α̅ + 2β + γ ≥ 2
 *   Widom:      γ = β(δ - 1)
 *   Fisher:     ν(2 - η) = γ */
bool                oep_pt_check_rushbrooke(double alpha, double beta, double gamma);
bool                oep_pt_check_widom(double gamma, double beta, double delta);
bool                oep_pt_check_fisher(double nu, double eta, double gamma);

/* Finite-size scaling:
 *   ξ_L(p) = L^{β/ν} * f(|p-p_c| * L^{1/ν}) */
double              oep_pt_finite_size_scaling(double L, double p, double p_c,
                                                double beta, double nu);

/* Universality class comparison */
bool                oep_pt_same_universality(const OEPPhaseTransition* a,
                                              const OEPPhaseTransition* b,
                                              double tolerance);

/* ── Nonequilibrium Criteria ───────────────────────────── */

/* Excess entropy production at phase transition:
 *   ΔP = ∫ dV [σ(p_c + ε) - σ(p_c - ε)] */
double              oep_pt_excess_entropy_production(const double* entropy_prod,
                                                      const double* control_params,
                                                      int n_points, double p_c);

/* Fluctuation enhancement near critical point:
 *   <δξ²> ∝ |p - p_c|^{-γ} */
double              oep_pt_fluctuation_enhancement(double p, double p_c, double gamma);

/* Critical slowing down:
 *   τ_corr ∝ |p - p_c|^{-νz} */
double              oep_pt_correlation_time(double tau0, double p, double p_c,
                                             double nu, double z);

/* Ginzburg criterion: region where mean-field theory breaks down */
double              oep_pt_ginzburg_criterion(double xi0, double p, double p_c,
                                               double nu, int dimension);

#endif
