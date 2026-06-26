#ifndef OEP_ORDER_PARAMETER_H
#define OEP_ORDER_PARAMETER_H
#include "oep_core.h"

/* ==============================================================
 * Order Parameter Structures and Dynamics
 * Based on: Haken (1977) "Synergetics", Haken (1983) "Advanced Synergetics"
 *
 * An order parameter ξ is a macroscopic variable that emerges from
 * the cooperative behavior of many microscopic degrees of freedom.
 * Near a critical point, the dynamics are governed by:
 *
 *   dξ/dt = α(p)ξ - βξ³ + √(Q) Γ(t)
 *
 * where α(p) = (pc - p)/pc controls the bifurcation parameter,
 * β > 0 ensures stability, and Γ(t) is white noise with intensity Q.
 * ============================================================== */

/* ── Order Parameter Type ──────────────────────────────── */

/* Scalar order parameter (single real-valued) */
typedef struct {
    double value;           /* Current magnitude */
    double critical_value;  /* Value at phase transition */
    double alpha;           /* Linear growth/decay coefficient α(p) */
    double beta;            /* Cubic saturation coefficient (> 0 for stability) */
    double noise_intensity; /* Q: strength of additive white noise */
    double relaxation_time; /* τ: characteristic relaxation time */
    bool   is_active;       /* Whether this OP is currently "enslaving" */
} OEPScalarOP;

/* Complex order parameter (for oscillatory instabilities) */
typedef struct {
    double amplitude;
    double phase;
    double frequency;       /* ωc: critical frequency */
    double growth_rate;     /* σ: linear growth/decay */
    double nonlinear_coeff; /* g: Landau constant (complex) */
    double noise_real;
    double noise_imag;
} OEPComplexOP;

/* Multi-component order parameter vector */
typedef struct {
    OEPVector*   values;        /* ξ₁, ξ₂, ..., ξₙ */
    OEPMatrix*   coupling;      /* n×n coupling matrix between OP components */
    OEPVector*   linear_coeffs; /* α₁,...,αₙ */
    OEPVector*   cubic_coeffs;  /* β₁,...,βₙ */
    int          n_components;
    double       noise_correlation; /* Shared noise correlation */
} OEPMultiOP;

/* Order parameter identification result */
typedef struct {
    double   significance;     /* Fraction of variance explained (0..1) */
    double   eigenvalue;       /* Covariance eigenvalue */
    OEPVector* direction;      /* Principal direction (loading vector) */
    double   relaxation_time;  /* τ = 1/|λ| from linear dynamics */
    double   nonlinear_coeff;  /* β from cubic fit */
    bool     is_critical;      /* Near bifurcation point */
} OEPIdentifiedOP;

/* ── Scalar Order Parameter API ────────────────────────── */

OEPScalarOP* oep_op_scalar_create(double alpha, double beta, double noise);
void         oep_op_scalar_free(OEPScalarOP* op);
void         oep_op_scalar_set_critical(OEPScalarOP* op, double pc, double p);
void         oep_op_scalar_set_bifurcation(OEPScalarOP* op, double control_param,
                                           double critical_point, double sensitivity);

/* Order parameter equation right-hand side:
 *   dξ/dt = α*ξ - β*ξ³ + noise
 * Used as callback for ODE solvers. */
double       oep_op_scalar_rhs(double xi, double alpha, double beta, double noise);

/* Deterministic drift: f(ξ) = αξ - βξ³ */
double       oep_op_scalar_drift(double xi, double alpha, double beta);

/* Potential function: V(ξ) = -αξ²/2 + βξ⁴/4 */
double       oep_op_scalar_potential(double xi, double alpha, double beta);

/* Equilibrium points of dξ/dt = αξ - βξ³ */
void         oep_op_scalar_equilibria(double alpha, double beta,
                                      double* eq0, double* eq1, double* eq2, int* n_equilibria);

/* Stability analysis at equilibrium */
double       oep_op_scalar_stability(double xi_eq, double alpha, double beta);
bool         oep_op_scalar_is_stable(double xi_eq, double alpha, double beta);

/* Relaxation time τ = -1/λ (λ = eigenvalue of linearized dynamics) */
double       oep_op_scalar_relaxation_time(double alpha, double xi_eq);

/* Susceptibility χ = ∂ξ/∂f |_{ext} */
double       oep_op_scalar_susceptibility(double alpha, double beta, double xi_eq);

/* Susceptibility divergence near criticality: χ ∝ |p-pc|^{-γ} */
double       oep_op_scalar_critical_exponent_susceptibility(double alpha);

/* Correlation length proxy from order parameter */
double       oep_op_scalar_correlation_length(double xi, double xi_sq_mean);

/* Time evolution via Runge-Kutta 4 */
void         oep_op_scalar_evolve(OEPScalarOP* op, double dt, int n_steps,
                                  double* trajectory);

/* Stochastic evolution (Euler-Maruyama):
 *   ξ_{n+1} = ξ_n + dt*(αξ_n - βξ_n³) + √(Q*dt) * N(0,1) */
void         oep_op_scalar_evolve_stochastic(OEPScalarOP* op, double dt, int n_steps,
                                             double* trajectory);

/* Stationary probability distribution (Fokker-Planck):
 *   P_s(ξ) ∝ exp(-2V(ξ)/Q) */
double       oep_op_scalar_stationary_pdf(double xi, double alpha, double beta, double Q);

/* Most probable amplitude (mode of stationary distribution) */
double       oep_op_scalar_most_probable(double alpha, double beta, double Q);

/* ── Complex Order Parameter API ───────────────────────── */

OEPComplexOP* oep_op_complex_create(double growth_rate, double frequency,
                                     double nonlinear_coeff, double noise_real,
                                     double noise_imag);
void          oep_op_complex_free(OEPComplexOP* op);

/* Complex Ginzburg-Landau:
 *   dA/dt = σA - g|A|²A
 *   with σ = growth_rate + i*frequency, g = nonlinear_coeff */
void          oep_op_complex_rhs(double amplitude, double phase, double sigma_re,
                                 double sigma_im, double g_re, double g_im,
                                 double* dA_re, double* dA_imag);
void          oep_op_complex_evolve(OEPComplexOP* op, double dt, int n_steps,
                                    double* re_traj, double* im_traj);

/* Limit cycle amplitude: |A|_eq = √(σ/g) */
double        oep_op_complex_equilibrium_amplitude(double sigma_re, double g_re);

/* Phase dynamics near limit cycle */
double        oep_op_complex_phase_velocity(double sigma_im, double g_im,
                                            double amplitude);

/* ── Multi-component Order Parameter API ───────────────── */

OEPMultiOP*   oep_op_multi_create(int n_components);
void          oep_op_multi_free(OEPMultiOP* op);
void          oep_op_multi_set_linear(OEPMultiOP* op, int i, double alpha);
void          oep_op_multi_set_cubic(OEPMultiOP* op, int i, double beta);
void          oep_op_multi_set_coupling(OEPMultiOP* op, int i, int j, double val);

/* Multi-component dynamics:
 *   dξᵢ/dt = αᵢξᵢ - βᵢξᵢ³ + Σⱼ C_{ij}ξⱼ + noise */
void          oep_op_multi_rhs(const OEPMultiOP* op, const OEPVector* state,
                               OEPVector* deriv);
void          oep_op_multi_evolve(OEPMultiOP* op, double dt, int n_steps,
                                  OEPMatrix* trajectory);

/* Normal modes of the linearized multi-component system */
void          oep_op_multi_normal_modes(const OEPMultiOP* op, OEPVector* eigenvalues,
                                        OEPMatrix* eigenvectors);

/* Adiabatic elimination: reduce n-component to m-component (m < n)
 * by eliminating fast-relaxing modes. Returns the reduced order. */
int           oep_op_multi_adiabatic_reduce(OEPMultiOP* op, double timescale_threshold,
                                            OEPVector* remaining_indices);

/* ── Order Parameter Identification API ─────────────────── */

/* Identify candidate order parameters from multivariate time series data.
 *   data: N_obs × n_vars matrix (row-major), each row is one observation
 *   returns array of identified OPs (caller frees each with oep_op_id_free) */
OEPIdentifiedOP** oep_op_identify_from_data(const OEPMatrix* data,
                                            int* n_identified);
void              oep_op_id_free(OEPIdentifiedOP* op);

/* Principal component analysis for order parameter discovery */
void              oep_op_pca(const OEPMatrix* data, OEPVector* eigenvalues,
                             OEPMatrix* loadings, OEPVector* explained_variance);

/* Test if a collective variable qualifies as order parameter:
 *   (a) large variance fraction
 *   (b) slow relaxation compared to other modes
 *   (c) symmetry-breaking near criticality */
bool              oep_op_is_valid_order_parameter(const OEPIdentifiedOP* op,
                                                  double var_threshold,
                                                  double relaxation_ratio);

/* Compute the circular causality measure: how much the OP
 * both enslaves and is supported by microscopic variables */
double            oep_op_circular_causality(const OEPMatrix* data,
                                            const OEPVector* op_direction,
                                            const OEPMatrix* microscopic);

#endif
