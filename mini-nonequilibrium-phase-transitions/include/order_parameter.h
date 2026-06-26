#ifndef NEPT_ORDER_PARAMETER_H
#define NEPT_ORDER_PARAMETER_H

#include <stdbool.h>
#include <stddef.h>
#include "phase_transition.h"

/* ============================================================================
 * Order Parameters and Correlation Functions
 *
 * The order parameter eta is the central concept in phase transition theory.
 * It measures the degree of symmetry breaking: eta = 0 in the symmetric
 * (disordered) phase, eta != 0 in the broken-symmetry (ordered) phase.
 *
 * Landau insight (1937): near a continuous transition, the free energy can
 * be expanded in powers of the order parameter:
 *   F(eta) = F0 + a*t*eta^2 + b*eta^4 + ...
 * where t = (T - Tc)/Tc is reduced temperature, a > 0, b > 0.
 *
 * For first-order transitions, the eta^3 term MUST be present (b < 0),
 * or equivalently, a negative quartic with a positive sextic term.
 *
 * Key references:
 *   L.D. Landau — On the Theory of Phase Transitions (1937)
 *   V.L. Ginzburg & L.D. Landau — On the Theory of Superconductivity (1950)
 *   K.G. Wilson — Renormalization Group and Critical Phenomena (1971)
 * ============================================================================ */

/* --- Order Parameter Types --- */
typedef enum {
    OP_SCALAR_REAL = 0,     /* Real scalar (Ising magnetization, density diff) */
    OP_SCALAR_COMPLEX = 1,  /* Complex scalar (superfluid, BEC wavefunction)   */
    OP_VECTOR_O2 = 2,       /* 2-component vector = complex (XY model)         */
    OP_VECTOR_O3 = 3,       /* 3-component vector (Heisenberg magnetization)   */
    OP_TENSOR_RANK2 = 4,    /* Rank-2 tensor (nematic liquid crystal)          */
    OP_FOURIER_MODE = 5,    /* Fourier amplitude (pattern formation, stripes)  */
    OP_CHEMICAL_CONC = 6    /* Chemical concentration difference (BZ reaction) */
} OrderParameterType;

/* --- Scalar Order Parameter --- */
typedef struct {
    double magnitude;           /* |eta|, the scalar order parameter value     */
    double phase;               /* For complex OP, the phase angle             */
    OrderParameterType type;    /* Classification of this order parameter      */
    bool is_conserved;          /* Model A (non-conserved) vs Model B (conserved) */
    double relaxation_rate;     /* Gamma in TDGL: d(eta)/dt = -Gamma * dF/d(eta) */
} OrderParameter;

/* --- Two-point Correlation Function --- */
typedef struct {
    double *values;             /* G(r) values at discrete r_i                */
    double *r;                  /* Radial distances r_i                       */
    int n_points;               /* Number of sampled points                   */
    double correlation_length;  /* xi: G(r) ~ exp(-r/xi) for large r          */
    double exponent_eta_fisher; /* Fisher exponent eta (not OP!): G(r) ~ r^{-(d-2+eta)} at Tc */
} CorrelationFunction;

/* --- Structure Factor (Fourier transform of G(r)) --- */
typedef struct {
    double *intensity;          /* S(k) = <eta(k) eta(-k)>                   */
    double *k_values;           /* Wavevector magnitudes                     */
    int n_points;
    double ornstein_zernike_amplitude; /* S(k) ~ A/(xi^{-2} + k^2) for small k */
} StructureFactor;

/* --- Symmetry Breaking Descriptor --- */
typedef struct {
    NEPT_SymmetryType broken_symmetry;    /* What symmetry is broken          */
    NEPT_SymmetryType residual_symmetry;  /* What symmetry remains            */
    int order_of_broken_group;            /* |G| of full symmetry group       */
    int order_of_residual_group;          /* |H| of unbroken subgroup         */
    int num_goldstone_modes;             /* = dim(G/H), massless excitations  */
} SymmetryBreaking;

/* --- Time-Dependent Ginzburg-Landau (TDGL) Parameters --- */
typedef struct {
    double a_coefficient;       /* a in F = a*t*eta^2 + b*eta^4             */
    double b_coefficient;       /* b > 0 for stability                       */
    double kappa;               /* Gradient term: (kappa/2)*|grad eta|^2     */
    double kinetic_coefficient; /* Gamma: relaxation rate                    */
    double noise_amplitude;     /* D in Langevin: sqrt(2*D*Gamma)*xi(t)      */
} TDGLParameters;

/* --- Correlation length scaling --- */
typedef struct {
    double xi0_plus;            /* Amplitude above Tc: xi = xi0_plus * |t|^{-nu} */
    double xi0_minus;           /* Amplitude below Tc                             */
    double nu_exponent;         /* Critical exponent for correlation length       */
} CorrelationLengthScaling;

/* --- Core API --- */

/* Initialize a scalar (real) order parameter with given magnitude */
OrderParameter nept_op_scalar_create(double magnitude);

/* Initialize a complex order parameter */
OrderParameter nept_op_complex_create(double real_part, double imag_part);

/* TDGL relaxation step: d(eta)/dt = -Gamma * (2*a*t*eta + 4*b*eta^3 - kappa*laplacian) */
double nept_tdgl_deterministic_step(OrderParameter *op,
                                     const TDGLParameters *params,
                                     double reduced_temp,
                                     double laplacian,
                                     double dt);

/* Compute Landau free energy density for a given order parameter */
double nept_landau_free_energy(double eta,
                                const TDGLParameters *params,
                                double reduced_temp);

/* Compute equilibrium order parameter from minimization of Landau F */
double nept_landau_equilibrium_op(const TDGLParameters *params,
                                   double reduced_temp,
                                   NEPT_TransitionType *trans_type_out);

/* Determine if transition is first-order or continuous from Landau coefficients */
NEPT_TransitionType nept_landau_classify(const TDGLParameters *params);

/* Allocate a correlation function with n discrete points */
CorrelationFunction *nept_correlation_create(int n_points);

/* Set G(r) = exp(-r/xi) for a simple exponential decay model */
void nept_correlation_exponential(CorrelationFunction *cf,
                                   double xi,
                                   double r_min,
                                   double r_max);

/* Compute the correlation length xi from G(r) data by fitting log(G) vs r */
double nept_correlation_extract_xi(const CorrelationFunction *cf);

/* Free correlation function */
void nept_correlation_free(CorrelationFunction *cf);

/* Compute structure factor from correlation function via numerical Fourier transform */
StructureFactor *nept_structure_factor_from_correlation(
    const CorrelationFunction *cf, int n_kpoints);

void nept_structure_factor_free(StructureFactor *sf);

/* Create symmetry breaking descriptor */
SymmetryBreaking nept_symmetry_breaking_create(
    NEPT_SymmetryType broken, NEPT_SymmetryType residual,
    int order_full, int order_sub);

/* Count Goldstone modes: dim(G/H) = dim(G) - dim(H) */
int nept_count_goldstone_modes(const SymmetryBreaking *sb);

/* Estimate correlation length exponent nu from data */
double nept_estimate_nu_exponent(const double *xi_values,
                                  const double *reduced_temps,
                                  int n_points,
                                  double *xi0_out);

/* Test hyperscaling relation: 2 - alpha = d * nu (d = spatial dimension) */
bool nept_check_hyperscaling(double alpha, double nu, int dimension,
                              double tolerance);

#endif /* NEPT_ORDER_PARAMETER_H */