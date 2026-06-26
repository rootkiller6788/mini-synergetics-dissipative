#ifndef HAKEN_ORDER_PARAM_H
#define HAKEN_ORDER_PARAM_H
#include "haken_core.h"

/* ==============================================================
 * haken_order_param.h — Order Parameter Identification & Analysis
 *
 * The order parameter is the key concept in Haken's synergetics.
 * It describes the macroscopic state that emerges from cooperative
 * behavior of microscopic components.
 *
 * Identification methods:
 *   1. Spectral: modes with near-zero eigenvalues
 *   2. Statistical: principal component analysis (PCA) of fluctuations
 *   3. Physical: amplitude of critical Fourier modes
 *   4. Ensemble: average over microscopic configurations
 *
 * Ginzburg-Landau theory: near critical point, free energy
 *   F(ξ) = F₀ + (a/2)ξ² + (b/4)ξ⁴
 * with equation of motion:
 *   dξ/dt = -∂F/∂ξ = -a·ξ - b·ξ³
 *
 * Reference:
 *   Haken (1977) Synergetics, §2, §5, §6
 *   Cross & Hohenberg (1993) Rev. Mod. Phys. 65, 851
 *   Ma (1976) Modern Theory of Critical Phenomena
 * ============================================================== */

/* -----------------------------------------------------------------
 * L1: Order parameter type definitions
 * ----------------------------------------------------------------- */

/** Types of order parameter identification */
typedef enum {
    HAKEN_OP_SPECTRAL       = 0,  /* Eigenmode amplitude */
    HAKEN_OP_PCA            = 1,  /* Principal component analysis */
    HAKEN_OP_MEAN_FIELD     = 2,  /* Average over all variables */
    HAKEN_OP_MODE_AMPLITUDE = 3,  /* Amplitude of critical Fourier mode */
    HAKEN_OP_CUSTOM         = 4
} Haken_OPMethod;

/** Ginzburg-Landau free energy parameters */
typedef struct {
    double   a;              /* Linear coefficient: a ∝ (α_c - α) */
    double   b;              /* Cubic coefficient (must be > 0 for stability) */
    double   critical_alpha; /* Critical control parameter value */
    double   xi0;            /* Bare correlation length */
    double   tau0;           /* Bare relaxation time */
    double   nu;             /* Correlation length exponent */
    double   z;              /* Dynamic critical exponent */
    double   beta;           /* Order parameter critical exponent */
} Haken_GLParameters;

/** Order parameter time series for fluctuation analysis */
typedef struct {
    int      n_samples;
    double*  data;           /* Time series, length n_samples */
    double*  power_spectrum; /* FFT power spectrum */
    double   mean;
    double   variance;
    double   skewness;       /* Third moment / σ³ */
    double   kurtosis;       /* Fourth moment / σ⁴ - 3 (excess) */
    double   autocorr_time;  /* Integrated autocorrelation time */
    double   critical_exponent; /* Estimate of β from scaling */
} Haken_OPTimeSeries;

/* -----------------------------------------------------------------
 * L5: Order parameter identification API
 * ----------------------------------------------------------------- */

/** Identify order parameters from spectral decomposition.
 *  Finds modes with eigenvalues closest to zero.
 *  @param sys     System with computed eigenmodes
 *  @param n_ops   Number of order parameters to identify
 *  @param eps     Real-part threshold for "near-zero"
 *  @param result  Output: filled order parameter descriptor */
void haken_identify_spectral(const Haken_System* sys, int n_ops,
                              double eps, Haken_OrderParam* result);

/** Identify order parameters from PCA of a trajectory.
 *  Uses covariance matrix of state fluctuations.
 *  @param trajectory  States over time, shape (n_steps × n_dim)
 *  @param n_steps     Number of time points
 *  @param n_dim       Dimension of state space
 *  @param n_ops       Number of OPs to extract
 *  @param result      Output: order parameters */
void haken_identify_pca(const double* trajectory, int n_steps,
                         int n_dim, int n_ops, Haken_OrderParam* result);

/** Compute the mean-field order parameter:
 *  ξ = (1/N) Σ_i q_i  (simple average)
 *  Useful for globally coupled systems. */
double haken_mean_field_op(const double* q, int n);

/** Compute the dominant Fourier mode amplitude:
 *  ξ_k = |(1/N) Σ_j q_j · exp(-2πijk/N)|²
 *  For pattern-forming systems on a lattice.
 *  @param q       Field values on N lattice sites
 *  @param n       Number of sites
 *  @param k       Wavenumber of interest (0 for uniform, 1 for first harmonic) */
double haken_fourier_mode_amplitude(const double* q, int n, int k);

/** Compute the spatial correlation function:
 *  C(r) = <q(x+r)·q(x)> / <q(x)²>
 *  @param field   1D field values
 *  @param n       Number of spatial points
 *  @param r       Maximum lag
 *  @param corr    Output: correlation values, length r+1 */
void haken_spatial_correlation(const double* field, int n, int r,
                                double* corr);

/** Estimate the correlation length from spatial correlation:
 *  C(r) ≈ exp(-r/ξ), fit log(C(r)) vs r. */
double haken_correlation_length(const double* correlation, int r_max);

/* -----------------------------------------------------------------
 * Ginzburg-Landau analysis
 * ----------------------------------------------------------------- */

/** Fit Ginzburg-Landau parameters from order parameter time series.
 *  The GL equation: dξ/dt = -a·ξ - b·ξ³
 *  Fit via least squares: minimize Σ_t (dξ/dt + aξ + bξ³)²
 *
 *  @param ts       Time series of order parameter
 *  @param dt       Time step between samples
 *  @param params   Output: fitted GL parameters */
void haken_fit_gl_parameters(const Haken_OPTimeSeries* ts, double dt,
                              Haken_GLParameters* params);

/** Compute the Ginzburg-Landau free energy:
 *  F(ξ) = -ʃ (a·ξ + b·ξ³) dξ = -(a/2)ξ² - (b/4)ξ⁴
 *  (sign convention: dξ/dt = +∂F/∂ξ for maximization in synergetics) */
double haken_gl_free_energy(double xi, const Haken_GLParameters* params);

/** Compute the GL potential landscape: F(ξ) for a range of ξ.
 *  Useful for visualizing attractors. */
void haken_gl_potential_landscape(const Haken_GLParameters* params,
                                   double xi_min, double xi_max,
                                   int n_points,
                                   double* xi_grid, double* F_grid);

/** Compute the stationary solutions of the GL equation:
 *  dξ/dt = 0 → -a·ξ - b·ξ³ = 0
 *  Returns number of real fixed points and their values.
 *  For a<0 (above threshold): 3 solutions (0, ±√(-a/b))
 *  For a>0 (below threshold): 1 solution (0)
 *  @param n_fixed  Output: number of fixed points (1 or 3)
 *  @param fixed_xi Output: fixed point values, at least 3 elements */
int haken_gl_fixed_points(const Haken_GLParameters* params,
                           double* fixed_xi);

/** Determine stability of GL fixed points.
 *  Linear stability: λ = -a - 3b·ξ₀²
 *  For ξ=0: λ = -a (stable if a>0, unstable if a<0)
 *  For ξ=±√(-a/b): λ = -a - 3b·(-a/b) = -a + 3a = 2a (stable if a<0)
 *  Returns array of stability indicators (negative = stable, positive = unstable) */
void haken_gl_stability(const Haken_GLParameters* params,
                         const double* fixed_xi, int n_fixed,
                         double* stability);

/** Compute critical exponents from data:
 *  ξ_eq ~ |α_c - α|^β
 *  τ ~ |α_c - α|^{-νz}
 *  χ ~ |α_c - α|^{-γ}
 *  Fits power-law to data from parameter sweep. */
void haken_critical_exponents(const Haken_ParameterSweep* sweep,
                               Haken_GLParameters* params);

/* -----------------------------------------------------------------
 * Order parameter validation
 * ----------------------------------------------------------------- */

/** Check if a variable qualifies as a bona fide order parameter:
 *  1. Vanishes in disordered phase (α < α_c)
 *  2. Grows continuously from zero (supercritical) or jumps (subcritical)
 *  3. Satisfies GL equation near threshold
 *  4. All other variables are slaved to it */
bool haken_validate_order_parameter(const Haken_System* sys,
                                     int candidate_mode,
                                     double tolerance);

/** Compute the fraction of total variance explained by
 *  the order parameters. Should be close to 1 for valid slaving. */
double haken_variance_explained(const Haken_System* sys,
                                 const Haken_OrderParam* op);

/** Compute power spectrum of a time series via periodogram.
 *  ps[k] = |Σ_t data[t]·exp(-2πikt/n)|² / n, for k in [0, n/2]. */
void haken_power_spectrum(double* data, int n, double* ps);

#endif /* HAKEN_ORDER_PARAM_H */
