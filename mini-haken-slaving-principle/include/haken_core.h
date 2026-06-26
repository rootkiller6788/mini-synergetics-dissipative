#ifndef HAKEN_CORE_H
#define HAKEN_CORE_H
#include <stdbool.h>
#include <stddef.h>

/* ==============================================================
 * haken_core.h — Core Types for Haken's Slaving Principle
 *
 * Hermann Haken (1975, 1983) "Synergetics: An Introduction"
 * The slaving principle (Versklavungsprinzip) states that near an
 * instability, the dynamics of a complex system are governed by a
 * few slow "order parameters" which enslave the fast-relaxing modes.
 *
 * Mathematical formulation:
 *   dq/dt = L(α)q + N(q,α)
 *   where α is the control parameter.
 *   Near α = α_c:  λ₁ ≈ 0 (order parameter mode)
 *                   λ_k < -γ < 0 (stable modes, k ≥ 2)
 *
 * Adiabatic elimination: For fast modes (k ≥ 2):
 *   0 ≈ λ_k q_k + N_k(q₁,...,q_N)  →  q_k = f_k(q₁)
 *
 * Closed order parameter equation:
 *   dq₁/dt = λ₁ q₁ + Ñ₁(q₁)
 *
 * References:
 *   Haken (1977) "Synergetics — An Introduction", Springer
 *   Haken (1983) "Advanced Synergetics", Springer
 *   Cross & Hohenberg (1993) Rev. Mod. Phys. 65, 851
 * ============================================================== */

/* -----------------------------------------------------------------
 * L1: Core type definitions
 * ----------------------------------------------------------------- */

/** Stability type of an eigenmode */
typedef enum {
    HAKEN_MODE_STABLE       = 0,  /* Re(λ) < -γ < 0, fast relaxing */
    HAKEN_MODE_CRITICAL     = 1,  /* Re(λ) ≈ 0, order parameter */
    HAKEN_MODE_UNSTABLE     = 2,  /* Re(λ) > 0, growing */
    HAKEN_MODE_MARGINAL     = 3   /* Re(λ) = 0 exactly */
} Haken_ModeType;

/** Nonlinear coupling type between modes */
typedef enum {
    HAKEN_COUPLING_QUADRATIC    = 0,  /* N(q) ~ q² */
    HAKEN_COUPLING_CUBIC        = 1,  /* N(q) ~ q³ (saturating) */
    HAKEN_COUPLING_QUARTIC      = 2,  /* N(q) ~ q⁴ */
    HAKEN_COUPLING_GENERAL      = 3   /* User-defined */
} Haken_CouplingType;

/** Bifurcation type classification */
typedef enum {
    HAKEN_BIF_PITCHFORK     = 0,  /* dq/dt = α q - q³ */
    HAKEN_BIF_HOPF          = 1,  /* Complex conjugate pair crosses axis */
    HAKEN_BIF_SADDLE_NODE   = 2,  /* Two fixed points collide */
    HAKEN_BIF_TRANSCRITICAL = 3,  /* Stability exchange */
    HAKEN_BIF_NONE          = 4
} Haken_BifurcationType;

/* -----------------------------------------------------------------
 * Mode descriptor: eigenvalues, eigenvectors, classification
 * ----------------------------------------------------------------- */

typedef struct {
    double   lambda_re;      /* Real part of eigenvalue */
    double   lambda_im;      /* Imaginary part (oscillatory modes) */
    double*  eigenvector_re;  /* Right eigenvector, real part */
    double*  eigenvector_im;  /* Right eigenvector, imag part */
    double*  left_eigenv_re;  /* Left eigenvector (adjoint) */
    double*  left_eigenv_im;
    double   relaxation_time; /* τ = -1/Re(λ), for stable modes */
    double   oscillation_period; /* T = 2π/Im(λ) */
    Haken_ModeType type;
    int      index;           /* Original index in state vector */
} Haken_Mode;

/* -----------------------------------------------------------------
 * Full system descriptor
 * ----------------------------------------------------------------- */

typedef struct {
    int      n_dim;            /* Total system dimension N */
    int      n_critical;       /* Number of critical (order param) modes */
    int      n_stable;         /* Number of stable (slaved) modes */
    int      n_unstable;       /* Number of unstable modes */

    double*  state;            /* Current state vector q, length n_dim */
    double   control_param;    /* Control parameter α */
    double   critical_param;   /* Critical value α_c */

    double*  linear_matrix;    /* Linear stability matrix L(α), n_dim×n_dim */
    Haken_Mode* modes;         /* Eigenmode array, length n_dim */

    int*     critical_indices; /* Indices of critical modes */
    int*     stable_indices;   /* Indices of stable (fast) modes */
    int*     unstable_indices; /* Indices of unstable modes */

    /* Nonlinear coupling function: N(q) evaluated at given state */
    void (*nonlinear_coupling)(const double* q, int n, double* Nq);
    Haken_CouplingType coupling_type;

    /* External parameter dependence: L(α) as function of α */
    void (*linear_matrix_func)(double alpha, int n, double* L_out);

    void*    context;          /* User auxiliary data */
} Haken_System;

/* -----------------------------------------------------------------
 * Order parameter descriptor
 * ----------------------------------------------------------------- */

typedef struct {
    int      n_op;             /* Dimension of order parameter space */
    double*  amplitude;        /* Order parameter amplitudes, length n_op */
    double*  projection_matrix;/* Projection onto critical modes, n_op×n_dim */

    /* Reduced dynamics: dξ/dt = λ ξ + f_nl(ξ) */
    double*  effective_lambda; /* Effective eigenvalues, length n_op */
    void (*reduced_dynamics)(const double* xi, int nop,
                             const void* params, double* dxi);
    void*    reduced_params;

    /* Ginzburg-Landau parameters */
    double   gl_potential_coeff;  /* Coefficient A in GL free energy */
    double   gl_cubic_coeff;      /* Coefficient B (saturation) */
    double   gl_correlation_len;  /* Correlation length ξ ~ (α_c-α)^{-1/2} */
    double   gl_relaxation_time;  /* Critical slowing down τ ~ (α_c-α)^{-1} */
} Haken_OrderParam;

/* -----------------------------------------------------------------
 * Slaving relation: q_fast = h(q_slow)
 * ----------------------------------------------------------------- */

typedef struct {
    int      n_fast;           /* Number of fast (slaved) variables */
    int      n_slow;           /* Number of slow (order param) variables */
    double*  slaving_coeff;    /* Coefficients of slaving expansion */
    double   accuracy;         /* Residual of adiabatic approximation */
    int      expansion_order;  /* Order of Taylor expansion in q_s */
    bool     is_valid;         /* Whether slaving relation converges */

    /* Functional form: q_fast = h(q_slow) */
    void (*slaving_func)(const double* q_slow, int ns,
                         double* q_fast, int nf,
                         const void* coeff, int order);
    void*    coeff_data;
} Haken_SlavingRelation;

/* -----------------------------------------------------------------
 * Complete analysis result
 * ----------------------------------------------------------------- */

typedef struct {
    Haken_System         sys;          /* System under analysis */
    Haken_OrderParam     order_param;  /* Identified order parameters */
    Haken_SlavingRelation slaving;     /* Slaving relation q_f= h(q_s) */
    Haken_BifurcationType bif_type;    /* Type of bifurcation */

    /* Error estimates */
    double   adiabatic_error;   /* Max error from adiabatic elimination */
    double   reduction_error;   /* Error from dimensional reduction */
    double   truncation_error;  /* Error from expansion truncation */

    /* Validation */
    bool     slaving_holds;     /* Whether slaving principle applies */
    double   spectral_gap;      /* Gap between slow and fast eigenvalues */
    double   critical_distance; /* |α - α_c| */
} Haken_SlavingResult;

/* -----------------------------------------------------------------
 * L2: Control parameter sweep descriptor
 * ----------------------------------------------------------------- */

typedef struct {
    double   alpha_start;
    double   alpha_end;
    double   alpha_step;
    int      n_points;
    double*  eigenvalues;     /* n_dim × n_points, eigenvalue traces */
    double*  order_param_amp; /* n_critical × n_points */
    double*  relaxation_times;/* n_stable × n_points */
} Haken_ParameterSweep;

/* ==============================================================
 * Core API — system lifecycle
 * ============================================================== */

/** Allocate and initialize a Haken system of dimension n.
 *  @param n_dim  Total system dimension N
 *  @param L_func Optional: function L(α) → matrix
 *  @param N_func Optional: nonlinear coupling N(q)
 *  @return       Initialized system (caller must free via haken_system_free) */
Haken_System* haken_system_create(int n_dim,
    void (*L_func)(double alpha, int n, double* L_out),
    void (*N_func)(const double* q, int n, double* Nq));

/** Deep-copy a Haken system */
Haken_System* haken_system_clone(const Haken_System* src);

/** Free all resources associated with a system */
void haken_system_free(Haken_System* sys);

/** Set control parameter α */
void haken_set_control_param(Haken_System* sys, double alpha);

/** Set critical parameter α_c */
void haken_set_critical_param(Haken_System* sys, double alpha_c);

/** Set the state vector from external array */
void haken_set_state(Haken_System* sys, const double* q);

/** Update the linear matrix L(α) for the current α */
void haken_update_linear_matrix(Haken_System* sys);

/** Compute eigenvalues and eigenvectors of L(α).
 *  Uses QR algorithm with Wilkinson shifts for real matrices.
 *  O(n³) time, O(n²) space. */
int haken_compute_eigenmodes(Haken_System* sys);

/** Classify modes into critical / stable / unstable based on eigenvalues.
 *  Critical: |Re(λ)| < ε_crit (near zero)
 *  Stable:   Re(λ) < -ε_stable
 *  Unstable: Re(λ) > ε_crit */
void haken_classify_modes(Haken_System* sys, double eps_crit, double eps_stable);

/** Project state onto critical subspace.
 *  q_s = P_critical · q */
void haken_project_critical(const Haken_System* sys, const double* q,
                             double* q_s);

/** Project state onto stable (fast) subspace.
 *  q_f = P_stable · q */
void haken_project_stable(const Haken_System* sys, const double* q,
                           double* q_f);

/** Compute the spectral gap of the system:
 *  γ = max(Re λ_crit) - max(Re λ_stable).
 *  If no gap exists, returns 0. */
double haken_spectral_gap(const Haken_System* sys);

#endif /* HAKEN_CORE_H */
