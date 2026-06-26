#ifndef HAKEN_MODE_DECOMP_H
#define HAKEN_MODE_DECOMP_H
#include "haken_core.h"

/* ==============================================================
 * haken_mode_decomp.h — Mode Decomposition & Spectral Analysis
 *
 * The separation of timescales is fundamental to the slaving
 * principle. Modes are decomposed by their eigenvalues:
 *
 *   λ_k:  Re(λ_k) ~ O(1)     → stable/fast (slaved)
 *         Re(λ_k) ~ O(ε)     → critical/slow (order parameters)
 *         Re(λ_k) ~ O(ε²)    → neutral (higher-order critical)
 *
 * Spectral gap condition (Haken 1977, §7.2):
 *   γ = max_k{Re(λ_critical)} - max_j{Re(λ_stable)} ≫ 0
 *
 * When the spectral gap exists, adiabatic elimination is justified.
 *
 * References:
 *   Haken (1977) Synergetics, §7.1-7.3
 *   Mori & Kuramoto (1998) Dissipative Structures and Chaos
 * ============================================================== */

/* -----------------------------------------------------------------
 * L2: Spectral decomposition types
 * ----------------------------------------------------------------- */

/** QR algorithm configuration */
typedef struct {
    int     max_iterations;    /* Max QR iterations (default 100) */
    double  tolerance;         /* Convergence tolerance (default 1e-12) */
    bool    compute_eigenvectors; /* Whether to compute eigenvectors */
    bool    sort_by_real_part; /* Sort eigenvalues by real part descending */
} Haken_QRConfig;

/** Hessenberg decomposition intermediate storage */
typedef struct {
    double* H;              /* Upper Hessenberg matrix, n×n */
    double* Q;              /* Orthogonal transformation, n×n */
    int     n;
} Haken_Hessenberg;

/** Real Schur form: block upper triangular */
typedef struct {
    double* T;              /* Quasi-triangular Schur form, n×n */
    double* U;              /* Orthogonal Schur vectors, n×n */
    int     n;
} Haken_SchurForm;

/** Eigenpair: eigenvalue + eigenvector */
typedef struct {
    double   lambda_re;
    double   lambda_im;
    double*  vector_re;
    double*  vector_im;
    double   residual;      /* |A v - λ v| */
} Haken_Eigenpair;

/* -----------------------------------------------------------------
 * L3: Mode decomposition results
 * ----------------------------------------------------------------- */

/** Full modal decomposition */
typedef struct {
    int      n_modes;
    Haken_Mode* modes;       /* Array of modes */
    double*  eigenvalues_re;  /* Real parts, length n_modes */
    double*  eigenvalues_im;  /* Imag parts, length n_modes */
    double   spectral_radius; /* max |λ_k| */
    double   spectral_abscissa;/* max Re(λ_k) */
    double   spectral_gap;    /* Re(λ_slow) - Re(λ_fast) */
    bool     has_gap;         /* Gap sufficient for slaving */
} Haken_ModeDecomp;

/* -----------------------------------------------------------------
 * Matrix operations for linear stability
 * ----------------------------------------------------------------- */

/** Balance a matrix to improve eigenvalue computation accuracy.
 *  Uses Osborne's method. Returns the scaling diagonal D.
 *  Reference: Osborne (1960) "On Preconditioning Matrices",
 *             J. ACM 7, 338-345. */
void haken_balance_matrix(const double* A, int n, double* balanced, double* D);

/** Reduce a square matrix to upper Hessenberg form:
 *  H = Q^T A Q, where Q is orthogonal.
 *  Uses Householder reflections. O(n³) time.
 *  Reference: Golub & Van Loan (2013) Matrix Computations, §7.4 */
void haken_hessenberg(const double* A, int n,
                      double* H, double* Q);

/** QR decomposition of a matrix A = Q·R.
 *  Uses Householder reflections.
 *  Reference: Golub & Van Loan (2013) §5.2 */
void haken_qr_decompose(const double* A, int rows, int cols,
                         double* Q, double* R);

/** Multiply two matrices C = A·B (all square n×n). */
void haken_matrix_multiply(const double* A, const double* B,
                            int n, double* C);

/** Matrix transpose of n×n matrix. */
void haken_matrix_transpose(const double* A, int n, double* At);

/** Compute eigenvalues of a real matrix using QR algorithm
 *  with Wilkinson shifts. Returns eigenvalues in order of
 *  decreasing real part.
 *  O(n³) per iteration, typically ~2n iterations.
 *  Reference: Wilkinson (1965) The Algebraic Eigenvalue Problem */
int haken_qr_eigenvalues(double* A, int n,
                          double* eigenvalues_re,
                          double* eigenvalues_im,
                          const Haken_QRConfig* config);

/** Compute all eigenpairs from the Schur decomposition.
 *  For each real eigenvalue: one eigenpair.
 *  For each complex conjugate pair: two eigenpairs.
 *  Reference: Golub & Van Loan (2013) §7.6 */
int haken_compute_eigenpairs(const double* A, int n,
                              Haken_Eigenpair* pairs,
                              const Haken_QRConfig* config);

/** Inverse iteration to refine an approximate eigenvector.
 *  Solves (A - μI) v_{k+1} = v_k.
 *  Reference: Wilkinson (1965) §9 */
void haken_inverse_iteration(const double* A, int n,
                              double lambda_re, double lambda_im,
                              double* eigenvector_re,
                              double* eigenvector_im,
                              int max_iter, double tol);

/** Compute left (adjoint) eigenvectors.
 *  Left eigenvectors satisfy: w^T A = λ w^T.
 *  Equivalent to right eigenvectors of A^T. */
void haken_left_eigenvectors(const double* A, int n,
                              const Haken_Eigenpair* right_pairs,
                              Haken_Eigenpair* left_pairs);

/** Compute the resolvent (A - λI)^{-1} for stability analysis.
 *  Uses LU decomposition. */
double haken_resolvent_norm(const double* A, int n,
                             double lambda_re, double lambda_im);

/** Compute pseudospectrum: the set of λ such that
 *  ||(A-λI)^{-1}|| ≥ ε^{-1}.
 *  Useful when A is non-normal.
 *  Reference: Trefethen & Embree (2005) Spectra and Pseudospectra */
int haken_pseudospectrum(const double* A, int n, double epsilon,
                          int grid_size, double range,
                          double* pseudo_grid);

/* -----------------------------------------------------------------
 * Spectral gap analysis
 * ----------------------------------------------------------------- */

/** Compute the spectral gap: difference between the smallest
 *  (in magnitude) critical eigenvalue and largest stable eigenvalue.
 *  Returns 0 if no gap exists. */
double haken_spectral_gap_decomp(const Haken_ModeDecomp* decomp,
                                  double eps_critical);

/** Check if the timescale separation condition holds:
 *  τ_fast ≪ τ_slow, i.e., spectral_gap ≫ ε.
 *  This justifies adiabatic elimination. */
bool haken_timescale_separation(const Haken_ModeDecomp* decomp,
                                 double min_gap_ratio);

/** Compute participation factors: how much each state variable
 *  participates in each mode.
 *  p_{ki} = u_{ki} · v_{ik} where u is left eigenvector,
 *  v is right eigenvector. (Normalized such that Σ_k p_{ki}=1) */
void haken_participation_factors(const Haken_Eigenpair* right,
                                  const Haken_Eigenpair* left,
                                  int n, double* p_matrix);

/** Sort modes by relaxation time (descending: slowest first) */
void haken_sort_modes_by_timescale(Haken_Mode* modes, int n);

/** Find the mode closest to marginal stability:
 *  argmin_k |Re(λ_k)| */
int haken_find_most_critical_mode(const Haken_Mode* modes, int n);

/** Compute the condition number of the eigenvector matrix:
 *  κ(V) = ||V||·||V^{-1}||.
 *  Large κ indicates nearly defective matrix — slaving may fail.
 *  Reference: Trefethen (1997) "Pseudospectra of Linear Operators" */
double haken_eigenvector_condition(const double* eigenvectors,
                                    int n);

#endif /* HAKEN_MODE_DECOMP_H */
