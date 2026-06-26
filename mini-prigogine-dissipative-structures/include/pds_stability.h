#ifndef PDS_STABILITY_H
#define PDS_STABILITY_H

#include "pds_core.h"

/* ============================================================================
 * PDS Stability Analysis - Linear and Nonlinear Stability
 * ============================================================================ */

typedef struct {
    bool is_stable;
    bool is_marginal;
    bool is_oscillatory;
    int n_unstable_modes;
    double* eigenvalues;
    int n_eigenvalues;
    double* eigenvectors;
    double critical_parameter;
    double oscillation_frequency;
} PDSStabilityResult;

/* Jacobian computation via finite differences */
void pds_compute_jacobian_fd(const PDSDissipativeSystem* sys,
                              void (*rhs_fn)(const double*, const double*, int, double*),
                              const double* params, int n_params,
                              PDSMatrix** J_out);

/* Eigenvalue computation for 2x2 and 3x3 systems */
void pds_eigenvalues_2x2(double a, double b, double c, double d,
                          double* lambda1_real, double* lambda1_imag,
                          double* lambda2_real, double* lambda2_imag);
void pds_eigenvalues_3x3(const double A[3][3],
                          double lambda_real[3], double lambda_imag[3]);

/* Linear stability analysis */
PDSStabilityResult* pds_linear_stability_analysis(const PDSMatrix* J);
void pds_stability_result_free(PDSStabilityResult* sr);

/* Hopf bifurcation detection */
bool pds_detect_hopf_bifurcation(const PDSMatrix* J, double* critical_frequency);

/* Lyapunov exponent computation from trajectory */
double pds_max_lyapunov_exponent(const double* const* trajectory,
                                  int n_points, int dim, double dt,
                                  int embedding_dim, int lag);

/* Entropy-based stability: d/dt (delta^2 S) criterion */
double pds_excess_entropy_production(const PDSDissipativeSystem* sys,
                                      const PDSVector* perturbation);

/* Print stability analysis results */
void pds_stability_print(const PDSStabilityResult* sr);

#endif /* PDS_STABILITY_H */