#ifndef SOFE_DYNAMICS_H
#define SOFE_DYNAMICS_H

#include "sofe_core.h"

typedef struct {
    double* state;
    double* derivative;
    double* previous_state;
    double* error_estimate;
    int n_variables;
    double time;
    double dt;
    int step_count;
} SOFEDynamicalState;

typedef struct {
    SOFEReactionNetwork* network;
    SOFEGrid* grid;
    SOFEField* field;
    double* diffusion_coefficients;
    bool has_diffusion;
    bool has_advection;
    double* advection_velocity;
    double* boundary_conditions;
    int bc_type;
} SOFEReactionDiffusionSystem;

typedef enum {
    SOFE_INT_EULER_EXPLICIT   = 0,
    SOFE_INT_EULER_IMPLICIT   = 1,
    SOFE_INT_RK4              = 2,
    SOFE_INT_RK45_ADAPTIVE    = 3,
    SOFE_INT_SEMI_IMPLICIT    = 4,
    SOFE_INT_CRANK_NICOLSON   = 5,
    SOFE_INT_IMEX             = 6
} SOFEIntegrationMethod;

typedef struct {
    double* eigenvalues_real;
    double* eigenvalues_imag;
    int n_eigenvalues;
    double max_real_part;
    double dominant_frequency;
    SOFEStabilityClass classification;
    double* eigenvectors_right;
    double* eigenvectors_left;
    bool is_hyperbolic;
    double condition_number;
} SOFEStabilityResult;

typedef struct {
    double parameter;
    double parameter_step;
    double* solution;
    int n_variables;
    SOFEBifurcationType detected_type;
    double bifurcation_parameter_value;
    double* null_vector;
    bool bifurcation_detected;
    int continuation_steps;
    double* parameter_history;
    double** solution_history;
    int history_length;
    int history_capacity;
} SOFEContinuation;

typedef struct {
    double* exponents;
    int n_exponents;
    double max_exponent;
    double sum_positive;
    double lyapunov_dimension;
    bool is_chaotic;
    int orthonormalization_interval;
} SOFELyapunovSpectrum;

SOFEDynamicalState* sofe_dyn_state_create(int n_variables);
void sofe_dyn_state_free(SOFEDynamicalState* ds);
void sofe_dyn_state_set(SOFEDynamicalState* ds, double* values);
void sofe_dyn_state_get(SOFEDynamicalState* ds, double* out);
void sofe_dyn_state_copy(SOFEDynamicalState* dst, SOFEDynamicalState* src);

SOFEReactionDiffusionSystem* sofe_rd_create(SOFEReactionNetwork* net,
                                              SOFEGrid* grid);
void sofe_rd_free(SOFEReactionDiffusionSystem* rd);
void sofe_rd_set_diffusion(SOFEReactionDiffusionSystem* rd,
                            double* diffusivities);
void sofe_rd_set_boundary_conditions(SOFEReactionDiffusionSystem* rd,
                                      int bc_type, double* values);
void sofe_rd_set_advection(SOFEReactionDiffusionSystem* rd,
                            double* velocity);

typedef void (*SOFERateFunction)(double* y, double t, double* params,
                                  int n_params, double* dydt);

void sofe_integrate_step(double* y, int n, double t, double dt,
                          SOFERateFunction f, double* params, int n_params,
                          SOFEIntegrationMethod method);
void sofe_integrate_adaptive(double* y, int n, double* t, double* dt,
                              SOFERateFunction f, double* params, int n_params,
                              double tolerance);
void sofe_integrate(double* y, int n, double t0, double t_end, double dt,
                     SOFERateFunction f, double* params, int n_params,
                     SOFEIntegrationMethod method,
                     double** trajectory, int* n_steps);

void sofe_rd_integrate_step(SOFEReactionDiffusionSystem* rd, double dt,
                              SOFEIntegrationMethod method);
void sofe_rd_integrate(SOFEReactionDiffusionSystem* rd, double t_end,
                        double dt, SOFEIntegrationMethod method);

void sofe_network_kinetics(SOFEReactionNetwork* net, double* dcdt);

SOFEStabilityResult* sofe_stability_analyze(double* jacobian, int n,
                                              double epsilon);
void sofe_stability_result_free(SOFEStabilityResult* sr);
void sofe_stability_result_print(SOFEStabilityResult* sr);
void sofe_stability_eigenvalues(double* A, int n,
                                 double* real_parts, double* imag_parts);

SOFEContinuation* sofe_continuation_create(int n_variables);
void sofe_continuation_free(SOFEContinuation* cont);
void sofe_continuation_set_initial(SOFEContinuation* cont,
                                    double parameter, double* solution);
void sofe_continuation_step(SOFEContinuation* cont,
                             SOFERateFunction f, double* params, int n_params,
                             int param_index);
bool sofe_continuation_detect_bifurcation(SOFEContinuation* cont,
                                           double* eigenvalues);

SOFELyapunovSpectrum* sofe_lyapunov_compute(
    double* y0, int n, double t_max, double dt,
    SOFERateFunction f, double* params, int n_params,
    int n_exponents);
void sofe_lyapunov_free(SOFELyapunovSpectrum* ls);
void sofe_lyapunov_print(SOFELyapunovSpectrum* ls);

void sofe_brusselator_rates(double* y, double t, double* p, int np, double* dydt);
void sofe_schlogl_rates(double* y, double t, double* p, int np, double* dydt);
void sofe_lotka_volterra_rates(double* y, double t, double* p, int np, double* dydt);
void sofe_oregonator_rates(double* y, double t, double* p, int np, double* dydt);
void sofe_ginzburg_landau_rates(double* y, double t, double* p, int np, double* dydt);

double sofe_rk4_step_component(SOFERateFunction f, double* y, double t,
                                double dt, double* params, int n_params,
                                int component_idx);
void sofe_finite_difference_jacobian(SOFERateFunction f, double* y,
                                      double t, double* params, int n_params,
                                      int n, double* jacobian, double eps);

#endif
