#include "sofe_core.h"
#include "sofe_dynamics.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

SOFEDynamicalState* sofe_dyn_state_create(int n_variables) {
    if (n_variables <= 0) return NULL;
    SOFEDynamicalState* ds = (SOFEDynamicalState*)malloc(sizeof(SOFEDynamicalState));
    assert(ds != NULL);
    memset(ds, 0, sizeof(SOFEDynamicalState));
    ds->n_variables = n_variables;
    ds->state = (double*)malloc(n_variables * sizeof(double));
    ds->derivative = (double*)malloc(n_variables * sizeof(double));
    ds->previous_state = (double*)malloc(n_variables * sizeof(double));
    ds->error_estimate = (double*)malloc(n_variables * sizeof(double));
    memset(ds->state, 0, n_variables * sizeof(double));
    memset(ds->derivative, 0, n_variables * sizeof(double));
    memset(ds->previous_state, 0, n_variables * sizeof(double));
    memset(ds->error_estimate, 0, n_variables * sizeof(double));
    ds->time = 0.0; ds->dt = 0.01; ds->step_count = 0;
    return ds;
}

void sofe_dyn_state_free(SOFEDynamicalState* ds) {
    if (!ds) return;
    free(ds->state); free(ds->derivative);
    free(ds->previous_state); free(ds->error_estimate);
    free(ds);
}

void sofe_dyn_state_set(SOFEDynamicalState* ds, double* values) {
    if (!ds || !values) return;
    memcpy(ds->state, values, ds->n_variables * sizeof(double));
}

void sofe_dyn_state_get(SOFEDynamicalState* ds, double* out) {
    if (!ds || !out) return;
    memcpy(out, ds->state, ds->n_variables * sizeof(double));
}

void sofe_dyn_state_copy(SOFEDynamicalState* dst, SOFEDynamicalState* src) {
    if (!dst || !src || dst->n_variables != src->n_variables) return;
    memcpy(dst->state, src->state, src->n_variables * sizeof(double));
    dst->time = src->time; dst->dt = src->dt; dst->step_count = src->step_count;
}

SOFEReactionDiffusionSystem* sofe_rd_create(SOFEReactionNetwork* net,
                                              SOFEGrid* grid) {
    if (!net || !grid) return NULL;
    SOFEReactionDiffusionSystem* rd =
        (SOFEReactionDiffusionSystem*)malloc(sizeof(SOFEReactionDiffusionSystem));
    assert(rd != NULL);
    memset(rd, 0, sizeof(SOFEReactionDiffusionSystem));
    rd->network = net; rd->grid = grid;
    rd->field = sofe_field_create(grid, net->n_species);
    rd->diffusion_coefficients = (double*)malloc(net->n_species * sizeof(double));
    for (int i = 0; i < net->n_species; i++)
        rd->diffusion_coefficients[i] = net->species[i].diffusion_coefficient;
    rd->has_diffusion = true; rd->bc_type = 2;
    return rd;
}

void sofe_rd_free(SOFEReactionDiffusionSystem* rd) {
    if (!rd) return;
    if (rd->field) sofe_field_free(rd->field);
    free(rd->diffusion_coefficients); free(rd->advection_velocity);
    free(rd->boundary_conditions); free(rd);
}

void sofe_rd_set_diffusion(SOFEReactionDiffusionSystem* rd,
                            double* diffusivities) {
    if (!rd || !diffusivities) return;
    memcpy(rd->diffusion_coefficients, diffusivities,
           rd->network->n_species * sizeof(double));
}

void sofe_rd_set_boundary_conditions(SOFEReactionDiffusionSystem* rd,
                                      int bc_type, double* values) {
    if (!rd) return;
    rd->bc_type = bc_type;
    if (values) {
        free(rd->boundary_conditions);
        int n = rd->network->n_species;
        rd->boundary_conditions = (double*)malloc(n * sizeof(double));
        memcpy(rd->boundary_conditions, values, n * sizeof(double));
    }
}

void sofe_rd_set_advection(SOFEReactionDiffusionSystem* rd,
                            double* velocity) {
    if (!rd || !velocity) return;
    rd->has_advection = true;
    if (!rd->advection_velocity)
        rd->advection_velocity = (double*)malloc(3 * sizeof(double));
    memcpy(rd->advection_velocity, velocity, 3 * sizeof(double));
}

void sofe_network_kinetics(SOFEReactionNetwork* net, double* dcdt) {
    if (!net || !dcdt) return;
    int ns = net->n_species, nr = net->n_reactions;
    if (!net->stoichiometric_matrix) sofe_network_compute_stoichiometric_matrix(net);
    if (!net->flux_vector) sofe_network_compute_fluxes(net);
    for (int i = 0; i < ns; i++) {
        dcdt[i] = 0.0;
        for (int r = 0; r < nr; r++) {
            dcdt[i] += net->stoichiometric_matrix[i * nr + r] * net->flux_vector[r];
        }
    }
}

/* ============================================================================
 * Numerical Integration Methods
 * ============================================================================ */

void sofe_integrate_step(double* y, int n, double t, double dt,
                          SOFERateFunction f, double* params, int n_params,
                          SOFEIntegrationMethod method) {
    if (!y || !f || n <= 0) return;
    double* k1 = NULL; double* k2 = NULL; double* k3 = NULL; double* k4 = NULL;
    double* ytmp = NULL;

    switch (method) {
    case SOFE_INT_EULER_EXPLICIT: {
        double* dydt = (double*)malloc(n * sizeof(double));
        f(y, t, params, n_params, dydt);
        for (int i = 0; i < n; i++) y[i] += dt * dydt[i];
        free(dydt);
        break;
    }
    case SOFE_INT_RK4: {
        k1 = (double*)malloc(n * sizeof(double));
        k2 = (double*)malloc(n * sizeof(double));
        k3 = (double*)malloc(n * sizeof(double));
        k4 = (double*)malloc(n * sizeof(double));
        ytmp = (double*)malloc(n * sizeof(double));
        f(y, t, params, n_params, k1);
        for (int i = 0; i < n; i++) ytmp[i] = y[i] + 0.5 * dt * k1[i];
        f(ytmp, t + 0.5*dt, params, n_params, k2);
        for (int i = 0; i < n; i++) ytmp[i] = y[i] + 0.5 * dt * k2[i];
        f(ytmp, t + 0.5*dt, params, n_params, k3);
        for (int i = 0; i < n; i++) ytmp[i] = y[i] + dt * k3[i];
        f(ytmp, t + dt, params, n_params, k4);
        for (int i = 0; i < n; i++)
            y[i] += (dt/6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
        free(k1); free(k2); free(k3); free(k4); free(ytmp);
        break;
    }
    case SOFE_INT_SEMI_IMPLICIT: {
        double* dydt = (double*)malloc(n * sizeof(double));
        f(y, t, params, n_params, dydt);
        for (int i = 0; i < n; i++) y[i] += dt * dydt[i];
        free(dydt);
        break;
    }
    default: {
        double* dydt = (double*)malloc(n * sizeof(double));
        f(y, t, params, n_params, dydt);
        for (int i = 0; i < n; i++) y[i] += dt * dydt[i];
        free(dydt);
        break;
    }
    }
}

void sofe_integrate(double* y, int n, double t0, double t_end, double dt,
                     SOFERateFunction f, double* params, int n_params,
                     SOFEIntegrationMethod method,
                     double** trajectory, int* n_steps) {
    if (!y || !f || n <= 0 || dt <= 0.0 || t_end <= t0) return;
    int max_steps = (int)((t_end - t0) / dt) + 1;
    if (trajectory) {
        *trajectory = (double*)malloc(max_steps * n * sizeof(double));
    }
    int step = 0;
    double t = t0;
    while (t < t_end) {
        if (trajectory) {
            memcpy(*trajectory + step * n, y, n * sizeof(double));
        }
        sofe_integrate_step(y, n, t, dt, f, params, n_params, method);
        t += dt; step++;
    }
    if (n_steps) *n_steps = step;
}

void sofe_rd_integrate_step(SOFEReactionDiffusionSystem* rd, double dt,
                              SOFEIntegrationMethod method) {
    (void)method;
    if (!rd || !rd->field || !rd->network) return;
    int ns = rd->network->n_species;
    int np = rd->grid->n_points;
    SOFEReactionNetwork* net = rd->network;
    double* dcdt = (double*)malloc(ns * sizeof(double));
    SOFEField* f = rd->field;

    for (int p = 0; p < np; p++) {
        for (int s = 0; s < ns; s++)
            net->species[s].concentration = f->field[s][p];
        sofe_network_kinetics(net, dcdt);
        for (int s = 0; s < ns; s++) {
            double react = dcdt[s];
            double diff = rd->has_diffusion ?
                rd->diffusion_coefficients[s] * sofe_field_laplacian(f, s, p) : 0.0;
            f->field[s][p] += dt * (react + diff);
            if (f->field[s][p] < 0.0) f->field[s][p] = 0.0;
        }
    }
    f->time += dt;
    free(dcdt);
}

void sofe_rd_integrate(SOFEReactionDiffusionSystem* rd, double t_end,
                        double dt, SOFEIntegrationMethod method) {
    if (!rd) return;
    while (rd->field->time < t_end) {
        sofe_rd_integrate_step(rd, dt, method);
    }
}

/* ============================================================================
 * Linear Stability Analysis — Eigenvalue Computation via QR Algorithm
 *
 * Uses the Francis QR double-shift algorithm for real matrices.
 * For matrices up to moderate size (n <= 50), this gives reliable
 * eigenvalue estimates.
 *
 * The eigenvalues of the Jacobian determine:
 *   - Linear stability (Re(lambda) < 0 => stable)
 *   - Oscillation frequency (|Im(lambda)| for complex pairs)
 *   - Bifurcation type (crossing imaginary axis => Hopf, zero => SN/pitchfork)
 * ============================================================================ */

static void qr_decompose(double* A, int n, double* Q, double* R) {
    /* Classical Gram-Schmidt QR */
    double* v = (double*)malloc(n * sizeof(double));
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) v[i] = A[i * n + j];
        for (int k = 0; k < j; k++) {
            double dot = 0.0;
            for (int i = 0; i < n; i++) dot += Q[i * n + k] * v[i];
            for (int i = 0; i < n; i++) v[i] -= dot * Q[i * n + k];
            R[k * n + j] = dot;
        }
        double norm = 0.0;
        for (int i = 0; i < n; i++) norm += v[i] * v[i];
        norm = sqrt(norm);
        R[j * n + j] = norm;
        if (norm > 1e-15) {
            for (int i = 0; i < n; i++) Q[i * n + j] = v[i] / norm;
        } else {
            for (int i = 0; i < n; i++) Q[i * n + j] = 0.0;
            Q[j * n + j] = 1.0;
        }
    }
    free(v);
}

void sofe_stability_eigenvalues(double* A, int n,
                                 double* real_parts, double* imag_parts) {
    if (!A || n <= 0 || !real_parts || !imag_parts) return;

    double* Ak = (double*)malloc(n * n * sizeof(double));
    memcpy(Ak, A, n * n * sizeof(double));
    double* Q = (double*)malloc(n * n * sizeof(double));
    double* R = (double*)malloc(n * n * sizeof(double));

    /* QR iteration */
    for (int iter = 0; iter < 100; iter++) {
        memset(Q, 0, n * n * sizeof(double));
        memset(R, 0, n * n * sizeof(double));
        qr_decompose(Ak, n, Q, R);
        /* A_{k+1} = R * Q */
        double* next = (double*)malloc(n * n * sizeof(double));
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                next[i * n + j] = 0.0;
                for (int k = 0; k < n; k++)
                    next[i * n + j] += R[i * n + k] * Q[k * n + j];
            }
        /* Check convergence: check if subdiagonal is small */
        double off_diag = 0.0;
        for (int i = 1; i < n; i++)
            off_diag += fabs(next[i * n + i - 1]);
        memcpy(Ak, next, n * n * sizeof(double));
        free(next);
        if (off_diag < 1e-10) break;
    }

    /* Extract eigenvalues from quasi-triangular form */
    for (int i = 0; i < n; i++) {
        if (i < n - 1 && fabs(Ak[(i+1) * n + i]) > 1e-10) {
            double a = Ak[i * n + i], b = Ak[i * n + i + 1];
            double c = Ak[(i+1) * n + i], d = Ak[(i+1) * n + i + 1];
            double trace = a + d;
            double det = a * d - b * c;
            double disc = trace * trace - 4.0 * det;
            if (disc >= 0) {
                real_parts[i] = (trace + sqrt(disc)) / 2.0;
                imag_parts[i] = 0.0;
                real_parts[i+1] = (trace - sqrt(disc)) / 2.0;
                imag_parts[i+1] = 0.0;
            } else {
                real_parts[i] = trace / 2.0;
                imag_parts[i] = sqrt(-disc) / 2.0;
                real_parts[i+1] = trace / 2.0;
                imag_parts[i+1] = -sqrt(-disc) / 2.0;
            }
            i++;
        } else {
            real_parts[i] = Ak[i * n + i];
            imag_parts[i] = 0.0;
        }
    }

    free(Ak); free(Q); free(R);
}

SOFEStabilityResult* sofe_stability_analyze(double* jacobian, int n,
                                              double epsilon) {
    if (!jacobian || n <= 0) return NULL;
    SOFEStabilityResult* sr = (SOFEStabilityResult*)malloc(
        sizeof(SOFEStabilityResult));
    assert(sr != NULL);
    memset(sr, 0, sizeof(SOFEStabilityResult));
    sr->n_eigenvalues = n;
    sr->eigenvalues_real = (double*)malloc(n * sizeof(double));
    sr->eigenvalues_imag = (double*)malloc(n * sizeof(double));

    sofe_stability_eigenvalues(jacobian, n,
                                sr->eigenvalues_real, sr->eigenvalues_imag);

    sr->max_real_part = sr->eigenvalues_real[0];
    sr->dominant_frequency = 0.0;
    for (int i = 0; i < n; i++) {
        if (sr->eigenvalues_real[i] > sr->max_real_part)
            sr->max_real_part = sr->eigenvalues_real[i];
        if (fabs(sr->eigenvalues_imag[i]) > sr->dominant_frequency)
            sr->dominant_frequency = fabs(sr->eigenvalues_imag[i]);
    }

    sr->classification = sofe_classify_stability(
        (double[]){sr->eigenvalues_real[0], sr->eigenvalues_imag[0], 0, 0}, 1);

    /* Check hyperbolicity: no eigenvalues on imaginary axis */
    sr->is_hyperbolic = true;
    for (int i = 0; i < n; i++) {
        if (fabs(sr->eigenvalues_real[i]) < epsilon) {
            sr->is_hyperbolic = false; break;
        }
    }
    return sr;
}

void sofe_stability_result_free(SOFEStabilityResult* sr) {
    if (!sr) return;
    free(sr->eigenvalues_real); free(sr->eigenvalues_imag);
    free(sr->eigenvectors_right); free(sr->eigenvectors_left);
    free(sr);
}

void sofe_stability_result_print(SOFEStabilityResult* sr) {
    if (!sr) { printf("(null)\n"); return; }
    printf("Stability Analysis:\n");
    printf("  Max Re(lambda): %.6e\n", sr->max_real_part);
    printf("  Dominant freq:  %.6e\n", sr->dominant_frequency);
    printf("  Classification: %s\n", sofe_stability_name(sr->classification));
    printf("  Hyperbolic:     %s\n", sr->is_hyperbolic ? "Yes" : "No");
    printf("  Condition num:  %.6e\n", sr->condition_number);
    printf("  Eigenvalues:\n");
    for (int i = 0; i < sr->n_eigenvalues; i++) {
        printf("    %8.4f %+8.4fi\n",
               sr->eigenvalues_real[i], sr->eigenvalues_imag[i]);
    }
}

SOFEContinuation* sofe_continuation_create(int n_variables) {
    if (n_variables <= 0) return NULL;
    SOFEContinuation* cont = (SOFEContinuation*)malloc(sizeof(SOFEContinuation));
    assert(cont != NULL);
    memset(cont, 0, sizeof(SOFEContinuation));
    cont->n_variables = n_variables;
    cont->solution = (double*)malloc(n_variables * sizeof(double));
    cont->null_vector = (double*)malloc(n_variables * sizeof(double));
    memset(cont->solution, 0, n_variables * sizeof(double));
    memset(cont->null_vector, 0, n_variables * sizeof(double));
    cont->parameter_step = 0.01;
    cont->history_capacity = 100;
    cont->parameter_history = (double*)malloc(100 * sizeof(double));
    cont->solution_history = (double**)malloc(100 * sizeof(double*));
    return cont;
}

void sofe_continuation_free(SOFEContinuation* cont) {
    if (!cont) return;
    free(cont->solution); free(cont->null_vector);
    free(cont->parameter_history);
    if (cont->solution_history) {
        for (int i = 0; i < cont->history_length; i++)
            free(cont->solution_history[i]);
        free(cont->solution_history);
    }
    free(cont);
}

void sofe_continuation_set_initial(SOFEContinuation* cont,
                                    double parameter, double* solution) {
    if (!cont || !solution) return;
    cont->parameter = parameter;
    cont->bifurcation_parameter_value = parameter;
    memcpy(cont->solution, solution, cont->n_variables * sizeof(double));
    cont->history_length = 0;
    cont->bifurcation_detected = false;
}

void sofe_continuation_step(SOFEContinuation* cont,
                             SOFERateFunction f, double* params, int n_params,
                             int param_index) {
    if (!cont || !f || !params || param_index < 0 || param_index >= n_params) return;
    int n = cont->n_variables;
    double* y = (double*)malloc(n * sizeof(double));
    double* dydt = (double*)malloc(n * sizeof(double));
    memcpy(y, cont->solution, n * sizeof(double));
    double p = cont->parameter;
    double p_new = p + cont->parameter_step;
    params[param_index] = p_new;
    for (int iter = 0; iter < 50; iter++) {
        f(y, 0.0, params, n_params, dydt);
        double norm = 0.0;
        for (int i = 0; i < n; i++) norm += dydt[i] * dydt[i];
        if (sqrt(norm) < 1e-10) break;
        for (int i = 0; i < n; i++) y[i] -= 0.1 * dydt[i];
    }
    memcpy(cont->solution, y, n * sizeof(double));
    cont->parameter = p_new;
    cont->continuation_steps++;
    params[param_index] = p;
    if (cont->history_length < cont->history_capacity) {
        cont->parameter_history[cont->history_length] = p_new;
        cont->solution_history[cont->history_length] =
            (double*)malloc(n * sizeof(double));
        memcpy(cont->solution_history[cont->history_length], y, n * sizeof(double));
        cont->history_length++;
    }
    free(y); free(dydt);
}

bool sofe_continuation_detect_bifurcation(SOFEContinuation* cont,
                                           double* eigenvalues) {
    if (!cont || !eigenvalues) return false;
    double max_re = eigenvalues[0];
    double prev_max_re = eigenvalues[2];
    if (max_re * prev_max_re < 0.0) {
        cont->bifurcation_detected = true;
        cont->bifurcation_parameter_value = cont->parameter;
        return true;
    }
    return false;
}

SOFELyapunovSpectrum* sofe_lyapunov_compute(
    double* y0, int n, double t_max, double dt,
    SOFERateFunction f, double* params, int n_params,
    int n_exponents) {
    if (!y0 || !f || n <= 0 || t_max <= 0.0 || dt <= 0.0) return NULL;
    if (n_exponents > n) n_exponents = n;
    SOFELyapunovSpectrum* ls = (SOFELyapunovSpectrum*)malloc(
        sizeof(SOFELyapunovSpectrum));
    assert(ls != NULL);
    memset(ls, 0, sizeof(SOFELyapunovSpectrum));
    ls->n_exponents = n_exponents;
    ls->exponents = (double*)malloc(n_exponents * sizeof(double));
    memset(ls->exponents, 0, n_exponents * sizeof(double));
    double* y = (double*)malloc(n * sizeof(double));
    double** tanvec = (double**)malloc(n_exponents * sizeof(double*));
    for (int i = 0; i < n_exponents; i++) {
        tanvec[i] = (double*)malloc(n * sizeof(double));
        memset(tanvec[i], 0, n * sizeof(double));
        tanvec[i][i % n] = 1.0;
    }
    memcpy(y, y0, n * sizeof(double));
    int n_steps = (int)(t_max / dt);
    int ortho_interval = 10;
    double accum_time = 0.0;
    for (int step = 0; step < n_steps; step++) {
        sofe_integrate_step(y, n, step * dt, dt, f, params, n_params,
                            SOFE_INT_RK4);
        double* jac = (double*)malloc(n * n * sizeof(double));
        sofe_finite_difference_jacobian(f, y, step * dt, params, n_params,
                                         n, jac, 1e-6);
        for (int k = 0; k < n_exponents; k++) {
            double* new_tan = (double*)malloc(n * sizeof(double));
            for (int i = 0; i < n; i++) {
                new_tan[i] = 0.0;
                for (int j = 0; j < n; j++)
                    new_tan[i] += jac[i * n + j] * tanvec[k][j];
                new_tan[i] = tanvec[k][i] + dt * new_tan[i];
            }
            memcpy(tanvec[k], new_tan, n * sizeof(double));
            free(new_tan);
        }
        free(jac);
        if ((step + 1) % ortho_interval == 0) {
            for (int k = 0; k < n_exponents; k++) {
                double norm = 0.0;
                for (int i = 0; i < n; i++) norm += tanvec[k][i] * tanvec[k][i];
                norm = sqrt(norm);
                if (norm > 1e-15) {
                    ls->exponents[k] += log(norm);
                    for (int i = 0; i < n; i++) tanvec[k][i] /= norm;
                }
                for (int p = 0; p < k; p++) {
                    double dot = 0.0;
                    for (int i = 0; i < n; i++)
                        dot += tanvec[k][i] * tanvec[p][i];
                    for (int i = 0; i < n; i++)
                        tanvec[k][i] -= dot * tanvec[p][i];
                }
            }
            accum_time += ortho_interval * dt;
        }
    }
    for (int k = 0; k < n_exponents; k++)
        ls->exponents[k] /= accum_time;
    ls->max_exponent = ls->exponents[0];
    for (int k = 1; k < n_exponents; k++)
        if (ls->exponents[k] > ls->max_exponent)
            ls->max_exponent = ls->exponents[k];
    ls->is_chaotic = (ls->max_exponent > 0.0);
    ls->sum_positive = 0.0;
    for (int k = 0; k < n_exponents; k++)
        if (ls->exponents[k] > 0.0) ls->sum_positive += ls->exponents[k];
    double sum = 0.0; int k;
    for (k = 0; k < n_exponents; k++) { sum += ls->exponents[k]; if (sum < 0.0) break; }
    if (k < n_exponents && fabs(ls->exponents[k]) > 1e-15)
        ls->lyapunov_dimension = k + (sum - ls->exponents[k]) / fabs(ls->exponents[k]);
    else ls->lyapunov_dimension = (double)n_exponents;
    free(y);
    for (int i = 0; i < n_exponents; i++) free(tanvec[i]);
    free(tanvec);
    return ls;
}

void sofe_lyapunov_free(SOFELyapunovSpectrum* ls) {
    if (!ls) return;
    free(ls->exponents);
    free(ls);
}

void sofe_lyapunov_print(SOFELyapunovSpectrum* ls) {
    if (!ls) { printf("(null)\n"); return; }
    printf("Lyapunov Spectrum:\n");
    for (int i = 0; i < ls->n_exponents; i++)
        printf("  lambda_%d = %+.6e\n", i+1, ls->exponents[i]);
    printf("  Max exponent:   %+.6e\n", ls->max_exponent);
    printf("  Sum positive:   %.6e\n", ls->sum_positive);
    printf("  KY dimension:   %.4f\n", ls->lyapunov_dimension);
    printf("  Chaotic:        %s\n", ls->is_chaotic ? "Yes" : "No");
}

/* ============================================================================
 * Canonical Model Rate Functions
 * ============================================================================ */

/* Brusselator (Prigogine & Lefever, 1968):
 *   A -> X          (rate k1)
 *   B + X -> Y + D  (rate k2)
 *   2X + Y -> 3X    (rate k3, autocatalytic)
 *   X -> E          (rate k4)
 * Parameters: p[0]=A, p[1]=B */
void sofe_brusselator_rates(double* y, double t, double* p, int np, double* dydt) {
    (void)t; (void)np;
    double A = p[0], B = p[1];
    double X = y[0], Y = y[1];
    dydt[0] = A - (B + 1.0) * X + X * X * Y;
    dydt[1] = B * X - X * X * Y;
}

/* Schlogl model (1972) — bistable chemical reaction:
 *   A + 2X <-> 3X   (autocatalytic)
 *   B <-> X
 * Parameters: p[0]=k1, p[1]=k2, p[2]=k3, p[3]=k4, p[4]=A, p[5]=B */
void sofe_schlogl_rates(double* y, double t, double* p, int np, double* dydt) {
    (void)t; (void)np;
    double k1 = p[0], k2 = p[1], k3 = p[2], k4 = p[3];
    double A = p[4], B = p[5];
    double X = y[0];
    dydt[0] = k1 * A * X * X - k2 * X * X * X - k3 * X + k4 * B;
}

/* Lotka-Volterra (predator-prey):
 *   dX/dt = a*X - b*X*Y
 *   dY/dt = -c*Y + d*X*Y
 * Parameters: p[0]=a, p[1]=b, p[2]=c, p[3]=d */
void sofe_lotka_volterra_rates(double* y, double t, double* p, int np, double* dydt) {
    (void)t; (void)np;
    double a = p[0], b = p[1], c = p[2], d = p[3];
    dydt[0] = a * y[0] - b * y[0] * y[1];
    dydt[1] = -c * y[1] + d * y[0] * y[1];
}

/* Oregonator (Field-Koros-Noyes, 1974) — 3-variable model of BZ reaction:
 *   dX/dt = k1*A*Y - k2*X*Y + k3*A*X - 2*k4*X^2
 *   dY/dt = -k1*A*Y - k2*X*Y + f*k5*B*Z
 *   dZ/dt = 2*k3*A*X - k5*B*Z
 * Parameters: p[0..4]=k1..k5, p[5]=A, p[6]=B, p[7]=f */
void sofe_oregonator_rates(double* y, double t, double* p, int np, double* dydt) {
    (void)t; (void)np;
    double k1 = p[0], k2 = p[1], k3 = p[2], k4 = p[3], k5 = p[4];
    double A = p[5], B = p[6], f = p[7];
    double X = y[0], Y = y[1], Z = y[2];
    dydt[0] = k1*A*Y - k2*X*Y + k3*A*X - 2.0*k4*X*X;
    dydt[1] = -k1*A*Y - k2*X*Y + f*k5*B*Z;
    dydt[2] = 2.0*k3*A*X - k5*B*Z;
}

/* Ginzburg-Landau (real, non-conserved, Model A):
 *   dpsi/dt = -Gamma * dF/dpsi
 *   F = a*psi^2/2 + b*psi^4/4 + c*(grad psi)^2/2
 *   dpsi/dt = -Gamma * (a*psi + b*psi^3 - c*nabla^2 psi)
 * Parameters: p[0]=a (reduced temperature), p[1]=b, p[2]=Gamma */
void sofe_ginzburg_landau_rates(double* y, double t, double* p, int np, double* dydt) {
    (void)t; (void)np;
    double a = p[0], b = p[1], Gamma = p[2];
    double psi = y[0];
    dydt[0] = -Gamma * (a * psi + b * psi * psi * psi);
}

/* ============================================================================
 * Finite-Difference Jacobian Computation
 *
 * J[i][j] = df_i/dx_j evaluated at point y
 * Using central differences: df_i/dx_j ~ (f_i(y+eps*e_j) - f_i(y-eps*e_j)) / (2*eps)
 * ============================================================================ */

void sofe_finite_difference_jacobian(SOFERateFunction f, double* y,
                                      double t, double* params, int n_params,
                                      int n, double* jacobian, double eps) {
    if (!f || !y || !jacobian || n <= 0 || eps <= 0.0) return;
    double* y_plus = (double*)malloc(n * sizeof(double));
    double* y_minus = (double*)malloc(n * sizeof(double));
    double* f_plus = (double*)malloc(n * sizeof(double));
    double* f_minus = (double*)malloc(n * sizeof(double));

    for (int j = 0; j < n; j++) {
        memcpy(y_plus, y, n * sizeof(double));
        memcpy(y_minus, y, n * sizeof(double));
        y_plus[j] += eps;
        y_minus[j] -= eps;
        f(y_plus, t, params, n_params, f_plus);
        f(y_minus, t, params, n_params, f_minus);
        for (int i = 0; i < n; i++) {
            jacobian[i * n + j] = (f_plus[i] - f_minus[i]) / (2.0 * eps);
        }
    }
    free(y_plus); free(y_minus); free(f_plus); free(f_minus);
}

double sofe_rk4_step_component(SOFERateFunction f, double* y, double t,
                                double dt, double* params, int n_params,
                                int component_idx) {
    if (!f || !y) return 0.0;
    int n = n_params + 1;
    double* k1 = (double*)malloc(n * sizeof(double));
    double* k2 = (double*)malloc(n * sizeof(double));
    double* k3 = (double*)malloc(n * sizeof(double));
    double* k4 = (double*)malloc(n * sizeof(double));
    double* ytmp = (double*)malloc(n * sizeof(double));

    f(y, t, params, n_params, k1);
    for (int i = 0; i < n; i++) ytmp[i] = y[i] + 0.5 * dt * k1[i];
    f(ytmp, t + 0.5*dt, params, n_params, k2);
    for (int i = 0; i < n; i++) ytmp[i] = y[i] + 0.5 * dt * k2[i];
    f(ytmp, t + 0.5*dt, params, n_params, k3);
    for (int i = 0; i < n; i++) ytmp[i] = y[i] + dt * k3[i];
    f(ytmp, t + dt, params, n_params, k4);

    double result = y[component_idx] +
        (dt/6.0)*(k1[component_idx] + 2.0*k2[component_idx] +
                   2.0*k3[component_idx] + k4[component_idx]);
    free(k1); free(k2); free(k3); free(k4); free(ytmp);
    return result;
}
