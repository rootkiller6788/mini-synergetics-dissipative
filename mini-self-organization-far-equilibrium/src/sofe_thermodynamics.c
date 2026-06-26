/* ============================================================================
 * sofe_thermodynamics.c -- Non-Equilibrium Thermodynamics for Self-Organization
 *
 * Implements:
 *   - Entropy production rate computation (direct, reaction-network, master eq)
 *   - Onsager reciprocal relations and transport matrix
 *   - Rayleigh-Onsager dissipation function
 *   - Glansdorff-Prigogine stability criterion
 *   - Fluctuation-Dissipation Theorem verification
 *   - Crooks Fluctuation Theorem and Jarzynski Equality
 *   - Thermodynamic Uncertainty Relations
 *   - Non-Equilibrium Steady State (NESS) characterization
 *
 * Key references:
 *   Onsager, Phys. Rev. 37:405 (1931); 38:2265 (1931)
 *   Prigogine, Thermodynamics of Irreversible Processes (1955)
 *   Glansdorff & Prigogine, Thermodynamic Theory of Structure... (1971)
 *   Crooks, Phys. Rev. E 60:2721 (1999)
 *   Jarzynski, Phys. Rev. Lett. 78:2690 (1997)
 *   Barato & Seifert, Phys. Rev. Lett. 114:158101 (2015)
 * ============================================================================ */

#include "sofe_core.h"
#include "sofe_thermodynamics.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

/* ============================================================================
 * Entropy Production from Thermodynamic Forces and Fluxes
 *
 * Fundamental formula: sigma = d_i S/dt = sum_k F_k * J_k >= 0
 * where F_k are thermodynamic forces (gradients) and J_k are conjugate fluxes.
 *
 * For chemical reactions: sigma = (1/T) * sum_r A_r * v_r
 * where A_r = -sum_i nu_{ir} * mu_i is the chemical affinity.
 * ============================================================================ */

double sofe_entropy_production_direct(SOFEForceFluxSystem* ffs) {
    if (!ffs) return 0.0;
    double sigma = 0.0;
    for (int i = 0; i < ffs->n_pairs; i++) {
        sigma += ffs->pairs[i].force * ffs->pairs[i].flux;
    }
    return sigma >= 0.0 ? sigma : 0.0; /* Second Law: sigma >= 0 */
}

double sofe_entropy_production_reactions(SOFEReactionNetwork* net, double T) {
    if (!net || T <= 0.0) return 0.0;
    if (!net->flux_vector) sofe_network_compute_fluxes(net);
    int nr = net->n_reactions;
    double sigma = 0.0;
    const double R = 8.314462618;
    for (int r = 0; r < nr; r++) {
        /* Affinity A_r = -sum_i N[i][r] * mu_i */
        double A_r = 0.0;
        for (int i = 0; i < net->n_species; i++) {
            double N_ir = net->stoichiometric_matrix ?
                net->stoichiometric_matrix[i * nr + r] : 0.0;
            double conc = net->species[i].concentration;
            double mu_i = 0.0;
            if (conc > 1e-15) mu_i = R * T * log(conc);
            A_r -= N_ir * mu_i;
        }
        sigma += A_r * net->flux_vector[r];
    }
    return sigma / T;
}

double sofe_entropy_production_master_equation(double* p, double** W,
                                                int n_states, double k_B) {
    if (!p || !W || n_states <= 0 || k_B <= 0.0) return 0.0;
    double sigma = 0.0;
    for (int i = 0; i < n_states; i++) {
        for (int j = 0; j < n_states; j++) {
            if (i == j) continue;
            double flux_ij = W[i][j] * p[j];
            double flux_ji = W[j][i] * p[i];
            if (flux_ij > 1e-300 && flux_ji > 1e-300) {
                sigma += (flux_ij - flux_ji) * log(flux_ij / flux_ji);
            }
        }
    }
    return 0.5 * k_B * sigma; /* Each pair counted twice */
}

double sofe_entropy_flow(double T, double heat_flow,
                          double* chem_potentials, double* particle_flows,
                          int n_species) {
    if (T <= 0.0) return 0.0;
    double d_e_S = heat_flow / T;
    if (chem_potentials && particle_flows) {
        for (int i = 0; i < n_species; i++) {
            d_e_S -= chem_potentials[i] * particle_flows[i] / T;
        }
    }
    return d_e_S;
}

/* ============================================================================
 * Onsager Reciprocal Relations
 *
 * L_ij = L_ji  (without magnetic field)
 * L_ij(B) = L_ji(-B)  (with magnetic field)
 *
 * The Onsager matrix relates thermodynamic forces X to fluxes J:
 *   J_i = sum_j L_ij * X_j
 *
 * The Rayleigh-Onsager dissipation function:
 *   Phi = (1/2) * sum_{i,j} L_ij * X_i * X_j >= 0
 * This is related to the entropy production: sigma = 2*Phi / T
 * ============================================================================ */

int sofe_onsager_create(SOFEOnsagerMatrix* ons, int n_forces) {
    if (!ons || n_forces <= 0) return -1;
    memset(ons, 0, sizeof(SOFEOnsagerMatrix));
    ons->n_forces = n_forces;
    ons->L = (double*)malloc(n_forces * n_forces * sizeof(double));
    ons->thermodynamic_forces = (double*)malloc(n_forces * sizeof(double));
    ons->thermodynamic_fluxes = (double*)malloc(n_forces * sizeof(double));
    assert(ons->L && ons->thermodynamic_forces && ons->thermodynamic_fluxes);
    memset(ons->L, 0, n_forces * n_forces * sizeof(double));
    memset(ons->thermodynamic_forces, 0, n_forces * sizeof(double));
    memset(ons->thermodynamic_fluxes, 0, n_forces * sizeof(double));
    ons->is_symmetric = true;
    ons->magnetic_field_present = false;
    return 0;
}

void sofe_onsager_free(SOFEOnsagerMatrix* ons) {
    if (!ons) return;
    free(ons->L);
    free(ons->thermodynamic_forces);
    free(ons->thermodynamic_fluxes);
    memset(ons, 0, sizeof(SOFEOnsagerMatrix));
}

void sofe_onsager_set_coefficient(SOFEOnsagerMatrix* ons, int i, int j,
                                   double value) {
    if (!ons || i < 0 || i >= ons->n_forces || j < 0 || j >= ons->n_forces)
        return;
    ons->L[i * ons->n_forces + j] = value;
}

double sofe_onsager_get_coefficient(SOFEOnsagerMatrix* ons, int i, int j) {
    if (!ons || i < 0 || i >= ons->n_forces || j < 0 || j >= ons->n_forces)
        return 0.0;
    return ons->L[i * ons->n_forces + j];
}

bool sofe_onsager_verify_reciprocity(SOFEOnsagerMatrix* ons, double tolerance) {
    if (!ons) return false;
    int n = ons->n_forces;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double L_ij = ons->L[i * n + j];
            double L_ji = ons->L[j * n + i];
            if (fabs(L_ij - L_ji) > tolerance) return false;
        }
    }
    return true;
}

void sofe_onsager_compute_fluxes(SOFEOnsagerMatrix* ons) {
    if (!ons) return;
    int n = ons->n_forces;
    for (int i = 0; i < n; i++) {
        ons->thermodynamic_fluxes[i] = 0.0;
        for (int j = 0; j < n; j++) {
            ons->thermodynamic_fluxes[i] +=
                ons->L[i * n + j] * ons->thermodynamic_forces[j];
        }
    }
}

double sofe_rayleigh_onsager_dissipation(SOFEOnsagerMatrix* ons) {
    if (!ons) return 0.0;
    int n = ons->n_forces;
    double Phi = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Phi += 0.5 * ons->L[i * n + j] *
                   ons->thermodynamic_forces[i] *
                   ons->thermodynamic_forces[j];
        }
    }
    return Phi >= 0.0 ? Phi : 0.0;
}

/* ============================================================================
 * Glansdorff-Prigogine Stability Criterion
 *
 * The GP criterion states that for a system to be stable:
 *   d_X P = sum_k J_k * dX_k <= 0  (excess entropy production)
 *
 * The excess entropy delta^2 S is a Lyapunov functional:
 *   delta^2 S = -(1/2) * sum_i C_i * (delta X_i)^2 < 0 for stable systems
 *
 * At a bifurcation: d_X P changes sign, and delta^2 S -> 0.
 *
 * Reference: Glansdorff & Prigogine (1971), Chapter 6-7
 * ============================================================================ */

void sofe_gp_init(SOFEGlansdorffPrigogineCriterion* gp) {
    if (!gp) return;
    memset(gp, 0, sizeof(SOFEGlansdorffPrigogineCriterion));
    gp->is_stable = true;
}

void sofe_gp_compute(SOFEGlansdorffPrigogineCriterion* gp,
                      SOFEReactionNetwork* net, double T,
                      double* perturbation, int n_species) {
    if (!gp || !net || !perturbation || T <= 0.0) return;

    /* Compute excess entropy: delta^2 S from perturbation */
    gp->excess_entropy = 0.0;
    for (int i = 0; i < n_species && i < net->n_species; i++) {
        double dx = perturbation[i];
        gp->excess_entropy -= 0.5 * dx * dx / T;
    }

    /* Compute excess entropy production */
    gp->excess_entropy_production = 0.0;
    if (!net->flux_vector) sofe_network_compute_fluxes(net);
    if (!net->jacobian) sofe_network_compute_jacobian(net);

    /* d_X P = sum_k J_k * dX_k */
    if (net->flux_vector && net->jacobian) {
        for (int i = 0; i < n_species && i < net->n_species; i++) {
            for (int j = 0; j < net->n_species; j++) {
                gp->excess_entropy_production +=
                    net->jacobian[i * net->n_species + j] *
                    perturbation[i] * perturbation[j];
            }
        }
    }

    /* Stability: d_X P <= 0 for stable systems */
    gp->is_stable = (gp->excess_entropy_production <= gp->critical_threshold);
    gp->stability_margin = -gp->excess_entropy_production;
}

bool sofe_gp_is_stable(SOFEGlansdorffPrigogineCriterion* gp) {
    if (!gp) return true;
    return gp->is_stable;
}

double sofe_excess_entropy(double* state_perturbation,
                            double* heat_capacities, int n) {
    if (!state_perturbation || !heat_capacities || n <= 0) return 0.0;
    double d2S = 0.0;
    for (int i = 0; i < n; i++) {
        double ci = heat_capacities[i] > 1e-15 ? heat_capacities[i] : 1.0;
        d2S -= 0.5 * state_perturbation[i] * state_perturbation[i] / ci;
    }
    return d2S;
}

/* ============================================================================
 * Fluctuation-Dissipation Theorem
 *
 * Connects equilibrium fluctuations to linear response:
 *   S(omega) = (2*kT/omega) * Im[chi(omega)]
 *
 * where S(omega) is the power spectrum of fluctuations and
 * chi(omega) is the susceptibility (linear response function).
 *
 * Einstein relation for fluctuation probability:
 *   P(delta_S) ~ exp(delta_S / k_B)
 * ============================================================================ */

void sofe_fd_init(SOFEFluctuationDissipation* fd, int n_timepoints,
                   double dt, double T) {
    if (!fd || n_timepoints <= 0) return;
    memset(fd, 0, sizeof(SOFEFluctuationDissipation));
    fd->n_timepoints = n_timepoints;
    fd->dt = dt;
    fd->temperature = T;
    fd->correlation_function = (double*)malloc(n_timepoints * sizeof(double));
    fd->response_function = (double*)malloc(n_timepoints * sizeof(double));
    assert(fd->correlation_function && fd->response_function);
    memset(fd->correlation_function, 0, n_timepoints * sizeof(double));
    memset(fd->response_function, 0, n_timepoints * sizeof(double));
}

void sofe_fd_free(SOFEFluctuationDissipation* fd) {
    if (!fd) return;
    free(fd->correlation_function);
    free(fd->response_function);
    memset(fd, 0, sizeof(SOFEFluctuationDissipation));
}

/* Compute auto-correlation function: C(tau) = <A(t+tau) * A(t)>
 * Using time-averaging over the available time series. */
void sofe_fd_compute_correlation(SOFEFluctuationDissipation* fd,
                                  double* timeseries, int n_samples) {
    if (!fd || !timeseries || n_samples <= 0) return;
    int max_lag = fd->n_timepoints;
    if (max_lag > n_samples) max_lag = n_samples;

    /* Compute mean */
    double mean = 0.0;
    for (int i = 0; i < n_samples; i++) mean += timeseries[i];
    mean /= n_samples;

    /* Compute autocorrelation for each lag */
    for (int lag = 0; lag < max_lag; lag++) {
        double sum = 0.0;
        int count = 0;
        for (int t = 0; t < n_samples - lag; t++) {
            sum += (timeseries[t + lag] - mean) * (timeseries[t] - mean);
            count++;
        }
        fd->correlation_function[lag] = (count > 0) ? sum / count : 0.0;
    }
}

/* Compute linear response function from perturbation experiment.
 * chi(t) ~ (<A(t)>_perturbed - <A(t)>_unperturbed) / epsilon
 * where epsilon is the perturbation strength. */
void sofe_fd_compute_response(SOFEFluctuationDissipation* fd,
                               double* perturbed, double* unperturbed,
                               int n_samples, double epsilon) {
    if (!fd || !perturbed || !unperturbed || n_samples <= 0 ||
        fabs(epsilon) < 1e-15) return;
    int n = fd->n_timepoints;
    if (n > n_samples) n = n_samples;
    for (int i = 0; i < n; i++) {
        fd->response_function[i] = (perturbed[i] - unperturbed[i]) / epsilon;
    }
}

/* Verify FDT: compare correlation and response via the relation
 * dC(t)/dt = -kT * chi(t)  for t > 0  (classical limit) */
bool sofe_fd_verify(SOFEFluctuationDissipation* fd, double tolerance) {
    if (!fd) return false;
    int n = fd->n_timepoints;
    if (n < 2) return false;
    double kT = 1.380649e-23 * fd->temperature;
    double max_err = 0.0;
    for (int i = 1; i < n; i++) {
        double dC_dt = (fd->correlation_function[i] - fd->correlation_function[i-1]) / fd->dt;
        double expected = -kT * fd->response_function[i-1];
        double err = fabs(dC_dt - expected);
        if (err > max_err) max_err = err;
    }
    fd->fluctuation_dissipation_ratio = max_err / (kT + 1e-30);
    fd->fdt_holds = (max_err < tolerance);
    return fd->fdt_holds;
}

double sofe_einstein_fluctuation_probability(double delta_entropy,
                                              double k_B) {
    if (k_B <= 0.0) return 0.0;
    return exp(delta_entropy / k_B);
}

/* ============================================================================
 * Crooks Fluctuation Theorem (1999)
 *
 * P_f(W) / P_r(-W) = exp((W - dF) / kT)
 *
 * where:
 *   P_f(W) = probability of work W in forward process
 *   P_r(-W) = probability of work -W in reverse process
 *   dF = free energy difference between equilibrium states
 *
 * This generalizes the Jarzynski Equality and provides a stronger
 * constraint on non-equilibrium fluctuations than the Second Law alone.
 *
 * Reference: Crooks, Phys. Rev. E 60:2721 (1999)
 * ============================================================================ */

double sofe_crooks_ratio(double work, double free_energy_diff, double kT) {
    if (kT <= 0.0) return 0.0;
    return exp((work - free_energy_diff) / kT);
}

/* ============================================================================
 * Jarzynski Equality (1997)
 *
 * <exp(-W/kT)> = exp(-dF/kT)
 *
 * This equality holds for ANY non-equilibrium process, regardless of
 * how far from equilibrium the system is driven. It provides a way to
 * estimate equilibrium free energy differences from non-equilibrium
 * measurements.
 *
 * Estimator: dF_est = -kT * log( (1/N) * sum_i exp(-W_i/kT) )
 *
 * Reference: Jarzynski, Phys. Rev. Lett. 78:2690 (1997)
 * ============================================================================ */

double sofe_jarzynski_estimator(double* work_values, int n_trajectories,
                                 double kT) {
    if (!work_values || n_trajectories <= 0 || kT <= 0.0) return 0.0;
    double sum_exp = 0.0;
    for (int i = 0; i < n_trajectories; i++) {
        sum_exp += exp(-work_values[i] / kT);
    }
    double avg_exp = sum_exp / n_trajectories;
    if (avg_exp <= 0.0) return 0.0;
    return -kT * log(avg_exp);
}

/* ============================================================================
 * Thermodynamic Uncertainty Relation (Barato & Seifert, 2015)
 *
 * For any current J in a NESS:
 *   Var[J] / <J>^2 >= 2 * k_B / sigma
 *
 * where sigma is the total entropy production rate.
 * This gives a fundamental bound on the precision of currents
 * in terms of the thermodynamic cost.
 *
 * Returns the TUR ratio = (Var[J]/<J>^2) / (2*k_B/sigma)
 * A value >= 1 indicates the TUR is satisfied.
 *
 * Reference: Barato & Seifert, Phys. Rev. Lett. 114:158101 (2015)
 * ============================================================================ */

double sofe_thermodynamic_uncertainty_ratio(double current_mean,
                                              double current_variance,
                                              double entropy_prod_rate,
                                              double k_B) {
    if (fabs(current_mean) < 1e-30 || entropy_prod_rate <= 0.0 || k_B <= 0.0)
        return 0.0;
    double lhs = current_variance / (current_mean * current_mean);
    double rhs = 2.0 * k_B / entropy_prod_rate;
    if (rhs <= 0.0) return 0.0;
    return lhs / rhs;
}

/* ============================================================================
 * Non-Equilibrium Steady State (NESS) Characterization
 *
 * Housekeeping heat Q_hk: the heat dissipated to maintain the NESS
 * against the violation of detailed balance.
 *
 * For a Markov process with steady-state probabilities p_i^ss and
 * transition rates W_ij:
 *   Q_hk = kT * sum_{i,j} p_j^ss * W_ij * ln(W_ij / W_ji)
 *
 * The excess heat Q_ex drives transitions between different NESS.
 *
 * Reference: Oono & Paniconi, Prog. Theor. Phys. Suppl. 130:29 (1998)
 * ============================================================================ */

double sofe_housekeeping_heat(double* steady_state_prob,
                               double* transition_rates,
                               int n_states, double kT) {
    if (!steady_state_prob || !transition_rates || n_states <= 0 || kT <= 0.0)
        return 0.0;
    double Q_hk = 0.0;
    for (int i = 0; i < n_states; i++) {
        for (int j = 0; j < n_states; j++) {
            if (i == j) continue;
            double W_ij = transition_rates[i * n_states + j];
            double W_ji = transition_rates[j * n_states + i];
            if (W_ij > 1e-300 && W_ji > 1e-300) {
                Q_hk += steady_state_prob[j] * W_ij * log(W_ij / W_ji);
            }
        }
    }
    return kT * Q_hk;
}

double sofe_excess_heat(double total_heat, double housekeeping_heat) {
    return total_heat - housekeeping_heat;
}

/* ============================================================================
 * Thermodynamic State and Utility Functions
 * ============================================================================ */

void sofe_thermo_state_print(SOFEThermoState* ts) {
    if (!ts) { printf("(null thermo state)\n"); return; }
    printf("Thermodynamic State:\n");
    printf("  S = %.6e  U = %.6e  F = %.6e  G = %.6e\n",
           ts->entropy, ts->internal_energy,
           ts->helmholtz_free_energy, ts->gibbs_free_energy);
    printf("  T = %.2f K  P = %.2f Pa  V = %.6e\n",
           ts->temperature, ts->pressure, ts->volume);
    printf("  mu = %.6e  N = %.6e\n",
           ts->chemical_potential, ts->particle_number);
}

/* Boltzmann entropy: S = k_B * ln(W) where W is the multiplicity */
double sofe_boltzmann_entropy(double multiplicity, double k_B) {
    if (multiplicity <= 0.0 || k_B <= 0.0) return 0.0;
    return k_B * log(multiplicity);
}

/* Gibbs (Shannon) entropy: S = -k_B * sum_i p_i * ln(p_i) */
double sofe_gibbs_entropy(double* probabilities, int n_states, double k_B) {
    if (!probabilities || n_states <= 0 || k_B <= 0.0) return 0.0;
    double S = 0.0;
    for (int i = 0; i < n_states; i++) {
        double p = probabilities[i];
        if (p > 1e-300) S -= p * log(p);
    }
    return k_B * S;
}

/* Chemical potential for ideal gas: mu(T,P) = mu_ref(T) + RT ln(P/P_ref) */
double sofe_ideal_gas_chemical_potential(double T, double P,
                                           double P_ref, double mu_ref) {
    if (T <= 0.0 || P_ref <= 0.0) return mu_ref;
    const double R = 8.314462618;
    return mu_ref + R * T * log(P / P_ref);
}
