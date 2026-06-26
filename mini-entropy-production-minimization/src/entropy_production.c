/**
 * @file entropy_production.c
 * @brief Implementation of core entropy production computations.
 *
 * Implements the fundamental bilinear form sigma = sum J_i X_i,
 * volumetric entropy production for continuous media, entropy
 * balance, and Second Law verification.
 *
 * References:
 *   - de Groot & Mazur (1984) "Non-Equilibrium Thermodynamics"
 *   - Kondepudi & Prigogine (2014) "Modern Thermodynamics"
 *   - Prigogine (1947) "Etude thermodynamique des phenomenes irreversibles"
 */

#include "entropy_production.h"
#include <math.h>
#include <string.h>
#include <stdio.h>

/* Physical constants */
#define R_GAS 8.314462618     /**< Universal gas constant [J/(mol*K)] */
#define K_B  1.380649e-23     /**< Boltzmann constant [J/K] */
#define N_A  6.02214076e23    /**< Avogadro's number [1/mol] */

/* ---------------------------------------------------------------------------
 * L2: Core Concepts -- Entropy production rate from forces and fluxes
 * ------------------------------------------------------------------------- */

double ep_compute_sigma(const double forces[], const double fluxes[], int n)
{
    double sigma = 0.0;
    int i;

    if (forces == NULL || fluxes == NULL || n <= 0) {
        return 0.0;
    }

    /* sigma = sum_{i=0}^{n-1} J_i * X_i */
    for (i = 0; i < n; i++) {
        sigma += forces[i] * fluxes[i];
    }

    return sigma;
}

double ep_volumetric_sigma(const ThermodynamicFlux *flux,
                           const ThermodynamicForce *force,
                           double temperature)
{
    double sigma = 0.0;
    int k, rho, i, j;

    if (flux == NULL || force == NULL || temperature <= 0.0) {
        return 0.0;
    }

    /* 1. Heat conduction contribution: sigma_q = J_q * grad(1/T) */
    sigma += flux->heat_flux * force->thermal_force;

    /* 2. Diffusion contribution: sigma_diff = sum_k J_k * (-grad(mu_k/T)) */
    for (k = 0; k < flux->num_species && k < force->num_species; k++) {
        sigma += flux->mass_flux[k] * force->diffusion_force[k];
    }

    /* 3. Chemical reaction contribution: sigma_react = sum_rho r_rho * A_rho/T */
    for (rho = 0; rho < flux->num_reactions && rho < force->num_reactions; rho++) {
        sigma += flux->reaction_rate[rho] * force->affinity[rho] / temperature;
    }

    /* 4. Viscous dissipation: sigma_visc = Pi : grad(v)/T = sum_{i,j} Pi_{ij} * X_v_{ij} */
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            sigma += flux->momentum_flux_tensor[i][j]
                   * force->viscous_force[i][j];
        }
    }

    /* 5. Electrical conduction: sigma_elec = i * (-grad(phi)/T) */
    sigma += flux->electric_current * force->electrical_force;

    return sigma;
}

int ep_entropy_balance(const ThermodynamicState *state,
                       double heat_flow,
                       const double mass_flow[],
                       double work_done,
                       double *d_eS,
                       double *d_iS)
{
    double dS_total, dS_exchange;
    int k;

    if (state == NULL || d_eS == NULL || d_iS == NULL) {
        return -1;
    }

    if (state->temperature <= 0.0) {
        return -1;
    }

    /* Entropy exchange with surroundings:
     * d_eS = dQ/T + sum_k s_k * dN_k
     * The partial molar entropy s_k is approximated as the
     * ideal gas value: s_k = s_k^0 - R ln(x_k P/P0)
     * For simplicity, we use the entropy already stored in state.
     */
    dS_exchange = heat_flow / state->temperature;

    if (mass_flow != NULL) {
        for (k = 0; k < state->num_species; k++) {
            /* Approximate partial molar entropy:
             * s_k approx = S / (sum_j N_j), i.e., average molar entropy */
            double total_moles = 0.0;
            int j;
            for (j = 0; j < state->num_species; j++) {
                total_moles += state->mole_numbers[j];
            }
            if (total_moles > 0.0) {
                double avg_molar_entropy = state->entropy / total_moles;
                dS_exchange += mass_flow[k] * avg_molar_entropy;
            }
        }
    }

    /* Total entropy change from First Law:
     * dU = dQ + dW --> dS_total = dQ/T + d_iS
     * We compute d_iS from the entropy production:
     * If the process is reversible, d_iS = 0.
     * For irreversible processes, d_iS > 0.
     */
    *d_eS = dS_exchange;

    /* The internal entropy production depends on the process.
     * For a general irreversible process, we estimate:
     * d_iS = dQ * (1/T_system - 1/T_surroundings) + dissipation from work
     * Here we use the dissipative work: d_iS = (dW_diss) / T
     * dW_diss is the portion of work that goes into entropy production.
     *
     * Without a detailed process specification, we use the Clausius
     * inequality: d_iS >= 0. We estimate it from energy conservation.
     */
    /* Conservative estimate: entropy production from heat flow alone */
    *d_iS = 0.0;
    if (heat_flow > 0.0) {
        /* Heat flows from hot to cold; entropy production positive */
        *d_iS = heat_flow / state->temperature;
    }

    return 0;
}

int ep_check_second_law(double sigma, double epsilon)
{
    /* sigma >= -epsilon means the Second Law is satisfied */
    return (sigma >= -epsilon) ? 1 : 0;
}

/* ---------------------------------------------------------------------------
 * L3: Mathematical Structures -- Bilinear and quadratic forms
 * ------------------------------------------------------------------------- */

double ep_rayleigh_dissipative(const double fluxes[], const double *R, int n)
{
    double Phi = 0.0;
    int i, j;

    if (fluxes == NULL || R == NULL || n <= 0) {
        return 0.0;
    }

    /* Phi(J, J) = (1/2) * sum_{i,j} J_i * R_{ij} * J_j */
    for (i = 0; i < n; i++) {
        double row_sum = 0.0;
        for (j = 0; j < n; j++) {
            row_sum += R[i * n + j] * fluxes[j];
        }
        Phi += fluxes[i] * row_sum;
    }

    return 0.5 * Phi;
}

double ep_quadratic_form_forces(const double forces[], const double *L, int n)
{
    double sigma = 0.0;
    int i, j;

    if (forces == NULL || L == NULL || n <= 0) {
        return 0.0;
    }

    /* sigma = X^T L X = sum_{i,j} L_{ij} X_i X_j */
    for (i = 0; i < n; i++) {
        double row_sum = 0.0;
        for (j = 0; j < n; j++) {
            row_sum += L[i * n + j] * forces[j];
        }
        sigma += forces[i] * row_sum;
    }

    return sigma;
}

double ep_quadratic_form_fluxes(const double fluxes[], const double *R, int n)
{
    /* sigma = J^T R J = 2 * Phi(J, J) */
    return 2.0 * ep_rayleigh_dissipative(fluxes, R, n);
}

double ep_cross_entropy_production(const double forces1[],
                                   const double forces2[],
                                   const double *L, int n)
{
    double sigma12 = 0.0;
    int i, j;

    if (forces1 == NULL || forces2 == NULL || L == NULL || n <= 0) {
        return 0.0;
    }

    /* sigma_{12} = X1^T L X2 = sum_{i,j} L_{ij} (X1)_i (X2)_j */
    for (i = 0; i < n; i++) {
        double row_sum = 0.0;
        for (j = 0; j < n; j++) {
            row_sum += L[i * n + j] * forces2[j];
        }
        sigma12 += forces1[i] * row_sum;
    }

    return sigma12;
}

/* ---------------------------------------------------------------------------
 * L4: Fundamental Laws -- Prigogine's Minimum Entropy Production Theorem
 * ------------------------------------------------------------------------- */

double ep_prigogine_dsigma_dt(const double forces[],
                              const double dforces_dt[],
                              const double *L, int n)
{
    double dsigma_dt = 0.0;
    int i, j;

    if (forces == NULL || dforces_dt == NULL || L == NULL || n <= 0) {
        return 0.0;
    }

    /* d(sigma)/dt = 2 * sum_{i,j} L_{ij} X_i (dX_j/dt) */
    for (i = 0; i < n; i++) {
        double Li_dot_dXdt = 0.0;
        for (j = 0; j < n; j++) {
            Li_dot_dXdt += L[i * n + j] * dforces_dt[j];
        }
        dsigma_dt += 2.0 * forces[i] * Li_dot_dXdt;
    }

    return dsigma_dt;
}

int ep_is_mep_steady_state(const double forces[],
                           const double *L,
                           const int constrained[],
                           int n, double tol)
{
    int i, j;

    if (forces == NULL || L == NULL || n <= 0) {
        return 0;
    }

    /* Check partial(sigma)/partial(X_i) = 2 sum_j L_{ij} X_j = 0
     * for all unconstrained forces. */
    for (i = 0; i < n; i++) {
        if (constrained != NULL && constrained[i]) {
            continue;  /* Constrained force -- skip check */
        }

        double partial_sigma = 0.0;
        for (j = 0; j < n; j++) {
            partial_sigma += 2.0 * L[i * n + j] * forces[j];
        }

        if (fabs(partial_sigma) > tol) {
            return 0;
        }
    }

    return 1;
}

int ep_find_mep_steady_state(double forces[],
                             const double *L,
                             const int constrained[],
                             const double fixed_values[],
                             int n, int max_iter, double tol)
{
    double grad[64];    /* Gradient: partial(sigma)/partial(X_i) */
    double step[64];    /* Newton step */
    int iter;
    int i, j;

    if (forces == NULL || L == NULL || n <= 0 || n > 64) {
        return -1;
    }

    /* Iterative minimization: apply Newton step for unconstrained forces.
     * For quadratic sigma = X^T L X, one Newton step suffices,
     * but we iterate for robustness against constraints. */
    for (iter = 0; iter < max_iter; iter++) {
        /* Compute gradient: g_i = 2 sum_j L_{ij} X_j */
        double grad_norm = 0.0;
        for (i = 0; i < n; i++) {
            if (constrained != NULL && constrained[i]) {
                grad[i] = 0.0;
                continue;
            }
            grad[i] = 0.0;
            for (j = 0; j < n; j++) {
                grad[i] += 2.0 * L[i * n + j] * forces[j];
            }
            grad_norm += grad[i] * grad[i];
        }
        grad_norm = sqrt(grad_norm);

        if (grad_norm < tol) {
            return 0;  /* Converged */
        }

        /* Fixed step size gradient descent for unconstrained forces */
        double alpha = 1.0 / (2.0 * n);  /* Conservative step */
        for (i = 0; i < n; i++) {
            if (constrained != NULL && constrained[i]) {
                if (fixed_values != NULL) {
                    forces[i] = fixed_values[i];
                }
                continue;
            }
            forces[i] -= alpha * grad[i];
        }
    }

    return -1;  /* Did not converge within max_iter */
}

/* ---------------------------------------------------------------------------
 * Utility functions
 * ------------------------------------------------------------------------- */

void ep_state_init(ThermodynamicState *state, double temperature,
                   double pressure, int num_species)
{
    int k;

    if (state == NULL) return;

    memset(state, 0, sizeof(ThermodynamicState));
    state->temperature = temperature;
    state->pressure = pressure;
    state->volume = 1.0;            /* Default volume */
    state->internal_energy = 0.0;
    state->entropy = 0.0;
    state->helmholtz_free_energy = 0.0;
    state->gibbs_free_energy = 0.0;
    state->enthalpy = 0.0;
    state->num_species = (num_species > 0 && num_species <= EP_MAX_SPECIES)
                         ? num_species : 0;

    for (k = 0; k < EP_MAX_SPECIES; k++) {
        state->chemical_potential[k] = 0.0;
        state->mole_numbers[k] = 0.0;
        state->concentrations[k] = 0.0;
    }
}

void ep_state_copy(ThermodynamicState *dest, const ThermodynamicState *src)
{
    if (dest == NULL || src == NULL) return;
    memcpy(dest, src, sizeof(ThermodynamicState));
}

double ep_ideal_gas_entropy(const ThermodynamicState *state,
                            const double std_entropy[])
{
    double total_entropy = 0.0;
    double total_moles = 0.0;
    int k;

    if (state == NULL) return 0.0;

    for (k = 0; k < state->num_species; k++) {
        total_moles += state->mole_numbers[k];
    }

    if (total_moles <= 0.0) return 0.0;

    /* Sackur-Tetrode generalized for mixtures:
     * S = sum_k N_k [s_k^0 + c_{p,k} ln(T/T0) - R ln(x_k P/P0)]
     * For simplicity, we use the standard molar entropies directly
     * and add the mixing entropy term. */
    for (k = 0; k < state->num_species; k++) {
        double x_k = state->mole_numbers[k] / total_moles;
        double s_k = 0.0;

        if (std_entropy != NULL) {
            s_k = std_entropy[k];
        }

        /* Mixing entropy: -R * N_k * ln(x_k) */
        if (x_k > 0.0) {
            s_k -= R_GAS * log(x_k);
        }

        total_entropy += state->mole_numbers[k] * s_k;
    }

    return total_entropy;
}

void ep_ideal_gas_chemical_potentials(ThermodynamicState *state,
                                      const double std_potential[],
                                      double R)
{
    double total_moles = 0.0;
    int k;

    if (state == NULL) return;

    for (k = 0; k < state->num_species; k++) {
        total_moles += state->mole_numbers[k];
    }

    if (total_moles <= 0.0 || state->pressure <= 0.0) return;

    for (k = 0; k < state->num_species; k++) {
        double x_k = state->mole_numbers[k] / total_moles;
        double mu_k0 = (std_potential != NULL) ? std_potential[k] : 0.0;

        /* mu_k = mu_k^0(T) + R T ln(x_k P / P0)
         * Using P0 = 1 bar = 1e5 Pa as standard pressure */
        if (x_k > 0.0) {
            double P0 = 1.0e5;  /* Standard pressure [Pa] */
            state->chemical_potential[k] = mu_k0
                + R * state->temperature * log(x_k * state->pressure / P0);
        } else {
            state->chemical_potential[k] = -1.0e10;  /* Very negative */
        }
    }
}

double ep_reaction_affinity(const ChemicalReaction *reaction,
                            const double chem_pot[], int num_species)
{
    double A = 0.0;
    int k;

    if (reaction == NULL || chem_pot == NULL) return 0.0;

    /* A_rho = -sum_k nu_{k,rho} * mu_k
     * Note: nu_k < 0 for reactants, > 0 for products
     * Positive A means the reaction proceeds forward. */
    for (k = 0; k < num_species && k < EP_MAX_SPECIES; k++) {
        A -= reaction->stoichiometry[k] * chem_pot[k];
    }

    return A;
}