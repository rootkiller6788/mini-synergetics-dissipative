#include "pds_thermo.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ============================================================================
 * PDS Thermodynamics ? Entropy Production, Onsager Relations, GP Stability
 *
 * Key References:
 *   Prigogine, I. (1947). Etude Thermodynamique des Phenomenes Irreversibles.
 *   Glansdorff, P. & Prigogine, I. (1971). Thermodynamic Theory of Structure,
 *     Stability and Fluctuations. Wiley-Interscience.
 *   Onsager, L. (1931). Reciprocal Relations in Irreversible Processes. I, II.
 *     Physical Review, 37, 405-426; 38, 2265-2279.
 *   De Donder, Th. & Van Rysselberghe, P. (1936). Thermodynamic Theory of
 *     Affinity. Stanford University Press.
 * ============================================================================ */

/* ============================================================================
 * Entropy Balance ? dS = d_i_S + d_e_S (Prigogine's Central Equation)
 *
 * L4: Second Law for Open Systems
 * The total entropy change dS is the sum of:
 *   d_i_S: entropy produced irreversibly inside the system (>= 0)
 *   d_e_S: entropy exchanged with the environment (can be + or -)
 *
 * For a dissipative structure to exist:
 *   d_e_S < 0 (negentropy import) and |d_e_S| > d_i_S
 *   => dS < 0 (system decreases its entropy while environment entropy increases)
 *   => "Order through fluctuations" (Prigogine & Stengers, 1984)
 * ============================================================================ */

void pds_entropy_balance_update(PDSDissipativeSystem* sys, double dt,
                                 double env_temperature,
                                 const double* env_chem_potentials) {
    if (!sys || dt <= 0.0 || env_temperature <= 0.0) return;

    /* Compute internal entropy production P = d_i_S/dt */
    double P = pds_entropy_production_rate(sys, env_temperature);
    double d_i_S = P * dt;

    /* Compute entropy exchange with environment.
     * d_e_S = -sum_k (J_s,k / T) * A_k * dt
     * where J_s,k is the flux of species k across the boundary.
     * Simplified model: entropy flow proportional to concentration gradients
     * between system and environment.
     */
    double d_e_S = 0.0;
    if (env_chem_potentials && sys->x && sys->x_equilibrium) {
        for (int i = 0; i < sys->n_species; i++) {
            double conc_gradient = pds_vector_get(sys->x, i) -
                                   pds_vector_get(sys->x_equilibrium, i);
            /* Entropy flow: - (mu_sys - mu_env)/T * flux_coeff * dt */
            double mu_diff = env_chem_potentials[i];
            d_e_S -= (mu_diff / env_temperature) * conc_gradient * 0.1 * dt;
        }
    }

    /* If no environment data, use a simplified negative entropy flow
     * to model the open system. A dissipative structure must have d_e_S < 0
     * large enough to overcome internal entropy production. */
    if (!env_chem_potentials) {
        d_e_S = -1.2 * d_i_S; /* Negentropy import: system exports entropy */
    }

    double dS = d_i_S + d_e_S;

    sys->entropy.d_i_S += d_i_S;
    sys->entropy.d_e_S += d_e_S;
    sys->entropy.S += dS;
    sys->entropy.dS_dt = (dt > 0.0) ? dS / dt : 0.0;
    sys->entropy.entropy_production_rate = P;
    sys->entropy.entropy_flow_rate = (dt > 0.0) ? d_e_S / dt : 0.0;
}

/* ============================================================================
 * Entropy Production Rate ? P = sum_i J_i * X_i >= 0
 *
 * L4: Second Law in Local Form
 * sigma = sum_i J_i * X_i >= 0  (entropy production per unit volume)
 *
 * The thermodynamic forces X_i drive the fluxes J_i.
 * Common force-flux pairs:
 *   - Heat flux J_q driven by grad(1/T)
 *   - Diffusion flux J_k driven by -grad(mu_k/T)
 *   - Chemical reaction rate v_r driven by affinity A_r/T
 *   - Viscous stress driven by velocity gradients
 * ============================================================================ */

double pds_entropy_production_rate(const PDSDissipativeSystem* sys,
                                    double temperature) {
    if (!sys || temperature <= 0.0) return 0.0;

    double P = 0.0;

    /* Chemical reaction contribution to entropy production.
     * For each reaction k:
     *   sigma_k = v_k * A_k / T
     * where v_k = reaction rate, A_k = -sum_i nu_{ik} mu_i = affinity.
     *
     * In simplified form: sigma_k = R * (v_f - v_r) * ln(v_f / v_r)
     * which is always >= 0 (Second Law).
     */
    for (int k = 0; k < sys->n_reactions; k++) {
        /* Reaction rate proportional to concentration product */
        double forward_rate = sys->rate_constants[k];
        for (int i = 0; i < sys->n_species; i++) {
            if (sys->stoichiometry[i][k] < 0) {
                /* Reactant ? negative stoichiometric coefficient */
                double conc = pds_vector_get(sys->x, i);
                if (conc < 0.0) conc = 0.0;
                forward_rate *= pow(conc, -sys->stoichiometry[i][k]);
            }
        }

        double reverse_rate = forward_rate * 0.01; /* Simplified equilibrium */

        /* Entropy production for this reaction:
         * sigma_k = R_gas * (v_f - v_r) * ln(v_f / v_r) >= 0 */
        if (forward_rate > 0.0 && reverse_rate > 0.0) {
            double net_rate = forward_rate - reverse_rate;
            double ratio = forward_rate / reverse_rate;
            if (ratio > 1e-10) {
                double R_gas = 8.314; /* J/(mol*K) */
                double sigma_k = fabs(net_rate) * R_gas * log(ratio);
                P += sigma_k;
            }
        }
    }

    /* Diffusion contribution to entropy production.
     * For species i: sigma_diff_i = D_i * (grad(c_i))^2 / (c_i * T)
     */
    if (sys->spatial_dims > 0 && sys->diffusion_coeffs) {
        for (int i = 0; i < sys->n_species; i++) {
            double ci = pds_vector_get(sys->x, i);
            /* Approximate |grad c| by (c - c_eq)/L */
            double grad_c = (ci - pds_vector_get(sys->x_equilibrium, i));
            if (fabs(ci) > 1e-10) {
                double sigma_diff = sys->diffusion_coeffs[i] *
                                    grad_c * grad_c / (fabs(ci) * temperature);
                P += sigma_diff;
            }
        }
    }

    return P / temperature; /* Normalize by temperature */
}

/* ============================================================================
 * Chemical Affinity (De Donder)
 *
 * L4: De Donder's Theorem
 * For reaction k: A_k = -sum_i nu_{ik} * mu_i
 * where nu_{ik} is the stoichiometric coefficient of species i in reaction k.
 *
 * The reaction proceeds in the forward direction iff A_k > 0.
 * At equilibrium: A_k = 0 for all reactions.
 * Entropy production: sigma = (1/T) * sum_k v_k * A_k >= 0
 * ============================================================================ */

double pds_chemical_affinity(const PDSDissipativeSystem* sys, int reaction_idx,
                              double temperature, const double* chem_potentials) {
    if (!sys || reaction_idx < 0 || reaction_idx >= sys->n_reactions) return 0.0;
    if (temperature <= 0.0) return 0.0;

    double A_k = 0.0;
    for (int i = 0; i < sys->n_species; i++) {
        /* nu_{ik}: negative for reactants, positive for products */
        double nu_ik = sys->stoichiometry[i][reaction_idx];
        double mu_i = chem_potentials ? chem_potentials[i] : 0.0;

        /* Standard chemical potential approximation:
         * mu_i = mu_i^0 + RT * ln(c_i / c_0) */
        if (!chem_potentials) {
            double ci = pds_vector_get(sys->x, i);
            if (ci < 1e-10) ci = 1e-10;
            double R_gas = 8.314;
            double c0 = 1.0; /* Standard concentration */
            mu_i = R_gas * temperature * log(ci / c0);
        }
        A_k -= nu_ik * mu_i; /* A_k = -sum_i nu_{ik} * mu_i */
    }
    return A_k;
}

/* ============================================================================
 * Onsager Reciprocal Relations
 *
 * L4: Onsager's Theorem (1931)
 * Near thermodynamic equilibrium, the phenomenological coefficients satisfy:
 *   L_{ij} = L_{ji}  (symmetry)
 *
 * This follows from microscopic reversibility and fluctuation theory.
 * For variables with different time-reversal parity (Casimir, 1945):
 *   L_{ij} = ?_i ?_j L_{ji}  where ?_i = +/-1 under time reversal.
 *
 * The entropy production is:
 *   sigma = sum_{i,j} L_{ij} X_i X_j >= 0
 * which requires L to be positive definite.
 * ============================================================================ */

PDSOnsagerMatrix* pds_onsager_matrix_create(int n_pairs) {
    if (n_pairs <= 0) return NULL;

    PDSOnsagerMatrix* om = (PDSOnsagerMatrix*)calloc(1, sizeof(PDSOnsagerMatrix));
    if (!om) return NULL;

    om->n_pairs = n_pairs;
    om->L = (double**)malloc(n_pairs * sizeof(double*));
    if (!om->L) { free(om); return NULL; }

    for (int i = 0; i < n_pairs; i++) {
        om->L[i] = (double*)calloc(n_pairs, sizeof(double));
        if (!om->L[i]) {
            for (int j = 0; j < i; j++) free(om->L[j]);
            free(om->L);
            free(om);
            return NULL;
        }
    }
    om->is_symmetric = true; /* Initialized to zero => symmetric */
    om->is_positive_definite = false;
    om->force_types = NULL;
    om->flux_types = NULL;
    return om;
}

void pds_onsager_matrix_free(PDSOnsagerMatrix* om) {
    if (!om) return;
    for (int i = 0; i < om->n_pairs; i++) {
        free(om->L[i]);
    }
    free(om->L);
    free(om->force_types);
    free(om->flux_types);
    free(om);
}

void pds_onsager_matrix_set(PDSOnsagerMatrix* om, int i, int j, double value) {
    if (!om || i < 0 || i >= om->n_pairs || j < 0 || j >= om->n_pairs) return;
    om->L[i][j] = value;
    /* Check Onsager symmetry */
    if (i != j) {
        om->is_symmetric = om->is_symmetric &&
                          (fabs(om->L[i][j] - om->L[j][i]) < 1e-12);
    }
}

void pds_onsager_matrix_check_symmetry(PDSOnsagerMatrix* om) {
    if (!om) return;
    om->is_symmetric = true;
    for (int i = 0; i < om->n_pairs && om->is_symmetric; i++) {
        for (int j = i + 1; j < om->n_pairs; j++) {
            if (fabs(om->L[i][j] - om->L[j][i]) > 1e-12) {
                om->is_symmetric = false;
                break;
            }
        }
    }
}

bool pds_onsager_matrix_is_positive_definite(const PDSOnsagerMatrix* om) {
    if (!om || om->n_pairs == 0) return false;

    /* Check positive definiteness using Sylvester's criterion:
     * All leading principal minors must be > 0.
     * For efficiency, we only check 1x1 and 2x2 cases numerically,
     * and use the trace-plus-determinant test for 2x2. */
    if (om->n_pairs == 1) {
        return om->L[0][0] > 0.0;
    }

    /* For n_pairs <= 3, check eigenvalues of the symmetric part.
     * L must be positive definite => (L + L^T)/2 pos def.
     * Since L is symmetric near equilibrium, check directly. */
    if (om->n_pairs == 2) {
        double a = om->L[0][0];
        double b = om->L[0][1];
        double d = om->L[1][1];
        /* For 2x2 symmetric: eigenvalues > 0 iff a>0 and a*d - b*b > 0 */
        double det = a * d - b * b;
        return (a > 0.0 && det > 0.0);
    }

    /* For larger matrices, use a numerical check:
     * For each vector e_i, check e_i^T L e_i > 0.
     * This is necessary but not sufficient, but practical. */
    bool diag_positive = true;
    for (int i = 0; i < om->n_pairs; i++) {
        if (om->L[i][i] <= 0.0) {
            diag_positive = false;
            break;
        }
    }
    return diag_positive;
}

void pds_onsager_matrix_compute_fluxes(const PDSOnsagerMatrix* om,
                                        const double* forces, double* fluxes) {
    if (!om || !forces || !fluxes) return;
    for (int i = 0; i < om->n_pairs; i++) {
        fluxes[i] = 0.0;
        for (int j = 0; j < om->n_pairs; j++) {
            fluxes[i] += om->L[i][j] * forces[j];
        }
    }
}

/* ============================================================================
 * Second Law Numerical Check
 * ============================================================================ */

bool pds_check_second_law(const PDSDissipativeSystem* sys, double tolerance) {
    if (!sys) return true; /* Vacuously true */
    /* Second Law: d_i_S must be non-negative (accounting for float errors) */
    return sys->entropy.d_i_S >= -tolerance;
}

/* ============================================================================
 * Minimum Entropy Production Theorem Check
 *
 * L4: Prigogine's Minimum Entropy Production Theorem (1945)
 *
 * Near equilibrium, steady states minimize the entropy production P.
 *   dP/dt <= 0  (P decreases over time)
 *   dP/dt = 0 at steady state
 *   d^2P/dt^2 >= 0 (P is at a minimum, not maximum)
 *
 * This theorem holds ONLY in the linear regime (near equilibrium).
 * FAR from equilibrium, the theorem breaks down ? dissipative structures
 * can and do produce MORE entropy than the corresponding equilibrium state.
 *
 * "Non-equilibrium may be a source of order" ? Prigogine (1977 Nobel Lecture)
 * ============================================================================ */

PDSMinEntropyProduction* pds_minimum_entropy_production_check(
    const PDSDissipativeSystem* sys) {
    if (!sys) return NULL;

    PDSMinEntropyProduction* mep =
        (PDSMinEntropyProduction*)calloc(1, sizeof(PDSMinEntropyProduction));
    if (!mep) return NULL;

    mep->entropy_production = sys->entropy.entropy_production_rate;
    mep->minimum_possible_P = 0.0; /* Equilibrium: P = 0 */

    /* At equilibrium, all forces X_i = 0 => P = 0 */
    mep->distance_from_minimum = mep->entropy_production;
    mep->is_steady_state = (sys->entropy.dS_dt < 1e-8);
    mep->is_minimum = false;

    /* Near equilibrium: if entropy production is small relative to
     * typical thermodynamic rates, the system is near minimum P. */
    if (mep->entropy_production < 1e-3) {
        mep->is_minimum = mep->is_steady_state;
    }

    /* Compute gradient dP/d(lambda_i) for each parameter */
    mep->n_gradient_components = sys->n_parameters;
    if (sys->n_parameters > 0) {
        mep->gradient_dP_dlambda =
            (double*)calloc(sys->n_parameters, sizeof(double));
        /* Approximate gradient via P's dependence on parameters.
         * Near equilibrium: P is quadratic in forces,
         * forces are linear in parameters => dP/dlambda approx linear. */
        for (int i = 0; i < sys->n_parameters; i++) {
            double lambda = sys->parameters[i];
            /* Simple model: P(lambda) ~ P_0 + alpha * (lambda - lambda_c)^2 */
            mep->gradient_dP_dlambda[i] = 2.0 * lambda * mep->entropy_production;
        }
    }

    mep->dP_dt = 0.0; /* Cannot determine without time derivatives */
    return mep;
}

void pds_min_entropy_production_free(PDSMinEntropyProduction* mep) {
    if (!mep) return;
    free(mep->gradient_dP_dlambda);
    free(mep);
}

/* ============================================================================
 * Glansdorff-Prigogine Stability Criterion
 *
 * L4: Glansdorff-Prigogine Universal Stability Criterion (1954, 1971)
 *
 * The excess entropy delta^2 S is always negative (maximum entropy at
 * equilibrium). The time derivative:
 *   (1/2) d/dt (delta^2 S) = delta^2 P = sum_k delta J_k delta X_k
 *
 * Stability condition:
 *   delta^2 P >= 0  =>  reference state is stable
 *   delta^2 P < 0   =>  reference state may become unstable (bifurcation)
 *
 * delta^2 S is a Lyapunov function for the thermodynamic branch:
 *   delta^2 S < 0  (negative definite)
 *   d/dt(delta^2 S) >= 0  =>  stable
 *   d/dt(delta^2 S) < 0  =>  unstable, bifurcation possible
 *
 * This is central to Prigogine's theory: when delta^2 P changes sign,
 * the thermodynamic branch becomes unstable and a dissipative structure
 * can emerge through a bifurcation.
 * ============================================================================ */

PDSGlansdorffPrigogineStability* pds_glansdorff_prigogine_stability(
    const PDSDissipativeSystem* sys) {
    if (!sys || sys->n_species == 0) return NULL;

    PDSGlansdorffPrigogineStability* gps =
        (PDSGlansdorffPrigogineStability*)calloc(1,
            sizeof(PDSGlansdorffPrigogineStability));
    if (!gps) return NULL;

    gps->n_variables = sys->n_species;

    /* Compute excess entropy: delta^2 S = -sum_i (1/c_i^eq) * (delta c_i)^2
     * For ideal mixtures, delta^2 S is always negative.
     * This ensures stability at equilibrium. */
    gps->excess_entropy = 0.0;
    for (int i = 0; i < sys->n_species; i++) {
        double ci = pds_vector_get(sys->x, i);
        double ci_eq = pds_vector_get(sys->x_equilibrium, i);
        if (ci_eq < 1e-10) ci_eq = 1e-10;
        double delta_ci = ci - ci_eq;
        /* delta^2 S = -(1/2) * sum_i (delta c_i)^2 / c_i^eq <= 0 */
        gps->excess_entropy -= 0.5 * delta_ci * delta_ci / ci_eq;
    }

    /* Compute excess entropy production:
     * delta^2 P = sum_{i,j} (dJ_i/dx_j)(dx_j) * (dX_i/dx_k)(dx_k)
     * Approximated via Jacobians of fluxes and forces. */
    gps->excess_entropy_production = 0.0;

    /* Build Jacobian of fluxes dJ_i/dx_j (simplified) */
    gps->jacobian_of_fluxes = (double**)malloc((size_t)sys->n_species * sizeof(double*));
    gps->jacobian_of_forces = (double**)malloc((size_t)sys->n_species * sizeof(double*));
    if (!gps->jacobian_of_fluxes || !gps->jacobian_of_forces) {
        free(gps->jacobian_of_fluxes);
        free(gps->jacobian_of_forces);
        free(gps);
        return NULL;
    }

    for (int i = 0; i < sys->n_species; i++) {
        gps->jacobian_of_fluxes[i] = (double*)calloc(sys->n_species, sizeof(double));
        gps->jacobian_of_forces[i] = (double*)calloc(sys->n_species, sizeof(double));
        if (!gps->jacobian_of_fluxes[i] || !gps->jacobian_of_forces[i]) {
            for (int j = 0; j <= i; j++) {
                free(gps->jacobian_of_fluxes[j]);
                free(gps->jacobian_of_forces[j]);
            }
            free(gps->jacobian_of_fluxes);
            free(gps->jacobian_of_forces);
            free(gps);
            return NULL;
        }
    }

    /* Simplify: assume fluxes proportional to concentration deviations.
     * delta J_i = sum_j L_ij * delta X_j
     * For chemical systems near equilibrium, this is linear. */
    for (int i = 0; i < sys->n_species; i++) {
        for (int j = 0; j < sys->n_species; j++) {
            /* Diagonal dominance: flux i mainly depends on its own force */
            gps->jacobian_of_fluxes[i][j] = (i == j) ? 1.0 : 0.1;
            gps->jacobian_of_forces[i][j] = (i == j) ? 1.0 : 0.05;
        }
    }

    /* Compute delta^2 P = sum_{i,j} delta_J_i * delta_X_j
     * Using the deviations from equilibrium. */
    for (int i = 0; i < sys->n_species; i++) {
        double delta_X_i = pds_vector_get(sys->x, i) -
                           pds_vector_get(sys->x_equilibrium, i);
        for (int j = 0; j < sys->n_species; j++) {
            double delta_J_i_j = gps->jacobian_of_fluxes[i][j];
            double delta_X_j = pds_vector_get(sys->x, j) -
                               pds_vector_get(sys->x_equilibrium, j);
            gps->excess_entropy_production += delta_J_i_j * delta_X_i * delta_X_j;
        }
    }

    /* GP stability criterion:
     * delta^2 S < 0 always.
     * If delta^2 P >= 0: stable (Lyapunov: d/dt(-delta^2 S) >= 0)
     * If delta^2 P < 0: possible instability => dissipative structure */
    gps->is_stable_by_gp = (gps->excess_entropy_production >= 0.0);
    gps->is_lyapunov_function = (gps->excess_entropy < 0.0);
    gps->gp_stability_margin = gps->excess_entropy_production;

    return gps;
}

void pds_gp_stability_free(PDSGlansdorffPrigogineStability* gps) {
    if (!gps) return;
    if (gps->jacobian_of_fluxes) {
        for (int i = 0; i < gps->n_variables; i++) {
            free(gps->jacobian_of_fluxes[i]);
        }
        free(gps->jacobian_of_fluxes);
    }
    if (gps->jacobian_of_forces) {
        for (int i = 0; i < gps->n_variables; i++) {
            free(gps->jacobian_of_forces[i]);
        }
        free(gps->jacobian_of_forces);
    }
    free(gps);
}

/* ============================================================================
 * Local Equilibrium Hypothesis Check
 *
 * L2: The local equilibrium hypothesis is a fundamental assumption of
 * irreversible thermodynamics (Prigogine, 1947). It states that locally,
 * thermodynamic state variables satisfy the same relations as at equilibrium,
 * even though global gradients drive irreversible processes.
 *
 * Validity condition: Knudsen number = mean free path / gradient length << 1
 * ============================================================================ */

bool pds_local_equilibrium_check(const PDSDissipativeSystem* sys,
                                  double max_gradient_norm) {
    if (!sys) return true;

    /* Check if spatial concentration gradients are small enough */
    double max_grad = 0.0;
    for (int i = 0; i < sys->n_species; i++) {
        double ci = pds_vector_get(sys->x, i);
        double ci_eq = pds_vector_get(sys->x_equilibrium, i);
        double grad = fabs(ci - ci_eq);
        if (grad > max_grad) max_grad = grad;
    }

    /* Local equilibrium holds if gradients are bounded */
    return (max_grad <= max_gradient_norm);
}

/* ============================================================================
 * Entropy Production Profile (Spatial)
 *
 * For spatially extended systems, compute sigma(x) = sum_i J_i(x) X_i(x)
 * at each spatial grid point. This identifies where entropy is produced
 * most intensely ? often at the boundaries or in the dissipative structure.
 * ============================================================================ */

PDSEntropyProductionProfile* pds_entropy_production_profile_compute(
    const PDSDissipativeSystem* sys, double temperature) {
    if (!sys || sys->n_grid_points <= 0) return NULL;

    PDSEntropyProductionProfile* epp =
        (PDSEntropyProductionProfile*)calloc(1,
            sizeof(PDSEntropyProductionProfile));
    if (!epp) return NULL;

    epp->n_cells = sys->n_grid_points;
    epp->cells = (PDSLocalEquilibriumCell*)calloc(
        epp->n_cells, sizeof(PDSLocalEquilibriumCell));
    epp->entropy_production_profile = (double*)calloc(epp->n_cells, sizeof(double));
    epp->entropy_flow_profile = (double*)calloc(epp->n_cells, sizeof(double));

    if (!epp->cells || !epp->entropy_production_profile) {
        pds_entropy_profile_free(epp);
        return NULL;
    }

    epp->total_entropy_production = 0.0;
    epp->max_entropy_production = 0.0;

    for (int i = 0; i < epp->n_cells; i++) {
        epp->cells[i].temperature = temperature;
        epp->cells[i].local_entropy_production = 0.0;
        epp->cells[i].in_local_equilibrium = true;

        /* Each cell's entropy production proportional to reactions
         * and diffusion at that spatial point. Simple model:
         * sigma(x_i) = sigma_0 * exp(-(x_i - x_center)^2 / width) */
        if (sys->spatial_grid) {
            double x_pos = sys->spatial_grid[i];
            double center = (sys->spatial_grid[0] + sys->spatial_grid[sys->n_grid_points - 1]) / 2.0;
            double width = (sys->spatial_grid[sys->n_grid_points - 1] - sys->spatial_grid[0]) / 4.0;
            if (width < 1e-10) width = 1.0;
            double dx = (x_pos - center) / width;
            double sigma = exp(-dx * dx) * 0.01;
            epp->entropy_production_profile[i] = sigma;
            epp->cells[i].local_entropy_production = sigma;
        } else {
            epp->entropy_production_profile[i] = 0.001;
        }

        epp->total_entropy_production += epp->entropy_production_profile[i];

        if (epp->entropy_production_profile[i] > epp->max_entropy_production) {
            epp->max_entropy_production = epp->entropy_production_profile[i];
            epp->max_production_cell_index = i;
        }
    }

    epp->mean_entropy_production = epp->total_entropy_production / epp->n_cells;
    return epp;
}

void pds_entropy_profile_free(PDSEntropyProductionProfile* epp) {
    if (!epp) return;
    free(epp->cells);
    free(epp->entropy_production_profile);
    free(epp->entropy_flow_profile);
    free(epp);
}

/* ============================================================================
 * Extended Irreversible Thermodynamics (EIT)
 *
 * L8: EIT (Jou, Casas-Vazquez, Lebon, 1993+) extends classical irreversible
 * thermodynamics by including fluxes as independent variables. This leads to:
 *
 *   tau * dJ/dt + J = L * X  (Maxwell-Cattaneo-Vernotte equation)
 *
 * Instead of the parabolic J = L * X, we get hyperbolic equations with
 * finite propagation speed for thermal and viscous signals.
 *
 * This resolves the "paradox of heat conduction" ? infinite propagation
 * speed in the Fourier law ? and is important for:
 *   - Nanoscale heat transport
 *   - Second sound in solids (helium II)
 *   - Fast chemical reactions
 *   - Relativistic thermodynamics
 * ============================================================================ */

PDSExtendedThermo* pds_extended_thermo_create(int n_extended_vars) {
    if (n_extended_vars <= 0) return NULL;

    PDSExtendedThermo* eit = (PDSExtendedThermo*)calloc(1, sizeof(PDSExtendedThermo));
    if (!eit) return NULL;

    eit->n_extended_vars = n_extended_vars;
    eit->extended_fluxes = (double*)calloc(n_extended_vars, sizeof(double));
    eit->relaxation_times = (double*)calloc(n_extended_vars, sizeof(double));

    if (!eit->extended_fluxes || !eit->relaxation_times) {
        free(eit->extended_fluxes);
        free(eit->relaxation_times);
        free(eit);
        return NULL;
    }

    /* Default relaxation times: fast fluxes have small tau */
    for (int i = 0; i < n_extended_vars; i++) {
        eit->relaxation_times[i] = 0.001 * (i + 1);
    }

    eit->generalized_entropy = 0.0;
    eit->generalized_entropy_production = 0.0;
    eit->hyperbolic_heat = false;
    return eit;
}

void pds_extended_thermo_free(PDSExtendedThermo* eit) {
    if (!eit) return;
    free(eit->extended_fluxes);
    free(eit->relaxation_times);
    free(eit);
}

void pds_extended_thermo_evolve(PDSExtendedThermo* eit, double dt) {
    if (!eit || dt <= 0.0) return;

    /* Relaxation evolution: tau * dJ/dt + J = L * X
     * Discretized: J(t+dt) = J(t) + (dt/tau) * (L*X - J(t))
     * This gives finite propagation speed c = sqrt(D/tau). */
    for (int i = 0; i < eit->n_extended_vars; i++) {
        double tau = eit->relaxation_times[i];
        if (tau < 1e-15) tau = 1e-15;
        double J_i = eit->extended_fluxes[i];
        /* Stationary flux = 1.0 (normalized), current J_i relaxes toward it */
        double J_stationary = 1.0;
        double dJ = (dt / tau) * (J_stationary - J_i);
        eit->extended_fluxes[i] += dJ;

        /* Generalized entropy production:
         * sigma_EIT = sigma_classical + sum_i (tau_i / 2L_ii) * (dJ_i/dt)^2 */
        eit->generalized_entropy_production +=
            (tau / 2.0) * dJ * dJ / (dt * dt);
    }

    /* Hyperbolic heat conduction condition:
     * tau_q > 0 => finite thermal wave speed */
    eit->hyperbolic_heat = (eit->n_extended_vars > 0 &&
                            eit->relaxation_times[0] > 1e-12);
}