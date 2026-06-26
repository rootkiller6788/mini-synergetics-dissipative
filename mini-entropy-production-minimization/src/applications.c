/**
 * @file applications.c
 * @brief Application examples of entropy production minimization.
 *
 * Real-world applications of the Minimum Entropy Production Principle:
 *   - Climate/Earth energy balance model
 *   - Metabolic network thermodynamics in biological cells
 *   - Industrial chemical reactor optimization
 *   - Thermoelectric generator efficiency
 *   - Building thermal management
 *
 * Each application demonstrates how MEPP provides a variational
 * principle for predicting steady states in non-equilibrium systems.
 *
 * L7: Applications (5 examples)
 */

#include "entropy_production.h"
#include "steady_state.h"
#include "onsager_relations.h"
#include "dissipative_systems.h"
#include "nonlinear_thermo.h"
#include <stdio.h>
#include <math.h>
#include <string.h>

/* ---------------------------------------------------------------------------
 * Application 1: Earth's climate energy balance (L7)
 *
 * The Earth-atmosphere system is a non-equilibrium steady state
 * driven by solar radiation. The MEPP predicts the meridional
 * temperature gradient that minimizes entropy production subject
 * to the constraint of total solar energy input.
 *
 * Two-box model: tropical (T_t) and polar (T_p) regions.
 * Heat flux: J = k (T_t - T_p) (Fourier)
 * Entropy production: sigma = J * (1/T_p - 1/T_t)
 *
 * Reference: Paltridge (1975) "Global dynamics and climate"
 *            Lorenz et al. (2001) "Maximum entropy production"
 * ------------------------------------------------------------------------- */

int app_climate_two_box(double solar_input, double albedo,
                        double *T_tropics, double *T_polar, double *sigma)
{
    /* Stefan-Boltzmann constant */
    const double sigma_sb = 5.670374419e-8;  /* [W/(m^2*K^4)] */

    /* Simple energy balance:
     * Tropics: Q_in * (1 - alpha) = sigma_sb * T_t^4 + J
     * Polar:   0                  = sigma_sb * T_p^4 - J
     *
     * where J = k (T_t - T_p) is the meridional heat transport.
     *
     * The MEPP predicts the value of k that maximizes sigma
     * (or equivalently, the T gradient that the circulation
     * organizes itself into).
     */

    double Q_absorbed = solar_input * (1.0 - albedo);

    /* Iterative solution for T_t and T_p given a transport coefficient k */
    double k_heat = 2.0;  /* Effective meridional heat transport coefficient */
    double T_t = 300.0;   /* Initial guess [K] */
    double T_p = 260.0;
    int iter;
    double tol = 1e-6;

    for (iter = 0; iter < 1000; iter++) {
        double J = k_heat * (T_t - T_p);

        double T_t_new = pow((Q_absorbed - J) / sigma_sb, 0.25);
        double T_p_new = pow(J / sigma_sb, 0.25);

        if (fabs(T_t_new - T_t) < tol && fabs(T_p_new - T_p) < tol) {
            T_t = T_t_new;
            T_p = T_p_new;
            break;
        }

        T_t = T_t_new;
        T_p = T_p_new;
    }

    *T_tropics = T_t;
    *T_polar = T_p;

    /* Compute entropy production:
     * sigma = J * (1/T_p - 1/T_t) */
    double J_final = k_heat * (T_t - T_p);
    *sigma = J_final * (1.0 / T_p - 1.0 / T_t);

    if (*sigma < 0.0) *sigma = 0.0;

    return 0;
}

/* ---------------------------------------------------------------------------
 * Application 2: Metabolic network thermodynamics (L7)
 *
 * Living cells maintain non-equilibrium steady states with continuous
 * entropy production. The MEPP has been applied to predict metabolic
 * fluxes in bacteria (E. coli) and mitochondria.
 *
 * Simplified model: glucose metabolism
 *   C6H12O6 + 6 O2 -> 6 CO2 + 6 H2O   (Delta G ~ -2880 kJ/mol)
 *
 * The entropy production is:
 *   sigma = r * (A/T)
 * where r is the reaction rate and A is the affinity.
 *
 * Reference: von Stockar et al. (2006) "Thermodynamics of microbial growth"
 * ------------------------------------------------------------------------- */

int app_metabolic_entropy_production(double glucose_uptake_rate,
                                     double temperature,
                                     double *sigma_metabolic,
                                     double *efficiency)
{
    /* Standard Gibbs free energy of glucose oxidation */
    const double delta_G0_glucose = -2880.0e3;  /* [J/mol] */
    /* Reaction affinity (simplified, ignoring concentration corrections):
     * A = -Delta_G = 2880 kJ/mol */
    double A = -delta_G0_glucose;

    /* Entropy production rate:
     * sigma = r * A / T */
    *sigma_metabolic = glucose_uptake_rate * A / temperature;

    /* Thermodynamic efficiency:
     * eta = (useful work) / (Gibbs free energy available)
     * For glycolysis, about 40% of Gibbs energy is captured as ATP. */
    double atp_per_glucose = 30.0;  /* ~30 ATP per glucose */
    double atp_energy = 50.0e3;     /* ~50 kJ/mol ATP under cellular conditions */
    double useful_work = atp_per_glucose * atp_energy;

    *efficiency = useful_work / (-delta_G0_glucose);

    return 0;
}

/* ---------------------------------------------------------------------------
 * Application 3: Chemical reactor optimization (L7)
 *
 * A continuous stirred-tank reactor (CSTR) with an exothermic reaction
 * operates at a non-equilibrium steady state. The MEPP predicts the
 * optimal operating temperature that balances reaction rate with
 * entropy production from heat transfer.
 *
 * Reaction: A -> B (first-order, exothermic)
 * Heat balance: Q_gen = r * Delta_H_r
 * Cooling: Q_cool = U A_cool (T - T_coolant)
 *
 * Entropy production:
 *   sigma = sigma_reaction + sigma_heat_transfer
 *         = r * (A/T) + U A_cool (T - T_coolant)^2 / (T * T_coolant)
 *
 * Reference: Johannessen & Kjelstrup (2004) "Minimum entropy production
 *            rate in plug flow reactors"
 * ------------------------------------------------------------------------- */

int app_cstr_optimization(double feed_conc, double residence_time,
                          double activation_energy, double pre_exponential,
                          double delta_H_r, double U_A, double T_coolant,
                          double *T_optimal, double *sigma_min)
{
    const double R_gas = 8.314;
    int iter;
    (void)residence_time;
    (void)U_A;
    double T = 350.0;  /* Initial guess [K] */
    double T_min = 300.0, T_max = 600.0;

    /* Scan temperature range to find minimum sigma */
    double best_sigma = 1e100;
    double best_T = 350.0;
    int n_steps = 100;

    for (iter = 0; iter < n_steps; iter++) {
        T = T_min + (T_max - T_min) * iter / (n_steps - 1);

        /* Reaction rate: r = k * C_A = k0 * exp(-Ea/RT) * C_A */
        double k_r = pre_exponential * exp(-activation_energy / (R_gas * T));
        double r = k_r * feed_conc;

        /* Reaction affinity (for A -> B with equimolar):
         * A = -Delta_G = -Delta_H + T * Delta_S_approx */
        double A_react = -delta_H_r;  /* Simplified */

        /* Reaction entropy production */
        double sigma_rxn = r * A_react / T;

        /* Heat transfer entropy production:
         * sigma_Q = Q * (1/T_coolant - 1/T) */
        double Q = r * delta_H_r;  /* Heat generated */
        double sigma_Q = fabs(Q) * fabs(1.0 / T_coolant - 1.0 / T);

        double sigma_total = sigma_rxn + sigma_Q;

        if (sigma_total < best_sigma) {
            best_sigma = sigma_total;
            best_T = T;
        }
    }

    *T_optimal = best_T;
    *sigma_min = best_sigma;

    return 0;
}

/* ---------------------------------------------------------------------------
 * Application 4: Thermoelectric generator efficiency (L7)
 *
 * A thermoelectric generator converts a temperature difference into
 * electrical power through the Seebeck effect. The entropy production
 * has two components: heat conduction and Joule heating.
 *
 * The MEPP can be used to find the optimal load resistance that
 * maximizes power output while minimizing total entropy production.
 *
 * sigma = (k/T_h T_c) (T_h - T_c)^2 + I^2 R / T
 *
 * where k is thermal conductance, I is current, R is resistance.
 *
 * Reference: Ioffe (1957) "Semiconductor Thermoelements and
 *            Thermoelectric Cooling"
 * ------------------------------------------------------------------------- */

int app_thermoelectric_optimization(double T_hot, double T_cold,
                                     double seebeck_coeff, double resistance,
                                     double thermal_conductance,
                                     double *optimal_current,
                                     double *max_power,
                                     double *efficiency)
{
    /* Open-circuit voltage: V_oc = S * (T_h - T_c) */

    double dT = T_hot - T_cold;
    double V_oc = seebeck_coeff * dT;

    /* Power output: P = I * (V_oc - I * R) = I * V_oc - I^2 * R
     * Maximum power at I = V_oc / (2R) (impedance matching) */
    double I_max_power = V_oc / (2.0 * resistance);
    double P_max = I_max_power * V_oc - I_max_power * I_max_power * resistance;

    /* Heat input at hot side (ignoring Joule and Thomson effects):
     * Q_h = k * dT + S * T_h * I */
    double Q_h = thermal_conductance * dT + seebeck_coeff * T_hot * I_max_power;

    /* Efficiency: eta = P / Q_h */
    double eta = (Q_h > 0.0) ? P_max / Q_h : 0.0;

    /* Carnot efficiency for comparison */
    (void)(dT / T_hot); /* eta_carnot - reference value */

    *optimal_current = I_max_power;
    *max_power = P_max;
    *efficiency = eta;

    (void)efficiency;  /* Used via pointer */

    return 0;
}

/* ---------------------------------------------------------------------------
 * Application 5: Rayleigh-Benard convection onset (L7)
 *
 * Rayleigh-Benard convection is the canonical example of a
 * dissipative structure: a fluid layer heated from below develops
 * convection rolls when the temperature gradient exceeds a
 * critical value (Rayleigh number Ra > Ra_c).
 *
 * The MEPP predicts that the convective state has HIGHER entropy
 * production than the conductive state at the same boundary
 * conditions, consistent with the idea that systems far from
 * equilibrium select states of higher dissipation.
 *
 * sigma_conductive = k (Delta T)^2 / (T_avg^2 d)
 * sigma_convective = sigma_conductive * Nu
 *
 * where Nu > 1 is the Nusselt number (enhanced heat transport).
 *
 * Reference: Chandrasekhar (1961) "Hydrodynamic and Hydromagnetic Stability"
 *            Cross & Hohenberg (1993) Rev. Mod. Phys. 65, 851
 * ------------------------------------------------------------------------- */

int app_rayleigh_benard_entropy(double delta_T, double layer_depth,
                                 double thermal_cond, double viscosity,
                                 double thermal_diffusivity, double expansion_coeff,
                                 double *Ra, double *Ra_critical,
                                 double *sigma_conductive, double *sigma_convective)
{
    /* Gravitational acceleration */
    const double g = 9.81;  /* [m/s^2] */

    /* Rayleigh number:
     * Ra = g * alpha * Delta_T * d^3 / (nu * kappa)
     * where alpha = thermal expansion coefficient
     *       nu = kinematic viscosity
     *       kappa = thermal diffusivity
     *       d = layer depth
     *       Delta_T = temperature difference */
    double Ra_val = g * expansion_coeff * delta_T
                   * layer_depth * layer_depth * layer_depth
                   / (viscosity * thermal_diffusivity);

    /* Critical Rayleigh number for rigid-rigid boundaries:
     * Ra_c ~ 1708 */
    double Ra_c = 1707.76;

    *Ra = Ra_val;
    *Ra_critical = Ra_c;

    /* Average temperature */
    double T_avg = 300.0;  /* Assumed average [K] */

    /* Conductive entropy production:
     * sigma_cond = k (Delta T)^2 / (T_avg^2 * d) */
    *sigma_conductive = thermal_cond * delta_T * delta_T
                      / (T_avg * T_avg * layer_depth);

    /* Convective enhancement: Nu ~ (Ra/Ra_c)^(1/3) for Ra > Ra_c */
    double Nu;
    if (Ra_val > Ra_c) {
        Nu = pow(Ra_val / Ra_c, 1.0 / 3.0);
        *sigma_convective = *sigma_conductive * Nu;
    } else {
        Nu = 1.0;
        *sigma_convective = *sigma_conductive;
    }

    return 0;
}

/* ---------------------------------------------------------------------------
 * Application summary printer
 * ------------------------------------------------------------------------- */

void app_print_summary(void)
{
    printf("============================================\n");
    printf("  Entropy Production Minimization\n");
    printf("  Application Examples Summary\n");
    printf("============================================\n\n");

    /* Climate example */
    double T_t, T_p, sigma_climate;
    app_climate_two_box(342.0, 0.3, &T_t, &T_p, &sigma_climate);
    printf("1. Climate Two-Box Model:\n");
    printf("   T_tropics = %.1f K, T_polar = %.1f K\n", T_t, T_p);
    printf("   sigma = %.6f W/(m^2*K)\n\n", sigma_climate);

    /* Metabolic example */
    double sigma_met, efficiency;
    app_metabolic_entropy_production(1.0e-3, 310.0, &sigma_met, &efficiency);
    printf("2. Metabolic Network (glucose oxidation):\n");
    printf("   sigma = %.3f W/K per mol/s glucose\n", sigma_met);
    printf("   ATP capture efficiency = %.1f%%\n\n", efficiency * 100.0);

    /* CSTR example */
    double T_opt, sigma_min;
    app_cstr_optimization(1.0, 3600.0, 80000.0, 1.0e10,
                           -200000.0, 500.0, 300.0,
                           &T_opt, &sigma_min);
    printf("3. CSTR Reactor Optimization:\n");
    printf("   Optimal temperature = %.1f K\n", T_opt);
    printf("   Minimum sigma = %.3f W/K\n\n", sigma_min);

    /* Thermoelectric example */
    double I_opt, P_max, eta_te;
    app_thermoelectric_optimization(500.0, 300.0, 200.0e-6,
                                     0.01, 1.0,
                                     &I_opt, &P_max, &eta_te);
    printf("4. Thermoelectric Generator:\n");
    printf("   Optimal current = %.6f A\n", I_opt);
    printf("   Max power = %.6f W\n", P_max);
    printf("   Efficiency = %.2f%%\n\n", eta_te * 100.0);

    /* Rayleigh-Benard example */
    double Ra, Ra_c, sigma_cond, sigma_conv;
    app_rayleigh_benard_entropy(20.0, 0.01, 0.6,
                                 1.0e-6, 1.4e-7, 2.0e-4,
                                 &Ra, &Ra_c, &sigma_cond, &sigma_conv);
    printf("5. Rayleigh-Benard Convection:\n");
    printf("   Ra = %.1f (Ra_c = %.1f)\n", Ra, Ra_c);
    printf("   sigma_conductive = %.6f W/K\n", sigma_cond);
    printf("   sigma_convective = %.6f W/K\n", sigma_conv);
    printf("   Entropy production ratio = %.2f\n\n",
           sigma_conv / (sigma_cond + 1e-15));
}