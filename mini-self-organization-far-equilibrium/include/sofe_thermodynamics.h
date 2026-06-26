#ifndef SOFE_THERMODYNAMICS_H
#define SOFE_THERMODYNAMICS_H

#include "sofe_core.h"

typedef struct {
    double entropy;
    double internal_energy;
    double helmholtz_free_energy;
    double gibbs_free_energy;
    double enthalpy;
    double temperature;
    double pressure;
    double volume;
    double chemical_potential;
    double particle_number;
} SOFEThermoState;

typedef struct {
    double* L;
    double* thermodynamic_forces;
    double* thermodynamic_fluxes;
    int n_forces;
    bool is_symmetric;
    bool magnetic_field_present;
} SOFEOnsagerMatrix;

typedef struct {
    double total_production;
    double entropy_flow;
    double internal_production;
    double thermal_conduction_term;
    double diffusion_term;
    double chemical_reaction_term;
    double viscous_term;
    double electrical_term;
    double radiation_term;
} SOFEEntropyProduction;

typedef struct {
    char* name;
    double force;
    double flux;
    double phenomenological_coeff;
    double* cross_coefficients;
    int n_cross;
    int conjugate_flux_index;
} SOFEForceFluxPair;

typedef struct {
    SOFEForceFluxPair* pairs;
    int n_pairs;
    int capacity;
} SOFEForceFluxSystem;

typedef struct {
    double excess_entropy_production;
    double excess_entropy;
    double time_derivative_excess_entropy;
    bool is_stable;
    double critical_threshold;
    double stability_margin;
} SOFEGlansdorffPrigogineCriterion;

typedef struct {
    double* correlation_function;
    double* response_function;
    int n_timepoints;
    double dt;
    double temperature;
    double fluctuation_dissipation_ratio;
    bool fdt_holds;
} SOFEFluctuationDissipation;

double sofe_entropy_production_direct(SOFEForceFluxSystem* ffs);
double sofe_entropy_production_reactions(SOFEReactionNetwork* net, double T);
double sofe_entropy_production_master_equation(double* p, double** W,
                                                int n_states, double k_B);
double sofe_entropy_flow(double T, double heat_flow,
                          double* chem_potentials, double* particle_flows,
                          int n_species);

int sofe_onsager_create(SOFEOnsagerMatrix* ons, int n_forces);
void sofe_onsager_free(SOFEOnsagerMatrix* ons);
void sofe_onsager_set_coefficient(SOFEOnsagerMatrix* ons, int i, int j,
                                   double value);
double sofe_onsager_get_coefficient(SOFEOnsagerMatrix* ons, int i, int j);
bool sofe_onsager_verify_reciprocity(SOFEOnsagerMatrix* ons, double tolerance);
void sofe_onsager_compute_fluxes(SOFEOnsagerMatrix* ons);
double sofe_rayleigh_onsager_dissipation(SOFEOnsagerMatrix* ons);

void sofe_gp_init(SOFEGlansdorffPrigogineCriterion* gp);
void sofe_gp_compute(SOFEGlansdorffPrigogineCriterion* gp,
                      SOFEReactionNetwork* net, double T,
                      double* perturbation, int n_species);
bool sofe_gp_is_stable(SOFEGlansdorffPrigogineCriterion* gp);
double sofe_excess_entropy(double* state_perturbation,
                            double* heat_capacities, int n);

void sofe_fd_init(SOFEFluctuationDissipation* fd, int n_timepoints,
                   double dt, double T);
void sofe_fd_free(SOFEFluctuationDissipation* fd);
void sofe_fd_compute_correlation(SOFEFluctuationDissipation* fd,
                                  double* timeseries, int n_samples);
void sofe_fd_compute_response(SOFEFluctuationDissipation* fd,
                               double* perturbed, double* unperturbed,
                               int n_samples, double perturbation_strength);
bool sofe_fd_verify(SOFEFluctuationDissipation* fd, double tolerance);
double sofe_einstein_fluctuation_probability(double delta_entropy,
                                              double k_B);

double sofe_crooks_ratio(double work, double free_energy_diff, double kT);
double sofe_jarzynski_estimator(double* work_values, int n_trajectories,
                                 double kT);

double sofe_thermodynamic_uncertainty_ratio(double current_mean,
                                              double current_variance,
                                              double entropy_prod_rate,
                                              double k_B);

double sofe_housekeeping_heat(double* steady_state_prob,
                               double* transition_rates,
                               int n_states, double kT);
double sofe_excess_heat(double total_heat, double housekeeping_heat);

void sofe_thermo_state_print(SOFEThermoState* ts);
double sofe_boltzmann_entropy(double multiplicity, double k_B);
double sofe_gibbs_entropy(double* probabilities, int n_states, double k_B);
double sofe_ideal_gas_chemical_potential(double T, double P,
                                           double P_ref, double mu_ref);

#endif
