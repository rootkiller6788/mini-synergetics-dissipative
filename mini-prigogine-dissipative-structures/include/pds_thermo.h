#ifndef PDS_THERMO_H
#define PDS_THERMO_H

#include "pds_core.h"

/* ============================================================================
 * PDS Thermodynamics - Entropy Production, Onsager, GP Stability
 *
 * Key Theorems:
 *   L4: Prigogine Minimum Entropy Production Theorem (1945)
 *   L4: Glansdorff-Prigogine Stability Criterion (1954, 1971)
 *   L4: Onsager Reciprocal Relations (1931)
 *   L4: Second Law for Open Systems (dS = d_i_S + d_e_S, d_i_S >= 0)
 *   L8: Extended Irreversible Thermodynamics (EIT)
 * ============================================================================ */

/* --- L1: Thermodynamic Force and Flux Types --- */

typedef enum {
    PDS_FORCE_TEMPERATURE_GRAD = 0,
    PDS_FORCE_CHEMICAL_AFFINITY = 1,
    PDS_FORCE_CONCENTRATION_GRAD = 2,
    PDS_FORCE_VELOCITY_GRAD = 3,
    PDS_FORCE_ELECTRIC_FIELD = 4,
    PDS_FORCE_GRAVITATIONAL = 5
} PDSThermoForceType;

typedef enum {
    PDS_FLUX_HEAT = 0,
    PDS_FLUX_DIFFUSION = 1,
    PDS_FLUX_CHEMICAL = 2,
    PDS_FLUX_MOMENTUM = 3,
    PDS_FLUX_ELECTRIC = 4
} PDSThermoFluxType;

/* --- L3: Thermodynamic Force-Flux Pair --- */

typedef struct {
    PDSThermoForceType force_type;
    PDSThermoFluxType flux_type;
    double force_value;
    double flux_value;
    double phenomenological_coeff;
    double entropy_production_contribution;
} PDSForceFluxPair;

/* --- L3: Onsager Coefficient Matrix --- */

typedef struct {
    double** L;
    int n_pairs;
    PDSThermoForceType* force_types;
    PDSThermoFluxType* flux_types;
    bool is_symmetric;
    bool is_positive_definite;
} PDSOnsagerMatrix;

/* --- L3: Local Equilibrium Cell --- */

typedef struct {
    double temperature;
    double pressure;
    double* chemical_potentials;
    double* concentrations;
    int n_species;
    double local_entropy;
    double local_entropy_production;
    double local_free_energy;
    bool in_local_equilibrium;
} PDSLocalEquilibriumCell;

/* --- L3: Entropy Production Profile (spatial) --- */

typedef struct {
    PDSLocalEquilibriumCell* cells;
    int n_cells;
    double total_entropy_production;
    double mean_entropy_production;
    double max_entropy_production;
    int max_production_cell_index;
    double* entropy_production_profile;
    double* entropy_flow_profile;
} PDSEntropyProductionProfile;

/* --- L4: Minimum Entropy Production State --- */

typedef struct {
    bool is_minimum;
    bool is_steady_state;
    double entropy_production;
    double dP_dt;
    double minimum_possible_P;
    double distance_from_minimum;
    double* gradient_dP_dlambda;
    int n_gradient_components;
} PDSMinEntropyProduction;

/* --- L8: Extended Irreversible Thermodynamics State --- */

typedef struct {
    double* extended_fluxes;
    double* relaxation_times;
    int n_extended_vars;
    double generalized_entropy;
    double generalized_entropy_production;
    bool hyperbolic_heat;
} PDSExtendedThermo;

/* --- L4: Glansdorff-Prigogine Stability --- */

typedef struct {
    double excess_entropy;
    double excess_entropy_production;
    bool is_stable_by_gp;
    bool is_lyapunov_function;
    double gp_stability_margin;
    double** jacobian_of_forces;
    double** jacobian_of_fluxes;
    int n_variables;
} PDSGlansdorffPrigogineStability;

/* ============================================================================
 * API Declarations
 * ============================================================================ */

void pds_entropy_balance_update(PDSDissipativeSystem* sys, double dt,
                                 double env_temperature,
                                 const double* env_chem_potentials);
double pds_entropy_production_rate(const PDSDissipativeSystem* sys,
                                    double temperature);
double pds_chemical_affinity(const PDSDissipativeSystem* sys, int reaction_idx,
                              double temperature, const double* chem_potentials);
PDSOnsagerMatrix* pds_onsager_matrix_create(int n_pairs);
void pds_onsager_matrix_free(PDSOnsagerMatrix* om);
void pds_onsager_matrix_set(PDSOnsagerMatrix* om, int i, int j, double value);
void pds_onsager_matrix_check_symmetry(PDSOnsagerMatrix* om);
bool pds_onsager_matrix_is_positive_definite(const PDSOnsagerMatrix* om);
void pds_onsager_matrix_compute_fluxes(const PDSOnsagerMatrix* om,
                                        const double* forces, double* fluxes);
bool pds_check_second_law(const PDSDissipativeSystem* sys, double tolerance);
PDSMinEntropyProduction* pds_minimum_entropy_production_check(
    const PDSDissipativeSystem* sys);
void pds_min_entropy_production_free(PDSMinEntropyProduction* mep);
PDSGlansdorffPrigogineStability* pds_glansdorff_prigogine_stability(
    const PDSDissipativeSystem* sys);
void pds_gp_stability_free(PDSGlansdorffPrigogineStability* gps);
bool pds_local_equilibrium_check(const PDSDissipativeSystem* sys,
                                  double max_gradient_norm);
PDSEntropyProductionProfile* pds_entropy_production_profile_compute(
    const PDSDissipativeSystem* sys, double temperature);
void pds_entropy_profile_free(PDSEntropyProductionProfile* epp);
PDSExtendedThermo* pds_extended_thermo_create(int n_extended_vars);
void pds_extended_thermo_free(PDSExtendedThermo* eit);
void pds_extended_thermo_evolve(PDSExtendedThermo* eit, double dt);

#endif /* PDS_THERMO_H */