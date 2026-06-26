/**
 * @file entropy_production.h
 * @brief Core definitions for entropy production in non-equilibrium thermodynamics.
 *
 * Defines the fundamental quantities: entropy production rate sigma,
 * thermodynamic forces X_i, fluxes J_i, and the bilinear form
 * sigma = sum_i J_i * X_i >= 0 (Second Law).
 *
 * References:
 *   - Prigogine (1947) "Etude thermodynamique des phenomenes irreversibles"
 *   - de Groot & Mazur (1984) "Non-Equilibrium Thermodynamics"
 *   - Kondepudi & Prigogine (2014) "Modern Thermodynamics"
 *   - Onsager (1931) "Reciprocal Relations in Irreversible Processes I, II"
 */

#ifndef ENTROPY_PRODUCTION_H
#define ENTROPY_PRODUCTION_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ---------------------------------------------------------------------------
 * L1: Core Definitions -- Fundamental thermodynamic quantities
 * ------------------------------------------------------------------------- */

/** Maximum number of independent thermodynamic forces/fluxes in the system */
#define EP_MAX_FORCES 64

/** Maximum number of chemical species in a reaction network */
#define EP_MAX_SPECIES 128

/** Maximum number of chemical reactions */
#define EP_MAX_REACTIONS 256

/**
 * @brief Thermodynamic state descriptor.
 *
 * Contains the extensive and intensive variables that define the local
 * thermodynamic state under the local equilibrium hypothesis.  This is
 * the foundational structure for all entropy production computations.
 *
 * L1: Definition of thermodynamic state -- the starting point for any
 *     non-equilibrium thermodynamic analysis.
 */
typedef struct {
    double temperature;           /**< Absolute temperature [K] */
    double pressure;              /**< Pressure [Pa] */
    double volume;                /**< Volume [m^3] */
    double internal_energy;       /**< Internal energy U [J] */
    double entropy;               /**< Entropy S [J/K] */
    double helmholtz_free_energy; /**< Helmholtz free energy F = U - TS [J] */
    double gibbs_free_energy;     /**< Gibbs free energy G = H - TS [J] */
    double enthalpy;              /**< Enthalpy H = U + PV [J] */
    double chemical_potential[EP_MAX_SPECIES]; /**< mu_k for each species [J/mol] */
    double mole_numbers[EP_MAX_SPECIES];       /**< N_k [mol] */
    double concentrations[EP_MAX_SPECIES];     /**< c_k = N_k/V [mol/m^3] */
    int    num_species;           /**< Number of chemical species present */
} ThermodynamicState;

/**
 * @brief Entropy production rate descriptor.
 *
 * The total entropy production rate sigma = d_iS/dt, where d_iS is the
 * entropy produced internally by irreversible processes.
 *
 * For a continuous system with local equilibrium:
 *   sigma = integral_V (sum_i J_i * X_i) dV  >= 0
 *
 * The component-wise decomposition allows identifying which irreversible
 * processes contribute most to the total entropy production.
 */
typedef struct {
    double sigma_total;           /**< Total entropy production rate [W/K] */
    double sigma_per_volume;      /**< Volumetric entropy production [W/(K*m^3)] */
    double sigma_components[EP_MAX_FORCES]; /**< Component sigma_i = J_i * X_i */
    double sigma_heat;            /**< Entropy production from heat conduction */
    double sigma_diffusion;       /**< Entropy production from mass diffusion */
    double sigma_reaction;        /**< Entropy production from chemical reactions */
    double sigma_viscous;         /**< Entropy production from viscous dissipation */
    double sigma_electrical;      /**< Entropy production from electrical conduction */
    int    num_components;        /**< Number of active irreversible processes */
    int    is_positive;           /**< 1 if sigma >= 0 (Second Law check passed) */
    double second_law_margin;     /**< sigma / |sigma_max_component| safety margin */
} EntropyProduction;

/**
 * @brief Thermodynamic flux descriptor.
 *
 * Fluxes J_i are the rates of irreversible processes.  They are the
 * "responses" to thermodynamic forces in the linear regime.
 *
 * L1: Definition of thermodynamic fluxes -- the kinetic quantities that,
 *     together with forces, determine entropy production.
 */
typedef struct {
    double heat_flux;             /**< Heat flux J_q [W/m^2] */
    double mass_flux[EP_MAX_SPECIES]; /**< Diffusion flux J_k [mol/(m^2*s)] */
    double reaction_rate[EP_MAX_REACTIONS]; /**< Reaction rate r_rho [mol/(m^3*s)] */
    double momentum_flux_tensor[3][3];   /**< Viscous pressure tensor Pi [Pa] */
    double electric_current;      /**< Electric current density i [A/m^2] */
    int    num_species;           /**< Number of diffusing species */
    int    num_reactions;         /**< Number of chemical reactions */
} ThermodynamicFlux;

/**
 * @brief Thermodynamic force (affinity) descriptor.
 *
 * Forces X_i are gradients of intensive thermodynamic variables that
 * drive irreversible processes.  Each force is conjugate to a flux.
 *
 * L1: Definition of thermodynamic forces -- the "causes" that drive
 *     irreversible processes away from equilibrium.
 *
 * Force definitions following de Groot & Mazur, Ch. IV:
 *   X_q  = grad(1/T)             thermal force
 *   X_k  = -grad(mu_k/T)         diffusion force for species k
 *   A_rho = -sum_k nu_{k,rho} mu_k / T   chemical affinity for reaction rho
 *   X_v  = -grad(v)/T            viscous force (tensor)
 *   X_e  = -grad(phi)/T          electrical force
 */
typedef struct {
    double thermal_force;         /**< grad(1/T) thermal driving force [K^{-1}*m^{-1}] */
    double thermal_force_grad[3]; /**< Vector components of grad(1/T) */
    double diffusion_force[EP_MAX_SPECIES]; /**< -grad(mu_k/T) for each species */
    double diffusion_force_grad[EP_MAX_SPECIES][3]; /**< Vector diffusion forces */
    double affinity[EP_MAX_REACTIONS];      /**< -sum_k nu_{k,rho} mu_k / T [J/(K*mol)] */
    double viscous_force[3][3];   /**< -grad(v)/T [s^{-1}/K] */
    double electrical_force;      /**< -grad(phi)/T [V/(K*m)] */
    double electrical_force_grad[3]; /**< Vector of -grad(phi)/T */
    int    num_species;
    int    num_reactions;
} ThermodynamicForce;

/**
 * @brief Chemical reaction descriptor with stoichiometry.
 *
 * L1: Definition of chemical reaction in thermodynamic context.
 * Each reaction rho: 0 = sum_k nu_{k,rho} M_k
 * where nu_{k,rho} > 0 for products, < 0 for reactants.
 */
typedef struct {
    double stoichiometry[EP_MAX_SPECIES]; /**< Stoichiometric coefficients nu_{k,rho} */
    double forward_rate_constant;    /**< k_f [appropriate units] */
    double reverse_rate_constant;    /**< k_r [appropriate units] */
    double activation_energy;        /**< E_a [J/mol] -- Arrhenius kinetics */
    double reaction_enthalpy;        /**< Delta H_r [J/mol] */
    int    num_species_involved;     /**< How many species participate */
    int    is_elementary;            /**< 1 if elementary step, 0 if overall */
} ChemicalReaction;

/* ---------------------------------------------------------------------------
 * L2: Core Concepts -- Entropy production rate computation
 * ------------------------------------------------------------------------- */

/**
 * @brief Compute entropy production rate from force-flux pairs.
 *
 * The fundamental bilinear form of irreversible thermodynamics:
 *   sigma = sum_i J_i * X_i
 *
 * Each term J_i * X_i represents the entropy production contribution
 * from a specific irreversible process i.  The sum must be non-negative
 * for any physically admissible combination of forces and fluxes.
 *
 * @param forces    Array of thermodynamic forces (length n)
 * @param fluxes    Array of thermodynamic fluxes (length n)
 * @param n         Number of force-flux pairs
 * @return          Entropy production rate sigma [W/K]
 *
 * Complexity: O(n)
 * Theorem: Second Law of Thermodynamics requires sigma >= 0.
 * Reference: de Groot & Mazur (1984), Ch. III, Eq. (21)
 */
double ep_compute_sigma(const double forces[], const double fluxes[], int n);

/**
 * @brief Compute volumetric entropy production rate for continuous media.
 *
 * sigma(x) = J_q * grad(1/T)
 *          + sum_k J_k * (-grad(mu_k/T))
 *          + sum_rho r_rho * (A_rho/T)
 *          + Pi : grad(v)/T
 *          + i * (-grad(phi)/T)
 *
 * Each term corresponds to a distinct irreversible process:
 *   - Heat conduction
 *   - Multicomponent diffusion
 *   - Chemical reactions
 *   - Viscous flow
 *   - Electrical conduction
 *
 * @param flux       Thermodynamic flux structure
 * @param force      Thermodynamic force structure
 * @param temperature Local temperature [K]
 * @return           Volumetric entropy production rate [W/(K*m^3)]
 *
 * Complexity: O(num_species + num_reactions)
 */
double ep_volumetric_sigma(const ThermodynamicFlux *flux,
                           const ThermodynamicForce *force,
                           double temperature);

/**
 * @brief Compute the total entropy change dS = d_eS + d_iS.
 *
 * The entropy balance equation:
 *   dS = d_eS + d_iS
 *
 * where:
 *   d_eS = dQ/T + sum_k s_k * dN_k   (entropy exchange with surroundings)
 *   d_iS >= 0                         (internal entropy production)
 *
 * @param state       Current thermodynamic state
 * @param heat_flow   Heat received from surroundings dQ [J]
 * @param mass_flow   Mass exchanged with surroundings dN_k [mol]
 * @param work_done   Work done on system dW [J]
 * @param d_eS        [out] Entropy exchanged with surroundings [J/K]
 * @param d_iS        [out] Entropy produced internally [J/K]
 * @return            0 on success, -1 on invalid state pointer
 */
int ep_entropy_balance(const ThermodynamicState *state,
                       double heat_flow,
                       const double mass_flow[],
                       double work_done,
                       double *d_eS,
                       double *d_iS);

/**
 * @brief Verify the Second Law: sigma >= 0.
 *
 * @param sigma     Computed entropy production rate
 * @param epsilon   Numerical tolerance (suggested: 1e-12)
 * @return          1 if sigma >= -epsilon, 0 otherwise
 */
int ep_check_second_law(double sigma, double epsilon);

/* ---------------------------------------------------------------------------
 * L3: Mathematical Structures -- Bilinear and quadratic forms
 * ------------------------------------------------------------------------- */

/**
 * @brief The Rayleigh-Onsager dissipative function.
 *
 * Phi(J, J) = (1/2) sum_{i,j} R_{ij} J_i J_j
 *
 * where R = L^{-1} is the resistance matrix. The Rayleigh dissipative
 * function is half the entropy production when expressed in fluxes:
 *   sigma = 2 * Phi(J, J)
 *
 * @param fluxes     Array of fluxes (length n)
 * @param R          Resistance matrix (n x n, row-major, symmetric)
 * @param n          System dimension
 * @return           Rayleigh dissipative function value Phi
 *
 * Reference: Rayleigh (1873), Onsager (1931)
 */
double ep_rayleigh_dissipative(const double fluxes[], const double *R, int n);

/**
 * @brief Entropy production as quadratic form in forces.
 *
 * sigma = X^T L X = sum_{i,j} L_{ij} X_i X_j
 *
 * @param forces     Array of forces (length n)
 * @param L          Onsager matrix (n x n, row-major, symmetric)
 * @param n          System dimension
 * @return           Entropy production rate sigma [W/K]
 */
double ep_quadratic_form_forces(const double forces[], const double *L, int n);

/**
 * @brief Entropy production as quadratic form in fluxes.
 *
 * sigma = J^T R J = sum_{i,j} R_{ij} J_i J_j = 2 * Phi(J, J)
 *
 * @param fluxes     Array of fluxes (length n)
 * @param R          Resistance matrix (n x n, row-major)
 * @param n          System dimension
 * @return           Entropy production rate sigma [W/K]
 */
double ep_quadratic_form_fluxes(const double fluxes[], const double *R, int n);

/**
 * @brief Compute the bilinear coupling between two sets of forces.
 *
 * sigma_{12} = X1^T L X2
 *
 * @param forces1    First force vector
 * @param forces2    Second force vector
 * @param L          Onsager coupling matrix
 * @param n          Dimension
 * @return           Cross entropy production contribution
 */
double ep_cross_entropy_production(const double forces1[],
                                   const double forces2[],
                                   const double *L, int n);

/* ---------------------------------------------------------------------------
 * L4: Fundamental Laws -- Prigogine Minimum Entropy Production Theorem
 * ------------------------------------------------------------------------- */

/**
 * @brief Compute d(sigma)/dt for Prigogine's minimum entropy production theorem.
 *
 * Prigogine's theorem: for linear constitutive relations with time-independent
 * boundary conditions, d(sigma)/dt <= 0 in the approach to steady state.
 *
 * d(sigma)/dt = 2 sum_{i,j} L_{ij} X_i (dX_j/dt)
 *
 * @param forces       Current force values (length n)
 * @param dforces_dt   Time derivatives of forces (length n)
 * @param L            Onsager coefficient matrix (n x n, row-major)
 * @param n            Dimension
 * @return             d(sigma)/dt [W/(K*s)]
 *
 * Reference: Prigogine (1945) Acad. Roy. Belg., Bull. Classe Sci., 31, 600-606
 */
double ep_prigogine_dsigma_dt(const double forces[],
                              const double dforces_dt[],
                              const double *L, int n);

/**
 * @brief Determine if a state satisfies the MEP steady-state conditions.
 *
 * Checks partial(sigma)/partial(X_k) = 0 for all unconstrained forces.
 *
 * @param forces        Force values at candidate steady state
 * @param L             Onsager matrix (n x n)
 * @param constrained   constrained[i]=1 if force X_i is fixed
 * @param n             Dimension
 * @param tol           Convergence tolerance
 * @return              1 if MEP steady state, 0 otherwise
 */
int ep_is_mep_steady_state(const double forces[],
                           const double *L,
                           const int constrained[],
                           int n, double tol);

/**
 * @brief Find the MEP steady state by constrained optimization.
 *
 * Minimizes sigma(X) = X^T L X subject to X_k = fixed_values[k] for
 * constrained forces.
 *
 * @param forces        [in/out] Initial guess, overwritten with solution
 * @param L             Onsager matrix (n x n)
 * @param constrained   constrained[i]=1 if force X_i is fixed
 * @param fixed_values  Fixed values for constrained forces
 * @param n             Dimension
 * @param max_iter      Maximum iterations
 * @param tol           Convergence tolerance
 * @return              0 on convergence, -1 on failure
 */
int ep_find_mep_steady_state(double forces[],
                             const double *L,
                             const int constrained[],
                             const double fixed_values[],
                             int n, int max_iter, double tol);

/* ---------------------------------------------------------------------------
 * Utility functions
 * ------------------------------------------------------------------------- */

void ep_state_init(ThermodynamicState *state, double temperature,
                   double pressure, int num_species);

void ep_state_copy(ThermodynamicState *dest, const ThermodynamicState *src);

double ep_ideal_gas_entropy(const ThermodynamicState *state,
                            const double std_entropy[]);

void ep_ideal_gas_chemical_potentials(ThermodynamicState *state,
                                      const double std_potential[],
                                      double R);

double ep_reaction_affinity(const ChemicalReaction *reaction,
                            const double chem_pot[], int num_species);

#ifdef __cplusplus
}
#endif

#endif /* ENTROPY_PRODUCTION_H */
