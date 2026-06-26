/**
 * @file dissipative_systems.h
 * @brief Dissipative structures, reaction-diffusion systems, and self-organization.
 *
 * Dissipative structures are organized states of matter that arise
 * far from thermodynamic equilibrium through the dissipation of energy
 * and matter. They are maintained by continuous entropy production
 * and exchange with the environment.
 *
 * Key concepts:
 *   - Turing instability (reaction-diffusion pattern formation)
 *   - Brusselator and Oregonator models
 *   - Bifurcation analysis from equilibrium
 *   - Order through fluctuations (Prigogine)
 *   - Synergetics (Haken): slaving principle, order parameters
 *
 * References:
 *   - Turing, A.M. (1952) "The Chemical Basis of Morphogenesis"
 *   - Prigogine & Lefever (1968) "Symmetry Breaking Instabilities in
 *     Dissipative Systems", J. Chem. Phys. 48, 1695
 *   - Haken, H. (1977) "Synergetics: An Introduction"
 *   - Nicolis & Prigogine (1977) "Self-Organization in Nonequilibrium Systems"
 *   - Cross & Hohenberg (1993) "Pattern formation outside of equilibrium"
 */

#ifndef DISSIPATIVE_SYSTEMS_H
#define DISSIPATIVE_SYSTEMS_H

#include "entropy_production.h"
#include "onsager_relations.h"

#ifdef __cplusplus
extern "C" {
#endif

#define DS_MAX_SPECIES 16
#define DS_MAX_GRID 256

/* ---------------------------------------------------------------------------
 * L1: Definitions -- Reaction-diffusion system
 * ------------------------------------------------------------------------- */

/**
 * @brief Reaction-diffusion system descriptor.
 *
 * A reaction-diffusion system describes the spatiotemporal evolution
 * of chemical concentrations c_k(x, t):
 *
 *   dc_k/dt = D_k grad^2 c_k + R_k({c_j})
 *
 * where D_k are diffusion coefficients and R_k are reaction rates.
 * Such systems can spontaneously form spatial patterns (Turing patterns)
 * when the homogeneous steady state becomes unstable to spatial
 * perturbations.
 */
typedef struct {
    int     num_species;                  /**< Number of chemical species */
    double  diffusion_coeff[DS_MAX_SPECIES]; /**< D_k [m^2/s] */
    double  concentrations[DS_MAX_SPECIES]; /**< Current concentrations [mol/m^3] */
    double  reaction_params[DS_MAX_SPECIES * 4]; /**< Reaction rate parameters */
    double  temperature;                  /**< System temperature [K] */
    int     grid_points;                  /**< Number of spatial grid points */
    double  domain_length;                /**< Physical domain length [m] */
    double  dx;                           /**< Grid spacing [m] */
} ReactionDiffusionSystem;

/**
 * @brief Spatial pattern descriptor for dissipative structures.
 *
 * Captures the spatial organization that emerges from a
 * reaction-diffusion instability. Patterns are characterized
 * by their dominant wavenumber, amplitude, and symmetry.
 */
typedef struct {
    int     num_species;
    int     grid_points;
    double *pattern[DS_MAX_SPECIES];      /**< Concentration field for each species */
    double  dominant_wavenumber;          /**< k_c, most unstable wavenumber [1/m] */
    double  amplitude;                    /**< Pattern amplitude */
    double  wavelength;                   /**< lambda = 2*pi/k_c [m] */
    int     symmetry_type;               /**< 0=homogeneous, 1=stripes, 2=hexagons, 3=spirals */
    double  entropy_production;           /**< sigma of the patterned state [W/K] */
    int     is_stable;                    /**< 1 if pattern is stable */
} DissipativePattern;

/**
 * @brief Bifurcation analysis descriptor.
 *
 * As a control parameter (e.g., temperature gradient, concentration
 * difference) is varied, the system may undergo bifurcations where
 * new steady states emerge or existing ones lose stability.
 */
typedef struct {
    double  control_parameter;            /**< Bifurcation parameter value */
    double  critical_value;               /**< Critical value where bifurcation occurs */
    double  eigenvalues[DS_MAX_SPECIES];   /**< Linear stability eigenvalues */
    int     num_modes;                    /**< Number of eigenmodes analyzed */
    int     bifurcation_type;            /**< 0=pitchfork, 1=Hopf, 2=saddle-node, 3=transcritical */
    double  order_parameter;             /**< Amplitude of the emerging mode */
    int     is_supercritical;            /**< 1 if supercritical (continuous), 0 if subcritical */
} BifurcationAnalysis;

/* ---------------------------------------------------------------------------
 * L2: Core Concepts -- Turing instability and pattern formation
 * ------------------------------------------------------------------------- */

/**
 * @brief Compute the homogeneous steady state of a reaction-diffusion system.
 *
 * Finds concentrations c_k^* such that R_k({c_j^*}) = 0 for all k,
 * i.e., the reaction rates balance at the steady state.
 *
 * Uses Newton-Raphson iteration with line search for robustness.
 *
 * @param rds           Reaction-diffusion system [in/out: concentrations updated]
 * @param max_iter      Maximum Newton iterations
 * @param tol           Convergence tolerance
 * @return              0 on success, -1 if no steady state found
 */
int ds_find_homogeneous_steady_state(ReactionDiffusionSystem *rds,
                                     int max_iter, double tol);

/**
 * @brief Compute the Turing instability condition.
 *
 * A homogeneous steady state is Turing-unstable if there exists a
 * wavenumber k > 0 such that the Jacobian J - k^2 D has an eigenvalue
 * with positive real part.
 *
 * The Turing conditions for a 2-species system:
 *   J_{11} + J_{22} < 0            (stable without diffusion)
 *   J_{11} J_{22} - J_{12} J_{21} > 0
 *   D_2 J_{11} + D_1 J_{22} > 0    (diffusion-driven instability)
 *   (D_2 J_{11} + D_1 J_{22})^2 >= 4 D_1 D_2 (J_{11} J_{22} - J_{12} J_{21})
 *
 * @param rds       Reaction-diffusion system at homogeneous steady state
 * @param jacobian  Jacobian matrix of reaction terms (n x n, row-major)
 * @param k_c       [out] Critical wavenumber of most unstable mode
 * @param lambda_max [out] Maximum growth rate at k_c
 * @return          1 if Turing unstable, 0 if stable
 *
 * Reference: Turing (1952), Murray (2003) "Mathematical Biology"
 */
int ds_check_turing_instability(const ReactionDiffusionSystem *rds,
                                const double *jacobian,
                                double *k_c, double *lambda_max);

/**
 * @brief Compute the dispersion relation lambda(k) for a reaction-diffusion system.
 *
 * The dispersion relation lambda(k) gives the growth rate of a
 * spatial perturbation with wavenumber k:
 *
 *   lambda(k) = eigenvalues of (J - k^2 D)
 *
 * where J is the Jacobian of the reaction terms and D is the
 * diagonal diffusion matrix.
 *
 * @param rds       Reaction-diffusion system
 * @param jacobian  Reaction Jacobian (n x n, row-major)
 * @param k         Wavenumber [1/m]
 * @param lambda    [out] Growth rates (n values, sorted by real part)
 * @return          0 on success
 */
int ds_dispersion_relation(const ReactionDiffusionSystem *rds,
                           const double *jacobian, double k,
                           double lambda[]);

/* ---------------------------------------------------------------------------
 * L3: Mathematical Structures -- Reaction-diffusion integrators
 * ------------------------------------------------------------------------- */

/**
 * @brief Compute the reaction terms R_k(c) for the Brusselator model.
 *
 * The Brusselator is a canonical model of chemical oscillations
 * and pattern formation:
 *
 *   A -> X           (constant supply)
 *   B + X -> Y + D   (bimolecular reaction)
 *   2X + Y -> 3X     (autocatalytic trimolecular step)
 *   X -> E           (decay)
 *
 * With scaled concentrations (u, v):
 *   du/dt = a - (b+1)u + u^2 v + D_u grad^2 u
 *   dv/dt = bu - u^2 v + D_v grad^2 v
 *
 * @param concentrations Current concentrations (length num_species)
 * @param reactions      [out] Reaction rates R_k (length num_species)
 * @param params         Model parameters: params[0]=a, params[1]=b
 * @param num_species    Number of species (should be 2)
 *
 * Reference: Prigogine & Lefever (1968)
 */
void ds_brusselator_reactions(const double concentrations[],
                              double reactions[],
                              const double params[], int num_species);

/**
 * @brief Compute the reaction terms for the Oregonator model.
 *
 * The Oregonator is a simplified model of the Belousov-Zhabotinsky
 * reaction, which produces chemical oscillations and waves:
 *
 *   A + Y -> X + P         rate = k1 a y
 *   X + Y -> 2P            rate = k2 x y
 *   A + X -> 2X + 2Z       rate = k3 a x
 *   2X -> A + P            rate = k4 x^2
 *   B + Z -> (f/2) Y       rate = kc b z
 *
 * @param concentrations [X, Y, Z, A, B, P] (6 species)
 * @param reactions      [out] Reaction rates
 * @param params         Rate constants [k1, k2, k3, k4, kc, f]
 * @param num_species    Should be 6
 *
 * Reference: Field & Noyes (1974) J. Chem. Phys. 60, 1877
 */
void ds_oregonator_reactions(const double concentrations[],
                             double reactions[],
                             const double params[], int num_species);

/**
 * @brief Perform one time step of reaction-diffusion integration.
 *
 * Uses operator splitting: reaction step (RK4) followed by
 * diffusion step (Crank-Nicolson implicit).
 *
 * @param rds    Reaction-diffusion system
 * @param dt     Time step [s]
 * @return       0 on success
 */
int ds_rxn_diffusion_step(ReactionDiffusionSystem *rds, double dt);

/**
 * @brief Simulate pattern formation from a small perturbation.
 *
 * Starts from the homogeneous steady state plus small random noise
 * and integrates until a stable pattern emerges or max_time is reached.
 *
 * @param rds        Reaction-diffusion system (at homogeneous SS)
 * @param noise_amp  Amplitude of initial perturbation
 * @param dt         Time step [s]
 * @param max_time   Maximum simulation time [s]
 * @param pattern    [out] Final spatial pattern
 * @return           0 on success
 */
int ds_simulate_pattern(ReactionDiffusionSystem *rds, double noise_amp,
                        double dt, double max_time,
                        DissipativePattern *pattern);

/* ---------------------------------------------------------------------------
 * L4: Fundamental Laws -- Order through fluctuations
 * ------------------------------------------------------------------------- */

/**
 * @brief Compute the entropy production of a patterned state.
 *
 * For a reaction-diffusion system, the local entropy production is:
 *   sigma(x) = sigma_diff + sigma_react
 *   sigma_diff(x) = R sum_k (D_k/c_k) |grad c_k|^2
 *   sigma_react(x) = sum_rho r_rho * (A_rho/T)
 *
 * The total entropy production is the spatial integral.
 *
 * @param rds    Reaction-diffusion system with concentrations set
 * @param sigma  [out] Total entropy production [W/K]
 * @return       0 on success
 */
int ds_pattern_entropy_production(const ReactionDiffusionSystem *rds,
                                  double *sigma);

/**
 * @brief Verify that patterned states have higher entropy production.
 *
 * A key prediction of non-equilibrium thermodynamics: structured
 * dissipative states produce more entropy than the corresponding
 * unstructured (homogeneous) state at the same boundary conditions.
 * This is the thermodynamic signature of self-organization.
 *
 * @param pattern        Patterned state
 * @param homogeneous    Homogeneous state at same constraints
 * @return               sigma_pattern - sigma_homogeneous (positive = verified)
 */
double ds_entropy_production_excess(const DissipativePattern *pattern,
                                    const ReactionDiffusionSystem *homogeneous);

/* ---------------------------------------------------------------------------
 * L5: Algorithms -- Bifurcation tracking and linear stability
 * ------------------------------------------------------------------------- */

/**
 * @brief Compute the Jacobian matrix of reaction terms.
 *
 * J_{kl} = partial(R_k) / partial(c_l) evaluated at the current
 * concentrations.
 *
 * Uses central finite differences for numerical differentiation.
 *
 * @param rds       Reaction-diffusion system
 * @param jacobian  [out] Jacobian matrix (n x n, row-major)
 * @return          0 on success
 */
int ds_reaction_jacobian(const ReactionDiffusionSystem *rds, double *jacobian);

/**
 * @brief Perform linear stability analysis of a homogeneous steady state.
 *
 * Computes eigenvalues of the Jacobian and determines stability.
 * For a stable focus, the system exhibits damped oscillations;
 * for a Hopf bifurcation, sustained oscillations emerge.
 *
 * @param rds      Reaction-diffusion system at steady state
 * @param analysis [out] Bifurcation analysis result
 * @return         0 on success
 */
int ds_linear_stability_analysis(const ReactionDiffusionSystem *rds,
                                 BifurcationAnalysis *analysis);

/**
 * @brief Track a bifurcation by continuation in the control parameter.
 *
 * Uses pseudo-arclength continuation to follow a branch of steady
 * states as the control parameter varies, detecting bifurcation points
 * where the stability changes.
 *
 * @param rds             Reaction-diffusion system
 * @param param_idx       Index of the control parameter
 * @param param_start     Initial parameter value
 * @param param_end       Final parameter value
 * @param n_steps         Number of continuation steps
 * @param bifurcations    [out] Array of detected bifurcations
 * @param max_bifurc      Maximum number of bifurcations to store
 * @return                Number of bifurcations found
 */
int ds_track_bifurcation(ReactionDiffusionSystem *rds,
                         int param_idx, double param_start,
                         double param_end, int n_steps,
                         BifurcationAnalysis bifurcations[],
                         int max_bifurc);

/* ---------------------------------------------------------------------------
 * L6: Canonical Models
 * ------------------------------------------------------------------------- */

/**
 * @brief Set up the Brusselator model with standard parameters.
 *
 * Initializes a 2-species reaction-diffusion system for the
 * Brusselator with parameters in the Turing-unstable regime:
 *   a = 2.0, b = 5.5 (Turing instability condition: b > 1 + a^2)
 *   D_u = 1.0, D_v = 16.0 (D_v/D_u > critical ratio)
 *
 * @param rds            [out] Initialized system
 * @param grid_points    Number of spatial grid points
 * @param domain_length  Physical length [m]
 */
void ds_setup_brusselator(ReactionDiffusionSystem *rds,
                          int grid_points, double domain_length);

/**
 * @brief Set up the Oregonator model for the BZ reaction.
 *
 * Initializes the Oregonator with parameters that produce
 * relaxation oscillations.
 *
 * @param rds            [out] Initialized system
 */
void ds_setup_oregonator(ReactionDiffusionSystem *rds);

/**
 * @brief Set up the Schnakenberg model (simplified activator-inhibitor).
 *
 * The Schnakenberg model is a minimal 2-species model for
 * pattern formation in biological systems:
 *   du/dt = gamma(a - u + u^2 v) + D_u grad^2 u
 *   dv/dt = gamma(b - u^2 v) + D_v grad^2 v
 *
 * @param rds            [out] Initialized system
 * @param grid_points    Number of spatial grid points
 * @param domain_length  Physical length [m]
 */
void ds_setup_schnakenberg(ReactionDiffusionSystem *rds,
                           int grid_points, double domain_length);

/* Utility */

void ds_init_rds(ReactionDiffusionSystem *rds, int num_species,
                 int grid_points, double domain_length, double temperature);
void ds_free_rds(ReactionDiffusionSystem *rds);
void ds_init_pattern(DissipativePattern *pattern, int num_species, int grid_points);
void ds_free_pattern(DissipativePattern *pattern);

#ifdef __cplusplus
}
#endif

#endif /* DISSIPATIVE_SYSTEMS_H */