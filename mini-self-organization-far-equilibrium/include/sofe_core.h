#ifndef SOFE_CORE_H
#define SOFE_CORE_H

#include <stdbool.h>
#include <stddef.h>

/* ============================================================================
 * Self-Organization Far-from-Equilibrium (SOFE) — Core Types & Definitions
 *
 * Foundational works:
 *   Ilya Prigogine       — Thermodynamics of Irreversible Processes (1955)
 *                         — Self-Organization in Non-Equilibrium Systems (1977)
 *   Hermann Haken        — Synergetics: An Introduction (1977)
 *                         — Advanced Synergetics (1983)
 *   Grégoire Nicolis     — Self-Organization in Nonequilibrium Systems (1977)
 *   Manfred Eigen        — The Hypercycle (1979)
 *   Peter Glansdorff     — Thermodynamic Theory of Structure, Stability
 *                          and Fluctuations (1971)
 *   Lars Onsager         — Reciprocal Relations in Irreversible Processes (1931)
 *   Alan Turing          — The Chemical Basis of Morphogenesis (1952)
 *   Per Bak              — How Nature Works: Self-Organized Criticality (1996)
 *   Stuart Kauffman      — At Home in the Universe (1995)
 * ============================================================================ */

/* --- L1: Core Enumerations (Definitions) --- */

/** Thermodynamic regime classification.
 *  The distance from equilibrium determines the qualitative behavior
 *  of the system. Near equilibrium, Onsager reciprocity and minimum
 *  entropy production hold. Far from equilibrium, new structures emerge. */
typedef enum {
    SOFE_REGIME_EQUILIBRIUM          = 0,  /* Isolated system, max entropy */
    SOFE_REGIME_NEAR_EQUILIBRIUM     = 1,  /* Linear irreversible thermodynamics */
    SOFE_REGIME_NONLINEAR_NONEQ      = 2,  /* Nonlinear region, bifurcations possible */
    SOFE_REGIME_FAR_EQUILIBRIUM      = 3,  /* Dissipative structures emerge */
    SOFE_REGIME_CRITICAL             = 4,  /* At bifurcation/critical point */
    SOFE_REGIME_CHAOTIC_NONEQ        = 5   /* Deterministic chaos far from equilibrium */
} SOFERegime;

/** Type of dissipative structure.
 *  Prigogine classified dissipative structures by their spatial and
 *  temporal characteristics. */
typedef enum {
    SOFE_DS_SPATIAL     = 0,  /* Stationary spatial pattern (Turing) */
    SOFE_DS_TEMPORAL    = 1,  /* Temporal oscillation (limit cycle) */
    SOFE_DS_SPATIOTEMP  = 2,  /* Traveling wave, spiral wave */
    SOFE_DS_CHAOTIC     = 3,  /* Deterministic chaos (strange attractor) */
    SOFE_DS_HIERARCHICAL = 4  /* Multi-level self-organized structure */
} SOFEDissipativeType;

/** Order parameter type.
 *  In Haken's synergetics, order parameters enslave the fast-relaxing
 *  variables near instability points. */
typedef enum {
    SOFE_ORDER_SCALAR         = 0,  /* Single scalar order parameter */
    SOFE_ORDER_VECTOR         = 1,  /* Multi-component order parameter */
    SOFE_ORDER_COMPLEX        = 2,  /* Complex amplitude (Hopf) */
    SOFE_ORDER_TENSOR         = 3,  /* Tensor order (liquid crystals) */
    SOFE_ORDER_CONSERVED      = 4,  /* Conserved order parameter (Model B) */
    SOFE_ORDER_NONCONSERVED   = 5   /* Non-conserved (Model A) */
} SOFEOrderType;

/** Bifurcation type in far-from-equilibrium systems. */
typedef enum {
    SOFE_BIF_NONE              = 0,
    SOFE_BIF_SADDLE_NODE       = 1,  /* Turning point, fold */
    SOFE_BIF_PITCHFORK_SUPER   = 2,  /* Supercritical pitchfork (Z2 breaking) */
    SOFE_BIF_PITCHFORK_SUB     = 3,  /* Subcritical pitchfork */
    SOFE_BIF_HOPF_SUPER        = 4,  /* Supercritical Hopf -> limit cycle */
    SOFE_BIF_HOPF_SUB          = 5,  /* Subcritical Hopf */
    SOFE_BIF_TRANSCRITICAL     = 6,  /* Exchange of stability */
    SOFE_BIF_TURING            = 7,  /* Turing instability (spatial) */
    SOFE_BIF_TAKENS_BOGDANOV   = 8,  /* Double-zero eigenvalue */
    SOFE_BIF_CODIM2            = 9   /* Higher codimension */
} SOFEBifurcationType;

/** Stability classification of steady states. */
typedef enum {
    SOFE_STABLE_NODE        = 0,
    SOFE_STABLE_FOCUS       = 1,
    SOFE_UNSTABLE_NODE      = 2,
    SOFE_UNSTABLE_FOCUS     = 3,
    SOFE_SADDLE             = 4,
    SOFE_CENTER             = 5,
    SOFE_STABLE_LIMIT_CYCLE = 6,
    SOFE_UNSTABLE_LIMIT_CYCLE = 7,
    SOFE_STRANGE_ATTRACTOR  = 8
} SOFEStabilityClass;

/** Symmetry breaking type (Landau theory classification). */
typedef enum {
    SOFE_SYMMETRY_NONE      = 0,
    SOFE_SYMMETRY_Z2        = 1,  /* Ising-type (up/down) */
    SOFE_SYMMETRY_O2        = 2,  /* Planar rotation */
    SOFE_SYMMETRY_SO3       = 3,  /* 3D rotation */
    SOFE_SYMMETRY_TRANS     = 4,  /* Translational (crystal) */
    SOFE_SYMMETRY_GAUGE     = 5,  /* Gauge symmetry */
    SOFE_SYMMETRY_CHIRAL    = 6,  /* Left/right handedness */
    SOFE_SYMMETRY_SCALE     = 7   /* Scale invariance breaking */
} SOFESymmetryType;

/** Fluctuation type classification. */
typedef enum {
    SOFE_FLUCT_THERMAL      = 0,  /* Thermal fluctuations (Johnson-Nyquist) */
    SOFE_FLUCT_SHOT          = 1,  /* Shot noise (discrete events) */
    SOFE_FLUCT_CRITICAL      = 2,  /* Critical fluctuations (diverging) */
    SOFE_FLUCT_NUCLEATION    = 3,  /* Nucleation events */
    SOFE_FLUCT_EXTERNAL      = 4   /* Environment-driven fluctuations */
} SOFEFluctuationType;

/** Feedback type in self-organizing systems. */
typedef enum {
    SOFE_FEEDBACK_POSITIVE   = 0,  /* Amplifying, autocatalytic */
    SOFE_FEEDBACK_NEGATIVE   = 1,  /* Stabilizing, homeostatic */
    SOFE_FEEDBACK_CROSS      = 2,  /* Cross-catalytic */
    SOFE_FEEDBACK_DELAYED    = 3,  /* Time-delayed feedback */
    SOFE_FEEDBACK_GLOBAL     = 4   /* Global (mean-field) coupling */
} SOFEFeedbackType;

/* --- L2: Core Numeric Types --- */

/** Complex number for Hopf bifurcation analysis. */
typedef struct {
    double real;
    double imag;
} SOFEComplex;

/** Dense matrix. */
typedef struct {
    double* data;
    int rows;
    int cols;
} SOFEMatrix;

/** Vector. */
typedef struct {
    double* components;
    int dim;
} SOFEVector;

/** Sparse matrix entry for reaction networks. */
typedef struct {
    int row;
    int col;
    double value;
} SOFESparseEntry;

/** Sparse matrix (COO format). */
typedef struct {
    SOFESparseEntry* entries;
    int n_entries;
    int capacity;
    int rows;
    int cols;
} SOFESparseMatrix;

/* --- L3: Reaction Network Type --- */

/** One chemical species in a far-from-equilibrium reaction network. */
typedef struct {
    char* name;
    double concentration;
    double diffusion_coefficient;
    double production_rate;
    double degradation_rate;
    double equilibrium_value;
    bool is_autocatalytic;
    bool is_inhibitor;
    bool is_activator;
    bool is_order_parameter;
} SOFESpecies;

/** One reaction in a network (stoichiometry + kinetics). */
typedef struct {
    char* name;
    int* reactant_indices;
    double* reactant_stoich;
    int n_reactants;
    int* product_indices;
    double* product_stoich;
    int n_products;
    double rate_constant;
    double activation_energy;
    double reaction_enthalpy;
    bool is_autocatalytic_step;
    bool is_rate_limiting;
    int kinetic_order;
    int kinetic_type;            /* 0=mass action, 1=Michaelis, 2=Hill */
    double hill_coefficient;     /* For Hill-type kinetics */
    double michaelis_constant;   /* For Michaelis-Menten */
} SOFEReaction;

/** Full reaction network (chemical reaction system far from equilibrium). */
typedef struct {
    char* name;
    SOFESpecies* species;
    int n_species;
    int species_capacity;
    SOFEReaction* reactions;
    int n_reactions;
    int reaction_capacity;
    double temperature;
    double volume;
    double total_mass;
    double* stoichiometric_matrix;   /* n_species x n_reactions */
    double* jacobian;                /* n_species x n_species */
    double* flux_vector;             /* length n_reactions */
    double* affinity_vector;         /* length n_reactions */
    bool jacobian_computed;
} SOFEReactionNetwork;

/* --- L3: Spatial Domain (for reaction-diffusion) --- */

typedef enum {
    SOFE_GRID_1D = 0,
    SOFE_GRID_2D = 1,
    SOFE_GRID_3D = 2,
    SOFE_GRID_HEXAGONAL = 3,
    SOFE_GRID_UNSTRUCTURED = 4
} SOFEGridType;

typedef struct {
    SOFEGridType type;
    int nx, ny, nz;
    double dx, dy, dz;
    int n_points;
    double* coordinates;
    int* neighbor_indices;
    int* neighbor_offsets;
    int max_neighbors;
    bool is_periodic_x;
    bool is_periodic_y;
    bool is_periodic_z;
} SOFEGrid;

/** Spatial concentration field on a grid. */
typedef struct {
    SOFEGrid* grid;
    double** field;          /* field[species_idx][grid_point] */
    int n_species;
    double time;
} SOFEField;

/* --- L3: Order Parameter Structure --- */

typedef struct {
    char* name;
    SOFEOrderType type;
    SOFESymmetryType broken_symmetry;
    int n_components;
    double* amplitude;
    double* phase;
    double critical_exponent_beta;
    double susceptibility;
    double correlation_length;
    bool is_slow_variable;
    double relaxation_time;
} SOFEOrderParameter;

/* --- Core System (the self-organizing system itself) --- */

typedef struct {
    char* name;
    SOFERegime regime;

    /* Reaction network (deterministic dynamics) */
    SOFEReactionNetwork* network;

    /* Spatial extension */
    SOFEField* concentration_field;

    /* Order parameters (Haken's synergetics) */
    SOFEOrderParameter* order_params;
    int n_order_params;
    int order_param_capacity;

    /* Thermodynamic state */
    double entropy;
    double entropy_production_rate;
    double entropy_flow_rate;
    double entropy_internal_rate;
    double free_energy;
    double free_energy_density;
    double temperature;
    double chemical_potential_sum;

    /* Stability and bifurcation */
    SOFEStabilityClass stability;
    SOFEBifurcationType last_bifurcation;
    double* eigenvalues;
    int n_eigenvalues;
    double lyapunov_exponent_max;

    /* Control parameters */
    double* control_params;
    char** control_param_names;
    int n_control_params;
    double bifurcation_parameter;

    /* History */
    double* history_entropy_prod;
    double* history_order_param;
    int history_len;
    int history_capacity;

    /* Statistics */
    double fluctuation_variance;
    double fluctuation_correlation_time;
    bool has_macroscopic_order;
} SOFESystem;

/* --- Core API --- */

SOFESystem* sofe_system_create(const char* name);
void sofe_system_free(SOFESystem* sys);
void sofe_system_set_regime(SOFESystem* sys, SOFERegime regime);
void sofe_system_set_temperature(SOFESystem* sys, double T);

SOFEReactionNetwork* sofe_network_create(const char* name);
void sofe_network_free(SOFEReactionNetwork* net);
int sofe_network_add_species(SOFEReactionNetwork* net, const char* name,
                              double init_conc, double diffusion_coeff);
int sofe_network_add_reaction(SOFEReactionNetwork* net, const char* name,
                               double rate_const, double activation_energy);
void sofe_network_set_reactants(SOFEReactionNetwork* net, int rxn_idx,
                                 int* species_indices, double* stoich, int n);
void sofe_network_set_products(SOFEReactionNetwork* net, int rxn_idx,
                                int* species_indices, double* stoich, int n);
void sofe_network_set_kinetics(SOFEReactionNetwork* net, int rxn_idx,
                                int kinetic_type, int order,
                                double hill_n, double km);
void sofe_network_set_autocatalytic(SOFEReactionNetwork* net, int rxn_idx,
                                     bool is_auto);
void sofe_network_compute_stoichiometric_matrix(SOFEReactionNetwork* net);
void sofe_network_compute_jacobian(SOFEReactionNetwork* net);
void sofe_network_compute_fluxes(SOFEReactionNetwork* net);
double sofe_network_reaction_rate(SOFEReactionNetwork* net, int rxn_idx);

void sofe_species_set_autocatalytic(SOFESpecies* sp, bool is_auto);
void sofe_species_set_diffusion(SOFESpecies* sp, double D);
void sofe_species_set_equilibrium(SOFESpecies* sp, double eq_val);

SOFEGrid* sofe_grid_create_1d(int nx, double dx, bool periodic);
SOFEGrid* sofe_grid_create_2d(int nx, int ny, double dx, double dy,
                               bool periodic_x, bool periodic_y);
SOFEGrid* sofe_grid_create_3d(int nx, int ny, int nz,
                               double dx, double dy, double dz,
                               bool periodic_x, bool periodic_y, bool periodic_z);
void sofe_grid_free(SOFEGrid* grid);
int sofe_grid_index(SOFEGrid* grid, int ix, int iy, int iz);
void sofe_grid_get_neighbors(SOFEGrid* grid, int point_idx,
                              int* neighbors, int* n_found);

SOFEField* sofe_field_create(SOFEGrid* grid, int n_species);
void sofe_field_free(SOFEField* field);
void sofe_field_set(SOFEField* field, int species, int point, double val);
double sofe_field_get(SOFEField* field, int species, int point);
double sofe_field_laplacian(SOFEField* field, int species, int point);
void sofe_field_copy(SOFEField* dst, SOFEField* src);

int sofe_system_add_order_parameter(SOFESystem* sys, const char* name,
                                     SOFEOrderType type,
                                     SOFESymmetryType broken_symmetry,
                                     int n_components);
void sofe_system_set_order_amplitude(SOFESystem* sys, int op_idx,
                                      int component, double amplitude);
double sofe_system_get_order_amplitude(SOFESystem* sys, int op_idx,
                                        int component);
void sofe_system_compute_order_parameters(SOFESystem* sys);

int sofe_system_add_control_param(SOFESystem* sys, const char* name,
                                   double initial_value);
void sofe_system_set_control_param(SOFESystem* sys, int idx, double value);
double sofe_system_get_control_param(SOFESystem* sys, int idx);

void sofe_system_record_state(SOFESystem* sys);
void sofe_system_clear_history(SOFESystem* sys);

SOFEStabilityClass sofe_classify_stability(double* eigenvalues, int n);
const char* sofe_stability_name(SOFEStabilityClass sc);
const char* sofe_bifurcation_name(SOFEBifurcationType bt);
const char* sofe_regime_name(SOFERegime reg);
const char* sofe_symmetry_name(SOFESymmetryType sym);

SOFEComplex sofe_complex_add(SOFEComplex a, SOFEComplex b);
SOFEComplex sofe_complex_mul(SOFEComplex a, SOFEComplex b);
double sofe_complex_abs(SOFEComplex z);

void sofe_matrix_print(SOFEMatrix* m);
void sofe_vector_print(SOFEVector* v);
double sofe_vector_norm(double* v, int n);
double sofe_vector_dot(double* a, double* b, int n);
void sofe_matrix_vector_mul(double* out, SOFEMatrix* A, double* x);
void sofe_system_print_summary(SOFESystem* sys);

#endif
