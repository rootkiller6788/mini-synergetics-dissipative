#ifndef PDS_CORE_H
#define PDS_CORE_H

#include <stdbool.h>
#include <stddef.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * Prigogine Dissipative Structures - Core Definitions
 *
 * Based on the foundational works of:
 *   Ilya Prigogine (Nobel Prize in Chemistry, 1977)
 *     - "Thermodynamic Theory of Structure, Stability and Fluctuations" (1971)
 *     - "Self-Organization in Non-Equilibrium Systems" (Nicolis and Prigogine, 1977)
 *     - "Order Out of Chaos" (Prigogine and Stengers, 1984)
 *   Paul Glansdorff and Ilya Prigogine
 *     - "Thermodynamic Theory of Structure, Stability and Fluctuations" (1971)
 *   Gregoire Nicolis
 *     - "Introduction to Nonlinear Science" (1995)
 *   Lars Onsager (Nobel Prize 1968)
 *     - Reciprocal relations in irreversible processes (1931)
 * ============================================================================ */

/* --- L1: Dissipative Structure Classification --- */

typedef enum {
    PDS_NEAR_EQUILIBRIUM = 0,
    PDS_WEAKLY_NONEQUIL = 1,
    PDS_FAR_EQUILIBRIUM = 2,
    PDS_CRITICAL_POINT = 3,
    PDS_DISSIPATIVE_STRUCTURE = 4
} PDSThermoRegime;

typedef enum {
    PDS_SPATIAL_PATTERN = 0,
    PDS_TEMPORAL_RHYTHM = 1,
    PDS_SPATIOTEMPORAL = 2,
    PDS_CHAOTIC_DISSIPATIVE = 3,
    PDS_STATIONARY_DISSIPATIVE = 4
} PDSStructureType;

typedef enum {
    PDS_STABLE_THERMO_BRANCH = 0,
    PDS_MARGINALLY_STABLE = 1,
    PDS_UNSTABLE_BRANCH = 2,
    PDS_NEW_BRANCH_EMERGED = 3,
    PDS_MULTIPLE_BRANCHES = 4
} PDSBranchState;

typedef enum {
    PDS_BENARD_CONVECTION = 0,
    PDS_BZ_REACTION = 1,
    PDS_BRUSSELATOR = 2,
    PDS_TURING_PATTERN = 3,
    PDS_LASER_THRESHOLD = 4,
    PDS_BIOLOGICAL_RHYTHM = 5,
    PDS_ECOSYSTEM_STRUCTURE = 6,
    PDS_SOCIAL_DISSIPATIVE = 7
} PDSCanonicalSystem;

/* --- L2: Entropy Balance (Prigogine central equation dS = d_i S + d_e S) --- */

typedef struct {
    double S;
    double d_i_S;
    double d_e_S;
    double dS_dt;
    double entropy_production_rate;
    double entropy_flow_rate;
    double excess_entropy;
    double excess_entropy_production;
} PDSEntropyBalance;

/* --- L3: Mathematical Structure - Vector --- */

typedef struct {
    double* data;
    int size;
    int capacity;
} PDSVector;

/* --- L3: Mathematical Structure - Matrix --- */

typedef struct {
    double** m;
    int rows;
    int cols;
} PDSMatrix;

/* --- L1: Dissipative Structure Core Object --- */

typedef struct {
    int id;
    char* name;
    PDSCanonicalSystem canonical_type;
    PDSThermoRegime regime;
    PDSStructureType structure_type;
    PDSBranchState branch_state;

    PDSVector* x;
    PDSVector* x_equilibrium;

    PDSEntropyBalance entropy;

    double* parameters;
    char** parameter_names;
    int n_parameters;

    PDSMatrix* L_matrix;

    int n_species;
    int n_reactions;
    double* rate_constants;
    double** stoichiometry;

    int spatial_dims;
    double* diffusion_coeffs;
    double* spatial_grid;
    int n_grid_points;

    double* eigenvalues;
    int n_eigenvalues;
    double max_real_eigenvalue;
    double lyapunov_exponent;

    PDSVector** history;
    int history_capacity;
    int history_count;
    double* time_history;

    double current_time;
    double creation_entropy;
    bool is_dissipative_structure;
} PDSDissipativeSystem;

/* --- Core API --- */

PDSDissipativeSystem* pds_system_create(const char* name,
                                         PDSCanonicalSystem canonical_type,
                                         int n_species, int n_reactions,
                                         int n_params);
void pds_system_free(PDSDissipativeSystem* sys);
void pds_system_set_parameter(PDSDissipativeSystem* sys,
                               int param_idx, double value);
double pds_system_get_parameter(const PDSDissipativeSystem* sys, int param_idx);
void pds_system_set_state(PDSDissipativeSystem* sys, const double* x);
void pds_system_get_state(const PDSDissipativeSystem* sys, double* x_out);
void pds_system_set_equilibrium(PDSDissipativeSystem* sys, const double* x_eq);
void pds_system_record_state(PDSDissipativeSystem* sys);
void pds_system_clear_history(PDSDissipativeSystem* sys);
void pds_system_print(const PDSDissipativeSystem* sys);

/* --- Vector Operations --- */

PDSVector* pds_vector_create(int size);
void pds_vector_free(PDSVector* v);
void pds_vector_set(PDSVector* v, int i, double value);
double pds_vector_get(const PDSVector* v, int i);
void pds_vector_copy(const PDSVector* src, PDSVector* dst);
double pds_vector_norm(const PDSVector* v);
double pds_vector_dot(const PDSVector* a, const PDSVector* b);
void pds_vector_scale(PDSVector* v, double alpha);
void pds_vector_add(PDSVector* a, const PDSVector* b);
void pds_vector_sub(PDSVector* a, const PDSVector* b);
void pds_vector_print(const PDSVector* v, const char* label);

/* --- Matrix Operations --- */

PDSMatrix* pds_matrix_create(int rows, int cols);
void pds_matrix_free(PDSMatrix* mat);
void pds_matrix_set(PDSMatrix* mat, int i, int j, double value);
double pds_matrix_get(const PDSMatrix* mat, int i, int j);
void pds_matrix_vec_mul(const PDSMatrix* A, const PDSVector* x, PDSVector* result);
void pds_matrix_print(const PDSMatrix* mat, const char* label);
void pds_matrix_copy(const PDSMatrix* src, PDSMatrix* dst);
void pds_matrix_add(PDSMatrix* A, const PDSMatrix* B, PDSMatrix* result);
void pds_matrix_scale(PDSMatrix* A, double alpha, PDSMatrix* result);
double pds_matrix_determinant_2x2(const PDSMatrix* A);
double pds_matrix_trace(const PDSMatrix* A);
void pds_matrix_transpose(const PDSMatrix* A, PDSMatrix* result);

/* --- Utility --- */

double pds_sign(double x);
double pds_clamp(double x, double lo, double hi);
bool pds_is_near(double a, double b, double tol);

#endif /* PDS_CORE_H */