#include "pds_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

/* ============================================================================
 * PDS Core ? Dissipative System Creation, Management, and Linear Algebra
 *
 * Implementation of Prigogine's dissipative structure framework:
 *   - System lifecycle management
 *   - Vector and matrix operations for state-space representation
 *   - History tracking for time-series analysis
 *
 * Reference: Nicolis & Prigogine, "Self-Organization in Non-Equilibrium
 *            Systems", Wiley (1977), Chapters 1-4.
 * ============================================================================ */

/* ============================================================================
 * Dissipative System Lifecycle
 * ============================================================================ */

PDSDissipativeSystem* pds_system_create(const char* name,
                                         PDSCanonicalSystem canonical_type,
                                         int n_species, int n_reactions,
                                         int n_params) {
    if (!name || n_species < 0 || n_reactions < 0 || n_params < 0) {
        return NULL;
    }

    PDSDissipativeSystem* sys = (PDSDissipativeSystem*)calloc(1, sizeof(PDSDissipativeSystem));
    if (!sys) return NULL;

    sys->id = (int)((size_t)sys % 100000);
    sys->name = (char*)malloc(strlen(name) + 1);
    if (!sys->name) { free(sys); return NULL; }
    strcpy(sys->name, name);

    sys->canonical_type = canonical_type;
    sys->regime = PDS_NEAR_EQUILIBRIUM;
    sys->structure_type = PDS_STATIONARY_DISSIPATIVE;
    sys->branch_state = PDS_STABLE_THERMO_BRANCH;

    sys->n_species = n_species;
    sys->n_reactions = n_reactions;
    sys->n_parameters = n_params;

    /* State vectors */
    sys->x = pds_vector_create(n_species);
    sys->x_equilibrium = pds_vector_create(n_species);

    /* Parameters */
    if (n_params > 0) {
        sys->parameters = (double*)calloc(n_params, sizeof(double));
        sys->parameter_names = (char**)calloc(n_params, sizeof(char*));
    }

    /* Reaction network */
    if (n_reactions > 0) {
        sys->rate_constants = (double*)calloc(n_reactions, sizeof(double));
        sys->stoichiometry = (double**)calloc(n_species, sizeof(double*));
        for (int i = 0; i < n_species; i++) {
            sys->stoichiometry[i] = (double*)calloc(n_reactions, sizeof(double));
        }
    }

    /* History */
    sys->history_capacity = 1024;
    sys->history = (PDSVector**)calloc(sys->history_capacity, sizeof(PDSVector*));
    sys->time_history = (double*)calloc(sys->history_capacity, sizeof(double));
    sys->history_count = 0;

    sys->current_time = 0.0;
    sys->creation_entropy = 0.0;
    sys->is_dissipative_structure = false;
    sys->max_real_eigenvalue = 0.0;
    sys->lyapunov_exponent = 0.0;
    sys->n_eigenvalues = 0;
    sys->eigenvalues = NULL;

    /* Spatial defaults */
    sys->spatial_dims = 0;
    sys->n_grid_points = 0;
    sys->diffusion_coeffs = NULL;
    sys->spatial_grid = NULL;

    sys->L_matrix = NULL;

    return sys;
}

void pds_system_free(PDSDissipativeSystem* sys) {
    if (!sys) return;
    free(sys->name);
    pds_vector_free(sys->x);
    pds_vector_free(sys->x_equilibrium);

    /* Free history */
    for (int i = 0; i < sys->history_count; i++) {
        pds_vector_free(sys->history[i]);
    }
    free(sys->history);
    free(sys->time_history);

    /* Free parameters */
    if (sys->parameter_names) {
        for (int i = 0; i < sys->n_parameters; i++) {
            free(sys->parameter_names[i]);
        }
        free(sys->parameter_names);
    }
    free(sys->parameters);

    /* Free reaction network */
    if (sys->stoichiometry) {
        for (int i = 0; i < sys->n_species; i++) {
            free(sys->stoichiometry[i]);
        }
        free(sys->stoichiometry);
    }
    free(sys->rate_constants);

    /* Free eigenvalues */
    free(sys->eigenvalues);

    /* Free spatial data */
    free(sys->diffusion_coeffs);
    free(sys->spatial_grid);

    /* Free Onsager matrix */
    if (sys->L_matrix) {
        pds_matrix_free(sys->L_matrix);
    }

    free(sys);
}

/* ============================================================================
 * System State Management
 * ============================================================================ */

void pds_system_set_parameter(PDSDissipativeSystem* sys,
                               int param_idx, double value) {
    if (!sys || param_idx < 0 || param_idx >= sys->n_parameters) return;
    sys->parameters[param_idx] = value;
}

double pds_system_get_parameter(const PDSDissipativeSystem* sys, int param_idx) {
    if (!sys || param_idx < 0 || param_idx >= sys->n_parameters) return 0.0;
    return sys->parameters[param_idx];
}

void pds_system_set_state(PDSDissipativeSystem* sys, const double* x) {
    if (!sys || !x) return;
    for (int i = 0; i < sys->n_species; i++) {
        pds_vector_set(sys->x, i, x[i]);
    }
}

void pds_system_get_state(const PDSDissipativeSystem* sys, double* x_out) {
    if (!sys || !x_out) return;
    for (int i = 0; i < sys->n_species; i++) {
        x_out[i] = pds_vector_get(sys->x, i);
    }
}

void pds_system_set_equilibrium(PDSDissipativeSystem* sys, const double* x_eq) {
    if (!sys || !x_eq) return;
    for (int i = 0; i < sys->n_species; i++) {
        pds_vector_set(sys->x_equilibrium, i, x_eq[i]);
    }
}

/* ============================================================================
 * History Recording ? Time Series for Stability/Entropy Analysis
 * ============================================================================ */

void pds_system_record_state(PDSDissipativeSystem* sys) {
    if (!sys || sys->history_count >= sys->history_capacity) return;

    PDSVector* snapshot = pds_vector_create(sys->n_species);
    pds_vector_copy(sys->x, snapshot);
    sys->history[sys->history_count] = snapshot;
    sys->time_history[sys->history_count] = sys->current_time;
    sys->history_count++;
}

void pds_system_clear_history(PDSDissipativeSystem* sys) {
    if (!sys) return;
    for (int i = 0; i < sys->history_count; i++) {
        pds_vector_free(sys->history[i]);
        sys->history[i] = NULL;
    }
    sys->history_count = 0;
}

/* ============================================================================
 * System Display
 * ============================================================================ */

void pds_system_print(const PDSDissipativeSystem* sys) {
    if (!sys) { printf("PDSDissipativeSystem: NULL\n"); return; }

    const char* regime_names[] = {
        "Near Equilibrium", "Weakly Nonequilibrium",
        "Far from Equilibrium", "Critical Point",
        "Dissipative Structure"
    };
    const char* branch_names[] = {
        "Stable Thermodynamic Branch", "Marginally Stable",
        "Unstable Branch", "New Branch Emerged",
        "Multiple Branches"
    };
    const char* canonical_names[] = {
        "Benard Convection", "BZ Reaction", "Brusselator",
        "Turing Pattern", "Laser Threshold", "Biological Rhythm",
        "Ecosystem Structure", "Social Dissipative"
    };

    printf("=== PDSDissipativeSystem [%d] '%s' ===\n", sys->id, sys->name);
    printf("  Canonical Type: %s\n", canonical_names[sys->canonical_type]);
    printf("  Regime: %s\n", regime_names[sys->regime]);
    printf("  Branch State: %s\n", branch_names[sys->branch_state]);
    printf("  Species: %d, Reactions: %d, Parameters: %d\n",
           sys->n_species, sys->n_reactions, sys->n_parameters);
    printf("  Time: %.6f\n", sys->current_time);
    printf("  Entropy: S=%.6f d_iS=%.6f d_eS=%.6f P=%.6f\n",
           sys->entropy.S, sys->entropy.d_i_S, sys->entropy.d_e_S,
           sys->entropy.entropy_production_rate);
    printf("  History: %d/%d snapshots\n", sys->history_count, sys->history_capacity);
    printf("  Max Re(lambda): %.6f\n", sys->max_real_eigenvalue);
    printf("  Max Lyapunov Exponent: %.6f\n", sys->lyapunov_exponent);
    printf("  Is Dissipative Structure: %s\n",
           sys->is_dissipative_structure ? "YES" : "no");
    printf("  Spatial Dims: %d, Grid Points: %d\n",
           sys->spatial_dims, sys->n_grid_points);

    if (sys->n_species > 0 && sys->x) {
        printf("  State: [");
        for (int i = 0; i < sys->n_species && i < 5; i++) {
            printf("%.4f%s", pds_vector_get(sys->x, i),
                   i < sys->n_species - 1 ? ", " : "");
        }
        if (sys->n_species > 5) printf("...");
        printf("]\n");
    }
}

/* ============================================================================
 * Vector Operations ? L3 Mathematical Structures
 *
 * PDSVector implements the state-space for dissipative systems.
 * The state vector x contains concentrations, temperatures, etc.
 * Vector operations enable:
 *   - Distance from equilibrium: ||x - x_eq||
 *   - Rate of change direction: dx/dt
 *   - Orthogonality checks for normal modes
 * ============================================================================ */

PDSVector* pds_vector_create(int size) {
    if (size <= 0) return NULL;

    PDSVector* v = (PDSVector*)malloc(sizeof(PDSVector));
    if (!v) return NULL;

    int cap = (size < 16) ? 16 : size * 2;
    v->data = (double*)calloc(cap, sizeof(double));
    v->size = size;
    v->capacity = cap;

    if (!v->data) { free(v); return NULL; }
    return v;
}

void pds_vector_free(PDSVector* v) {
    if (!v) return;
    free(v->data);
    free(v);
}

void pds_vector_set(PDSVector* v, int i, double value) {
    if (!v || i < 0 || i >= v->size) return;
    v->data[i] = value;
}

double pds_vector_get(const PDSVector* v, int i) {
    if (!v || i < 0 || i >= v->size) return 0.0;
    return v->data[i];
}

void pds_vector_copy(const PDSVector* src, PDSVector* dst) {
    if (!src || !dst || src->size != dst->size) return;
    memcpy(dst->data, src->data, src->size * sizeof(double));
}

double pds_vector_norm(const PDSVector* v) {
    if (!v || v->size == 0) return 0.0;
    double sum_sq = 0.0;
    for (int i = 0; i < v->size; i++) {
        sum_sq += v->data[i] * v->data[i];
    }
    return sqrt(sum_sq);
}

double pds_vector_dot(const PDSVector* a, const PDSVector* b) {
    if (!a || !b || a->size != b->size) return 0.0;
    double dot = 0.0;
    for (int i = 0; i < a->size; i++) {
        dot += a->data[i] * b->data[i];
    }
    return dot;
}

void pds_vector_scale(PDSVector* v, double alpha) {
    if (!v) return;
    for (int i = 0; i < v->size; i++) {
        v->data[i] *= alpha;
    }
}

void pds_vector_add(PDSVector* a, const PDSVector* b) {
    if (!a || !b || a->size != b->size) return;
    for (int i = 0; i < a->size; i++) {
        a->data[i] += b->data[i];
    }
}

void pds_vector_sub(PDSVector* a, const PDSVector* b) {
    if (!a || !b || a->size != b->size) return;
    for (int i = 0; i < a->size; i++) {
        a->data[i] -= b->data[i];
    }
}

void pds_vector_print(const PDSVector* v, const char* label) {
    if (!v) { printf("%s: NULL\n", label ? label : "Vector"); return; }
    printf("%s [%d]: ", label ? label : "Vector", v->size);
    for (int i = 0; i < v->size && i < 8; i++) {
        printf("%.4f ", v->data[i]);
    }
    if (v->size > 8) printf("...");
    printf("\n");
}

/* ============================================================================
 * Matrix Operations ? L3 Mathematical Structures
 *
 * PDSMatrix implements:
 *   - Jacobian matrix for linear stability analysis
 *   - Onsager coefficient matrix L_ij
 *   - Diffusion matrix for reaction-diffusion systems
 *   - State transition matrices
 *
 * Matrix-vector multiplication is essential for:
 *   - dx/dt = J * (x - x_s) (linearized dynamics)
 *   - J_i = sum_j L_ij * X_j (Onsager flux-force relations)
 * ============================================================================ */

PDSMatrix* pds_matrix_create(int rows, int cols) {
    if (rows <= 0 || cols <= 0) return NULL;

    PDSMatrix* mat = (PDSMatrix*)malloc(sizeof(PDSMatrix));
    if (!mat) return NULL;

    mat->m = (double**)malloc(rows * sizeof(double*));
    if (!mat->m) { free(mat); return NULL; }

    for (int i = 0; i < rows; i++) {
        mat->m[i] = (double*)calloc(cols, sizeof(double));
        if (!mat->m[i]) {
            for (int j = 0; j < i; j++) free(mat->m[j]);
            free(mat->m);
            free(mat);
            return NULL;
        }
    }
    mat->rows = rows;
    mat->cols = cols;
    return mat;
}

void pds_matrix_free(PDSMatrix* mat) {
    if (!mat) return;
    for (int i = 0; i < mat->rows; i++) {
        free(mat->m[i]);
    }
    free(mat->m);
    free(mat);
}

void pds_matrix_set(PDSMatrix* mat, int i, int j, double value) {
    if (!mat || i < 0 || i >= mat->rows || j < 0 || j >= mat->cols) return;
    mat->m[i][j] = value;
}

double pds_matrix_get(const PDSMatrix* mat, int i, int j) {
    if (!mat || i < 0 || i >= mat->rows || j < 0 || j >= mat->cols) return 0.0;
    return mat->m[i][j];
}

void pds_matrix_vec_mul(const PDSMatrix* A, const PDSVector* x, PDSVector* result) {
    if (!A || !x || !result) return;
    if (A->cols != x->size || A->rows != result->size) return;

    for (int i = 0; i < A->rows; i++) {
        double sum = 0.0;
        for (int j = 0; j < A->cols; j++) {
            sum += A->m[i][j] * x->data[j];
        }
        result->data[i] = sum;
    }
}

void pds_matrix_print(const PDSMatrix* mat, const char* label) {
    if (!mat) { printf("%s: NULL\n", label ? label : "Matrix"); return; }
    printf("%s [%dx%d]:\n", label ? label : "Matrix", mat->rows, mat->cols);
    int max_rows = (mat->rows > 6) ? 6 : mat->rows;
    int max_cols = (mat->cols > 6) ? 6 : mat->cols;
    for (int i = 0; i < max_rows; i++) {
        printf("  ");
        for (int j = 0; j < max_cols; j++) {
            printf("%8.4f ", mat->m[i][j]);
        }
        if (mat->cols > 6) printf("...");
        printf("\n");
    }
    if (mat->rows > 6) printf("  ... (%d more rows)\n", mat->rows - 6);
}

void pds_matrix_copy(const PDSMatrix* src, PDSMatrix* dst) {
    if (!src || !dst) return;
    if (src->rows != dst->rows || src->cols != dst->cols) return;
    for (int i = 0; i < src->rows; i++) {
        memcpy(dst->m[i], src->m[i], src->cols * sizeof(double));
    }
}

void pds_matrix_add(PDSMatrix* A, const PDSMatrix* B, PDSMatrix* result) {
    if (!A || !B || !result) return;
    if (A->rows != B->rows || A->cols != B->cols) return;
    if (A->rows != result->rows || A->cols != result->cols) return;
    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->cols; j++) {
            result->m[i][j] = A->m[i][j] + B->m[i][j];
        }
    }
}

void pds_matrix_scale(PDSMatrix* A, double alpha, PDSMatrix* result) {
    if (!A || !result) return;
    if (A->rows != result->rows || A->cols != result->cols) return;
    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->cols; j++) {
            result->m[i][j] = A->m[i][j] * alpha;
        }
    }
}

double pds_matrix_determinant_2x2(const PDSMatrix* A) {
    if (!A || A->rows < 2 || A->cols < 2) return 0.0;
    return A->m[0][0] * A->m[1][1] - A->m[0][1] * A->m[1][0];
}

double pds_matrix_trace(const PDSMatrix* A) {
    if (!A) return 0.0;
    int n = (A->rows < A->cols) ? A->rows : A->cols;
    double tr = 0.0;
    for (int i = 0; i < n; i++) {
        tr += A->m[i][i];
    }
    return tr;
}

void pds_matrix_transpose(const PDSMatrix* A, PDSMatrix* result) {
    if (!A || !result) return;
    if (A->rows != result->cols || A->cols != result->rows) return;
    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->cols; j++) {
            result->m[j][i] = A->m[i][j];
        }
    }
}

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

double pds_sign(double x) {
    if (x > 0.0) return 1.0;
    if (x < 0.0) return -1.0;
    return 0.0;
}

double pds_clamp(double x, double lo, double hi) {
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

bool pds_is_near(double a, double b, double tol) {
    return fabs(a - b) <= tol;
}