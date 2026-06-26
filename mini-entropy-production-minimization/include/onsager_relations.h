/**
 * @file onsager_relations.h
 * @brief Onsager reciprocal relations and coefficient matrices.
 *
 * The Onsager reciprocal relations are the cornerstone of linear
 * irreversible thermodynamics. They state that the matrix of
 * phenomenological coefficients L_{ij} relating fluxes to forces
 * is symmetric:  L_{ij} = L_{ji}  (Onsager, 1931).
 *
 * For systems with magnetic fields B or angular velocity omega,
 * the generalized Casimir-Onsager relations hold:
 *   L_{ij}(B) = L_{ji}(-B)
 *
 * References:
 *   - Onsager, L. (1931) Phys. Rev. 37, 405; 38, 2265
 *   - Casimir, H.B.G. (1945) Rev. Mod. Phys. 17, 343
 *   - de Groot & Mazur (1984), Ch. IV-V
 */

#ifndef ONSAGER_RELATIONS_H
#define ONSAGER_RELATIONS_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define ONSAGER_MAX_DIM 64

/**
 * @brief Onsager phenomenological coefficient matrix.
 *
 * The matrix L_{ij} relates thermodynamic fluxes J_i to forces X_j
 * in the linear regime: J_i = sum_j L_{ij} X_j.
 *
 * Properties (per Onsager reciprocity and Second Law):
 *   1. Symmetry: L_{ij} = L_{ji}
 *   2. Positive semidefiniteness: all eigenvalues >= 0
 *   3. The diagonal elements L_{ii} are direct coefficients
 *   4. The off-diagonal elements are coupling coefficients
 *
 * Stored in row-major order: L[i * dim + j] = L_{ij}
 */
typedef struct {
    double *data;
    int     dim;
    int     is_symmetric;
    int     is_psd;
} OnsagerMatrix;

/**
 * @brief Resistance matrix R = L^{-1}.
 *
 * Expresses forces in terms of fluxes: X_i = sum_j R_{ij} J_j.
 * It is the inverse of the Onsager matrix and inherits symmetry
 * and positive definiteness from L.
 */
typedef struct {
    double *data;
    int     dim;
    int     is_symmetric;
    int     is_pd;
} ResistanceMatrix;

/**
 * @brief Coupling block structure for multi-process systems.
 *
 * For systems with multiple types of irreversible processes
 * (e.g., heat + diffusion + reaction), the Onsager matrix has
 * block structure where diagonal blocks represent direct effects
 * and off-diagonal blocks represent cross-effects (Soret, Dufour,
 * Seebeck, Peltier effects).
 */
typedef struct {
    int     num_blocks;
    int     block_sizes[8];
    double *block_matrices[8][8];
} OnsagerBlockStructure;

/* ---- L2: Core Concepts -- Matrix operations ---- */

OnsagerMatrix *onsager_alloc(int dim);
void onsager_free(OnsagerMatrix *L);
int  onsager_set(OnsagerMatrix *L, int i, int j, double value);
double onsager_get(const OnsagerMatrix *L, int i, int j);
int  onsager_init_diagonal(OnsagerMatrix *L, const double diag[]);
int  onsager_init_fourier(OnsagerMatrix *L, double thermal_cond, double temperature);
int  onsager_init_diffusion(OnsagerMatrix *L, const double *D, const double conc[], double T, int n);
int  onsager_init_two_process(OnsagerMatrix *L, double L11, double L22, double L12);
int  onsager_init_chemical(OnsagerMatrix *L, const double eq_rates[], const double *stoichiometry, int r, int s, double T, double R_gas);

/* ---- L3: Mathematical Structures -- Matrix verification ---- */

int    onsager_verify_reciprocity(OnsagerMatrix *L, double tol);
int    onsager_verify_psd(const OnsagerMatrix *L, double tol);
int    onsager_eigenvalues(const OnsagerMatrix *L, double eigenvalues[], int max_iter);
double onsager_condition_number(const OnsagerMatrix *L);

/* ---- L4: Fundamental Laws -- Onsager reciprocity theorem ---- */

double onsager_reciprocity_deviation(const OnsagerMatrix *L);
void   onsager_magnetic_reversal(OnsagerMatrix *L, const int parity[], int n);

/* ---- L5: Algorithms -- Matrix decomposition and solving ---- */

int onsager_ldlt_decompose(double *A, double *D, int n);
int onsager_solve(const double *L, const double *b, double *x, int n);
int onsager_invert(const double *L, double *R, int n);
int onsager_compute_fluxes(const OnsagerMatrix *L, const double forces[], double fluxes[]);
int onsager_compute_forces(const OnsagerMatrix *L, const double fluxes[], double forces[]);
int onsager_estimate_from_data(const double *forces, const double *fluxes, int m, int n, double *L);
void onsager_print(const OnsagerMatrix *L, const char *name);

#ifdef __cplusplus
}
#endif

#endif /* ONSAGER_RELATIONS_H */
