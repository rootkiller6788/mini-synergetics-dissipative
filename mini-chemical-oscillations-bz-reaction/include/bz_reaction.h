/**
 * @file bz_reaction.h
 * @brief Core BZ reaction definitions: Field-Koros-Noyes (FKN) mechanism,
 *        mass-action kinetics, and reaction network topology.
 *
 * Knowledge coverage:
 *   L1: BZ reaction definition, FKN mechanism, autocatalysis, oscillation
 *   L2: Positive/negative feedback loops, excitability, bistability
 *   L3: Stoichiometric matrix, rate-law ODE system, reaction network
 *
 * References:
 *   - Field, Koros, Noyes (1972) JACS 94, 8649-8664
 *   - Field & Noyes (1974) J. Chem. Phys. 60, 1877-1884
 *   - Scott (1994) "Oscillations, Waves, and Chaos in Chemical Kinetics"
 */

#ifndef BZ_REACTION_H
#define BZ_REACTION_H

#include <stddef.h>

/* ---------------------------------------------------------------------------
 * L1: Core definitions -- FKN mechanism species and rate constants
 * --------------------------------------------------------------------------- */

#define BZ_FKN_NUM_SPECIES  7

/** Indices for the 7 key FKN species */
typedef enum {
    BZ_SPECIES_BR    = 0,   /* Br-      bromide ion */
    BZ_SPECIES_BRO2  = 1,   /* BrO2.    bromine dioxide radical */
    BZ_SPECIES_BRO3  = 2,   /* BrO3-    bromate ion */
    BZ_SPECIES_HBRO2 = 3,   /* HBrO2    bromous acid */
    BZ_SPECIES_CE4   = 4,   /* Ce4+     cerium(IV) */
    BZ_SPECIES_CE3   = 5,   /* Ce3+     cerium(III) */
    BZ_SPECIES_MA    = 6    /* MA       malonic acid */
} bz_fkn_species_t;

#define BZ_FKN_NUM_REACTIONS  10

/** The three major processes of the FKN mechanism */
typedef enum {
    BZ_PROCESS_A = 0,  /* Br- consumption process */
    BZ_PROCESS_B = 1,  /* Autocatalytic HBrO2 production */
    BZ_PROCESS_C = 2   /* Ce4+ reduction by malonic acid */
} bz_fkn_process_t;

/* ---------------------------------------------------------------------------
 * L1: Core state definitions
 * --------------------------------------------------------------------------- */

/** Concentrations of all relevant BZ species (mol/L) */
typedef struct {
    double br;       /* [Br-]      bromide */
    double bro2;     /* [BrO2.]    bromine dioxide radical */
    double bro3;     /* [BrO3-]    bromate */
    double hbro2;    /* [HBrO2]    bromous acid */
    double ce4;      /* [Ce4+]     cerium(IV) */
    double ce3;      /* [Ce3+]     cerium(III) */
    double ma;       /* [MA]       malonic acid */
    double h;        /* [H+]       proton (held constant) */
    double hobr;     /* [HOBr]     hypobromous acid intermediate */
    double br2;      /* [Br2]      bromine intermediate */
} bz_state_t;

/** BZ reaction parameters (rate constants, T=25 C unless noted) */
typedef struct {
    double k_a1;     /* BrO3- + Br- + 2H+ -> HBrO2 + HOBr */
    double k_a2;     /* HBrO2 + Br- + H+ -> 2HOBr */
    double k_b1;     /* BrO3- + HBrO2 + H+ -> 2BrO2. + H2O */
    double k_b2;     /* BrO2. + Ce3+ + H+ -> HBrO2 + Ce4+ */
    double k_b3;     /* 2HBrO2 -> BrO3- + HOBr + H+ */
    double k_c1;     /* Ce4+ + MA -> products */
    double k_c2;     /* MA + Br2 -> BrMA + Br- + H+ */
    double k_c3;     /* BrMA + Ce4+ -> Ce3+ + Br- + products */
    double f;        /* Stoichiometric factor, typically 0.5 to 3.0 */
    double q;        /* Rate constant ratio */
    double epsilon;  /* Time-scale separation parameter */
    double epsilon_p;/* Second time-scale separation parameter */
} bz_params_t;

/* ---------------------------------------------------------------------------
 * L2: Core concepts -- dynamical regimes of the BZ reaction
 * --------------------------------------------------------------------------- */

typedef enum {
    BZ_MODE_STEADY_STATE    = 0,  /* Approach equilibrium, no oscillation */
    BZ_MODE_OSCILLATORY     = 1,  /* Sustained limit-cycle oscillation */
    BZ_MODE_EXCITABLE       = 2,  /* Quiescent but excitable */
    BZ_MODE_BISTABLE        = 3,  /* Two stable steady states coexist */
    BZ_MODE_CHAOTIC         = 4,  /* Deterministic chaos */
    BZ_MODE_BURSTING        = 5   /* Periodic bursting */
} bz_mode_t;

/** Feedback loop type in the reaction network */
typedef enum {
    BZ_FEEDBACK_POSITIVE = 0,  /* Autocatalysis: HBrO2 -> 2 HBrO2 */
    BZ_FEEDBACK_NEGATIVE = 1,  /* Inhibition: Br- production */
    BZ_FEEDBACK_DELAYED  = 2   /* Delayed feedback from slow process */
} bz_feedback_type_t;

/* ---------------------------------------------------------------------------
 * L3: Mathematical structures -- stoichiometry and reaction network
 * --------------------------------------------------------------------------- */

/** Stoichiometric matrix: S[i][j] = change of species i in reaction j */
typedef struct {
    int  n_species;
    int  n_reactions;
    int *data;       /* Row-major layout */
} bz_stoichiometric_t;

/** Reaction rate vector for mass-action kinetics */
typedef struct {
    double *rates;   /* rates[j] = k_j * Prod_i [X_i]^{order_{i,j}} */
    int     n_reactions;
} bz_rate_vector_t;

/** BZ reaction network adjacency representation */
typedef struct {
    int     n_species;
    int     n_edges;
    int    *from;    /* Source species index */
    int    *to;      /* Target species index */
    double *weight;  /* Reaction flux magnitude */
} bz_reaction_network_t;

/* ---------------------------------------------------------------------------
 * L3: Reduced models -- Oregonator (Field-Noyes 1974)
 * --------------------------------------------------------------------------- */

/** Oregonator 3-variable state (dimensionless variables) */
typedef struct {
    double x;        /* [HBrO2] -- autocatalytic species, fast */
    double y;        /* [Br-]   -- control intermediate, medium */
    double z;        /* [Ce4+]  -- catalyst redox state, slow */
} bz_oregonator_state_t;

/** Oregonator dimensionless parameters */
typedef struct {
    double f;        /* Stoichiometric factor */
    double q;        /* Rate constant ratio */
    double epsilon;  /* Time-scale ratio, ~1e-2 */
    double epsilon_p;/* Second time-scale ratio, ~1e-4 */
} bz_oregonator_params_t;

/* ---------------------------------------------------------------------------
 * API function declarations
 * --------------------------------------------------------------------------- */

/* Initialization */

void bz_state_init_typical(bz_state_t *s);

void bz_oregonator_init_oscillatory(bz_oregonator_state_t *s,
                                     const bz_oregonator_params_t *p);

void bz_params_init_default(bz_params_t *p);

void bz_oregonator_params_init(bz_oregonator_params_t *p, double f);

/* Rate functions */

double bz_fkn_rhs(bz_state_t *dsdt, const bz_state_t *s, const bz_params_t *p);

void bz_oregonator_rhs(bz_oregonator_state_t *dsdt,
                        const bz_oregonator_state_t *s,
                        const bz_oregonator_params_t *p);

void bz_tyson_2var_rhs(double *dudt, double *dvdt,
                        double u, double v,
                        const bz_oregonator_params_t *p);

/* Steady-state analysis */

int bz_oregonator_steady_states(bz_oregonator_state_t *ss,
                                 const bz_oregonator_params_t *p);

void bz_oregonator_nullclines(double *y_x_nullcline,
                               double *y_y_nullcline,
                               double *z_nullcline,
                               double x, double z_val,
                               const bz_oregonator_params_t *p);

#endif /* BZ_REACTION_H */
