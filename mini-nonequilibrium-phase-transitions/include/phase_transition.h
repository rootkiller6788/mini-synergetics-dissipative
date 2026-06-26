#ifndef NEPT_PHASE_TRANSITION_H
#define NEPT_PHASE_TRANSITION_H

#include <stdbool.h>
#include <stddef.h>

/* ============================================================================
 * Non-Equilibrium Phase Transitions — Core Definitions
 *
 * References:
 *   H. Haken  — Synergetics: An Introduction (1977)
 *   I. Prigogine — Self-Organization in Nonequilibrium Systems (1977)
 *   G. Nicolis & I. Prigogine — Self-Organization in Nonequilibrium Systems (1977)
 *   L.D. Landau & E.M. Lifshitz — Statistical Physics, Part 1 (1980)
 *   H.E. Stanley — Introduction to Phase Transitions and Critical Phenomena (1971)
 *   M.C. Cross & P.C. Hohenberg — Pattern Formation Outside Equilibrium (1993)
 *   N.G. van Kampen — Stochastic Processes in Physics and Chemistry (2007)
 *   C.W. Gardiner — Handbook of Stochastic Methods (2004)
 *   P.M. Chaikin & T.C. Lubensky — Principles of Condensed Matter Physics (1995)
 *   S.K. Ma — Modern Theory of Critical Phenomena (1976)
 * ============================================================================ */

/* --- Phase Transition Type --- */
typedef enum {
    NEPT_FIRST_ORDER = 0,       /* Discontinuous: latent heat, hysteresis      */
    NEPT_CONTINUOUS = 1,        /* Second-order: continuous order parameter     */
    NEPT_INFINITE_ORDER = 2,    /* Berezinskii-Kosterlitz-Thouless transition   */
    NEPT_QUANTUM = 3,           /* Driven by quantum fluctuations at T=0        */
    NEPT_DYNAMICAL = 4,         /* Transient behavior, non-steady criticality   */
    NEPT_NOISE_INDUCED = 5      /* Purely fluctuation-driven transition         */
} NEPT_TransitionType;

/* --- System Classification (far from equilibrium) --- */
typedef enum {
    NEPT_EQUILIBRIUM = 0,       /* Boltzmann-Gibbs thermal equilibrium          */
    NEPT_NEAR_EQUILIBRIUM = 1,  /* Linear regime, Onsager relations hold        */
    NEPT_FAR_EQUILIBRIUM = 2,   /* Nonlinear regime, dissipative structures     */
    NEPT_DRIVEN_STEADY = 3,     /* Constant external driving, steady state NESS */
    NEPT_QUENCHED = 4           /* Rapidly quenched, glassy dynamics            */
} NEPT_SystemRegime;

/* --- Order of the phase transition (Ehrenfest classification) --- */
typedef enum {
    NEPT_EHRENFEST_NONE = 0,    /* No singularity in free energy derivatives   */
    NEPT_EHRENFEST_FIRST = 1,   /* Discontinuity in first derivative of F       */
    NEPT_EHRENFEST_SECOND = 2,  /* Discontinuity in second derivative (chi, Cv)   */
    NEPT_EHRENFEST_THIRD = 3    /* Higher-order singularities                    */
} NEPT_EhrenfestClass;

/* --- Symmetry Breaking Pattern --- */
typedef enum {
    NEPT_SYMMETRY_NONE = 0,        /* No symmetry change                  */
    NEPT_SYMMETRY_DISCRETE_Z2 = 1, /* Ising-like Z2 breaking              */
    NEPT_SYMMETRY_DISCRETE_ZN = 2, /* Potts-like ZN breaking              */
    NEPT_SYMMETRY_CONTINUOUS = 3,  /* O(N) rotational breaking (Heisenberg, XY) */
    NEPT_SYMMETRY_TRANSLATIONAL = 4, /* Crystal, pattern, stripe formation */
    NEPT_SYMMETRY_GAUGE = 5        /* Gauge symmetry (superconductors etc.) */
} NEPT_SymmetryType;

/* --- Control Parameter Structures --- */
typedef struct {
    double value;               /* Current value of control parameter     */
    double critical;            /* Critical value at which transition occurs */
    double sweep_rate;          /* Rate of change (for quench dynamics)   */
    const char *name;           /* Human-readable identifier              */
} ControlParameter;

/* --- Phase Region Descriptor --- */
typedef struct {
    int region_index;              /* Ordinal label within phase diagram  */
    const char *phase_name;        /* Descriptive name (e.g. "disordered") */
    NEPT_SymmetryType symmetry;    /* Symmetry in this phase              */
    bool is_stable;                /* Thermodynamic/dynamic stability     */
    double free_energy_density;    /* Free energy per unit volume         */
} PhaseRegion;

/* --- Bifurcation Point as Phase Transition --- */
typedef struct {
    double control_value;             /* Critical value of the control parameter */
    NEPT_TransitionType trans_type;   /* Type of transition at this point       */
    int num_stable_branches_before;   /* Number of stable steady states before   */
    int num_stable_branches_after;    /* Number of stable steady states after    */
    double order_jump;                /* Magnitude of discontinuity (1st-order)  */
} TransitionPoint;

/* --- Phase Diagram — the central organizing structure --- */
#define NEPT_MAX_REGIONS 32
#define NEPT_MAX_TRANSITIONS 64

typedef struct {
    PhaseRegion regions[NEPT_MAX_REGIONS];
    int num_regions;

    TransitionPoint transitions[NEPT_MAX_TRANSITIONS];
    int num_transitions;

    /* Parameter axes */
    ControlParameter param_x;         /* Primary control parameter        */
    ControlParameter param_y;         /* Secondary (for 2D phase diagram) */
    bool is_two_dimensional;

    /* Global properties */
    double tricritical_point;         /* Location if tricritical exists   */
    bool has_critical_endpoint;       /* Whether a critical point exists  */
    double critical_endpoint_value;
} PhaseDiagram;

/* --- Non-Equilibrium Steady State (NESS) --- */
typedef struct {
    double *state_variables;       /* Array of state variable values        */
    int num_variables;
    double entropy_production;     /* sigma = d_iS/dt >= 0 (Prigogine)     */
    double entropy_flux;           /* d_eS/dt exchange with environment     */
    double lyapunov_time;          /* Characteristic time scale             */
    bool is_dissipative_structure; /* True if spatial/temporal order exists */
} NonEquilibriumSteadyState;

/* --- Transition Observation / Measurement --- */
typedef struct {
    double order_parameter_value;  /* Measured OP at this snapshot         */
    double susceptibility;         /* Response function chi                */
    double correlation_length;     /* xi, diverges at criticality          */
    double time;                   /* Observation time                     */
    double control_parameter;      /* Control parameter at measurement     */
} TransitionObservation;

typedef struct {
    TransitionObservation *observations;
    int num_observations;
    int capacity;
    double critical_point_estimate;
    double critical_exponent_estimate;
} TransitionSeries;

/* --- Core API --- */

/* Initialize a phase diagram with a single parameter axis */
PhaseDiagram nept_phasediagram_create(const char *param_name,
                                       double param_critical);

/* Add a phase region to the diagram */
int nept_phasediagram_add_region(PhaseDiagram *pd,
                                  const char *name,
                                  NEPT_SymmetryType sym,
                                  double free_energy);

/* Add a transition point between phases */
int nept_phasediagram_add_transition(PhaseDiagram *pd,
                                      double control_val,
                                      NEPT_TransitionType type,
                                      int branches_before,
                                      int branches_after,
                                      double jump);

/* Classify the transition by analyzing order parameter data */
NEPT_TransitionType nept_classify_transition(const double *order_param,
                                              const double *control,
                                              int n_points,
                                              double *jump_out);

/* Determine Ehrenfest class from free energy derivatives */
NEPT_EhrenfestClass nept_ehrenfest_classify(double dF1_below,
                                             double dF1_above,
                                             double dF2_below,
                                             double dF2_above,
                                             double tolerance);

/* Compute the entropy production for a given steady state */
double nept_entropy_production(const double *fluxes,
                                const double *forces,
                                int n_pairs);

/* Check Prigogine minimum entropy production principle
 * (valid in the linear near-equilibrium regime) */
bool nept_prigogine_min_entropy_check(const double *fluxes,
                                       const double *forces,
                                       int n_pairs,
                                       double *entropy_out);

/* Glansdorff-Prigogine evolution criterion (nonlinear regime) */
double nept_glansdorff_prigogine_criterion(const double *state,
                                            const double *state_deriv,
                                            const double *flux_deriv,
                                            int n_vars);

/* Determine system regime from dimensionless numbers */
NEPT_SystemRegime nept_regime_classify(double reynolds_or_analogue,
                                        double critical_value);

/* Allocate and initialize a transition observation series */
TransitionSeries *nept_series_create(int capacity);

/* Push an observation to the series */
int nept_series_push(TransitionSeries *ts,
                      double order_param,
                      double susceptibility,
                      double corr_length,
                      double time,
                      double control);

/* Free the transition series */
void nept_series_free(TransitionSeries *ts);

/* Estimate the critical point from the series */
double nept_series_estimate_critical_point(const TransitionSeries *ts);

/* --- NESS API --- */
NonEquilibriumSteadyState *nept_ness_create(int num_variables);
void nept_ness_free(NonEquilibriumSteadyState *ness);
void nept_ness_set_state(NonEquilibriumSteadyState *ness,
                          const double *values);
bool nept_ness_is_dissipative(const NonEquilibriumSteadyState *ness,
                               double threshold);

#endif /* NEPT_PHASE_TRANSITION_H */