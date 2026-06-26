#ifndef NEPT_BIFURCATION_H
#define NEPT_BIFURCATION_H

#include <stdbool.h>
#include "phase_transition.h"

/* ============================================================================
 * Bifurcation Theory as Non-Equilibrium Phase Transitions
 *
 * In nonlinear dynamics, a bifurcation is a qualitative change in the
 * behavior of a dynamical system as a control parameter crosses a critical
 * value. This is the dynamical systems analog of a phase transition.
 *
 * The Landau analogy: the bifurcation parameter mu plays the role of
 * reduced temperature t, and the amplitude of the emerging mode plays
 * the role of the order parameter.
 *
 * Normal forms capture the universal structure near bifurcation points:
 *   Saddle-node:  dx/dt = mu + a*x^2       (discontinuous, hysteresis)
 *   Transcritical: dx/dt = mu*x - x^2      (exchange of stability)
 *   Pitchfork:    dx/dt = mu*x - b*x^3     (b>0: supercritical = continuous)
 *   Hopf:         dz/dt = (mu + i*omega)*z - (1+i*c)*|z|^2*z
 *
 * References:
 *   Yu.A. Kuznetsov — Elements of Applied Bifurcation Theory (2004)
 *   J. Guckenheimer & P. Holmes — Nonlinear Oscillations... (1983)
 *   S.H. Strogatz — Nonlinear Dynamics and Chaos (2015)
 *   H. Haken — Synergetics (1977)
 * ============================================================================ */

/* --- Bifurcation Types --- */
typedef enum {
    BIF_NONE = 0,
    BIF_SADDLE_NODE = 1,        /* Fold: creation/annihilation of fixed points */
    BIF_TRANSCRITICAL = 2,      /* Exchange of stability between two branches   */
    BIF_PITCHFORK_SUPER = 3,    /* Supercritical pitchfork (continuous)         */
    BIF_PITCHFORK_SUB = 4,     /* Subcritical pitchfork (discontinuous)         */
    BIF_HOPF_SUPER = 5,        /* Supercritical Hopf: stable limit cycle        */
    BIF_HOPF_SUB = 6,         /* Subcritical Hopf: unstable limit cycle         */
    BIF_HOMOCLINIC = 7,        /* Homoclinic bifurcation (global)               */
    BIF_HETEROCLINIC = 8,      /* Heteroclinic cycle bifurcation                */
    BIF_PERIOD_DOUBLING = 9,   /* Flip: route to chaos                         */
    BIF_TORUS = 10,            /* Neimark-Sacker: 2-torus birth                 */
    BIF_SNIPER = 11            /* Saddle-node on invariant circle               */
} BifurcationType;

/* --- Normal Form Structures --- */

/* Saddle-node: dx/dt = mu + a*x^2 */
typedef struct {
    double mu;                  /* Control parameter (bifurcation parameter)     */
    double a;                   /* Quadratic coefficient (a != 0)                */
    double x_eq_stable;         /* Stable branch: x = -sqrt(-mu/a) for mu/a < 0 */
    double x_eq_unstable;       /* Unstable branch: x = +sqrt(-mu/a)            */
    bool folded_exists;         /* Whether real equilibria exist (mu/a <= 0)    */
} SaddleNodeNF;

/* Pitchfork: dx/dt = mu*x - b*x^3 (b > 0 = supercritical, b < 0 = subcritical) */
typedef struct {
    double mu;                  /* Bifurcation parameter                         */
    double b;                   /* Cubic coefficient                             */
    bool is_supercritical;      /* b > 0 => supercritical (continuous)           */
    double x_eq_trivial;        /* x = 0 branch                                 */
    double x_eq_nontrivial;     /* x = +/- sqrt(mu/b) for mu*b > 0              */
} PitchforkNF;

/* Hopf: dz/dt = (mu + i*omega)*z - (c1 + i*c2)*|z|^2*z */
typedef struct {
    double mu;                  /* Real part of linear eigenvalue                */
    double omega;               /* Imaginary part (oscillation frequency)        */
    double c1;                  /* Real part of cubic (c1 > 0 => supercritical) */
    double c2;                  /* Imag part of cubic (frequency shift)          */
    double radius_equilibrium;  /* |z_eq| = sqrt(mu/c1) for mu*c1 > 0           */
    double limit_cycle_freq;    /* omega_eq = omega + c2*|z_eq|^2               */
    bool is_supercritical;      /* c1 > 0                                       */
} HopfNF;

/* --- Bifurcation Diagram Data --- */
#define BIF_MAX_BRANCHES 16
#define BIF_MAX_POINTS_PER_BRANCH 512

typedef struct {
    double *control_values;
    double *state_values;
    int n_points;
    bool is_stable;
    int branch_id;
} BifurcationBranch;

typedef struct {
    BifurcationBranch branches[BIF_MAX_BRANCHES];
    int n_branches;
    double control_min, control_max;
    BifurcationType bif_type;
    double critical_parameter;
} BifurcationDiagram;

/* --- Linear Stability Analysis --- */
typedef struct {
    int dimension;              /* State space dimension                         */
    double *jacobian;           /* Jacobian matrix at fixed point (dim x dim)   */
    double *eigenvalues_real;   /* Real parts of eigenvalues                    */
    double *eigenvalues_imag;   /* Imaginary parts of eigenvalues               */
    double max_real_part;       /* Leading eigenvalue real part (stability)     */
    double max_imag_part;       /* Leading eigenvalue imaginary part            */
    bool is_stable;             /* True if all Re(lambda) < 0                   */
    bool has_hopf_candidate;    /* True if complex pair crosses imag axis       */
} LinearStability;

/* --- Centre Manifold Reduction --- */
typedef struct {
    int dim_centre;             /* Dimension of the centre eigenspace           */
    double *centre_eigenvectors;/* Basis for centre manifold                    */
    double *reduced_jacobian;   /* Jacobian restricted to centre manifold       */
    double *normal_form_coeffs; /* Coefficients of the reduced normal form      */
} CentreManifold;

/* --- Core API --- */

/* Saddle-node normal form evaluation */
double nept_sn_eval(const SaddleNodeNF *nf, double x);

/* Compute saddle-node equilibrium branches */
void nept_sn_equilibria(SaddleNodeNF *nf);

/* Pitchfork normal form evaluation */
double nept_pitchfork_eval(const PitchforkNF *nf, double x);

/* Compute pitchfork equilibrium branches */
void nept_pitchfork_equilibria(PitchforkNF *nf);

/* Hopf normal form evaluation (returns d|z|/dt) */
double nept_hopf_amplitude_ode(const HopfNF *nf, double radius);

/* Compute Hopf limit cycle amplitude and frequency */
void nept_hopf_cycle_properties(const HopfNF *nf,
                                 double *amplitude_out,
                                 double *frequency_out);

/* Classify bifurcation type from eigenvalue crossing pattern */
BifurcationType nept_classify_bifurcation(double real_part_before,
                                           double real_part_after,
                                           double imag_part,
                                           double nonlinear_coeff);

/* Perform linear stability analysis of a fixed point */
LinearStability *nept_linear_stability(const double *jacobian, int dim);
void nept_linear_stability_free(LinearStability *ls);

/* Create a bifurcation diagram structure */
BifurcationDiagram *nept_bifdiagram_create(double ctrl_min, double ctrl_max);
void nept_bifdiagram_free(BifurcationDiagram *bd);

/* Add a computed branch to the diagram */
int nept_bifdiagram_add_branch(BifurcationDiagram *bd,
                                const double *ctrl, const double *state,
                                int n_points, bool is_stable);

/* Compute normal form coefficients from explicit ODE at bifurcation point */
CentreManifold *nept_centre_manifold_reduce(const double *jacobian, int dim,
                                              const double *nonlinear_terms,
                                              int n_nonlinear);
void nept_centre_manifold_free(CentreManifold *cm);

/* Determine if a bifurcation is supercritical or subcritical from cubic term */
bool nept_is_supercritical(double cubic_coefficient);

/* Compute the leading Lyapunov exponent via variational equation (discrete) */
double nept_max_lyapunov_exponent(
    void (*ode)(double t, const double *x, double *dxdt, void *params),
    const double *x0, int dim, void *params,
    double dt, int n_steps);

/* Compute 1D bifurcation diagram by numerical continuation (simple) */
int nept_continue_branch(
    double (*f)(double x, double mu),
    double mu_start, double mu_end, int n_steps,
    double x_guess,
    double *ctrl_out, double *state_out);

/* Check for eigenvalue crossing (simple 1D case) */
bool nept_eigenvalue_crossing(double lambda_real, double tolerance);

#endif /* NEPT_BIFURCATION_H */