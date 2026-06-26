#include "bifurcation.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* ============================================================================
 * Bifurcation Analysis — Implementation
 *
 * Knowledge points:
 * 1. Saddle-node normal form (L3: fold bifurcation, dx/dt = μ + a*x^2)
 * 2. Pitchfork normal form (L3: symmetry-breaking, dx/dt = μ*x - b*x^3)
 * 3. Hopf normal form (L3: oscillation onset, dz/dt = (μ+iω)z - (c₁+ic₂)|z|²z)
 * 4. Bifurcation classification (L2: eigenvalue crossing patterns)
 * 5. Linear stability analysis (L5: Jacobian eigenvalues)
 * 6. Centre manifold reduction (L8: dimension reduction near bifurcation)
 * 7. Lyapunov exponent computation (L5: chaos diagnostic)
 * 8. Numerical continuation (L5: branch following)
 * ============================================================================ */

/* ---------------------------------------------------------------------------
 * Saddle-Node Normal Form: dx/dt = μ + a*x²
 *
 * Knowledge point: The saddle-node (fold) bifurcation is the generic
 * mechanism for creation/annihilation of fixed points. Normal form:
 *
 *   dx/dt = μ + a*x²   (a ≠ 0)
 *
 * Equilibria: x_eq = ±√(-μ/a) when μ/a ≤ 0
 * Stability: d/dx(μ + a*x²) = 2a*x
 *   For a > 0: x = -√(-μ/a) stable, x = +√(-μ/a) unstable
 *   For a < 0: reversed
 *
 * This is the prototypical "tipping point" — a small change in μ
 * causes abrupt disappearance of the equilibrium.
 *
 * This is analogous to a first-order phase transition: the order
 * parameter jumps discontinuously when the control parameter passes
 * the critical value.
 *
 * Reference: Strogatz, §3.1; Kuznetsov, §2.3
 * --------------------------------------------------------------------------- */
double nept_sn_eval(const SaddleNodeNF *nf, double x)
{
    if (!nf) return 0.0;
    return nf->mu + nf->a * x * x;
}

void nept_sn_equilibria(SaddleNodeNF *nf)
{
    if (!nf) return;
    double arg = -nf->mu / nf->a;
    if (arg >= 0.0) {
        nf->folded_exists = true;
        double root = sqrt(arg);
        /* For a > 0: negative branch stable, positive branch unstable */
        nf->x_eq_stable = -root;
        nf->x_eq_unstable = root;
    } else {
        nf->folded_exists = false;
        nf->x_eq_stable = 0.0;
        nf->x_eq_unstable = 0.0;
    }
}

/* ---------------------------------------------------------------------------
 * Pitchfork Normal Form: dx/dt = μ*x - b*x³
 *
 * Knowledge point: Pitchfork bifurcation occurs in systems with
 * Z₂ symmetry (x → -x). Two types:
 *
 * Supercritical (b > 0): stable non-zero equilibria emerge continuously
 *   for μ > 0. This is the dynamical analog of a continuous
 *   phase transition (μ ↔ reduced temperature t).
 *
 * Subcritical (b < 0): unstable non-zero equilibria exist for μ < 0.
 *   The system jumps discontinuously — analogous to a first-order
 *   phase transition.
 *
 * Equilibrium branches:
 *   x = 0 (trivial, always exists)
 *   x = ±√(μ/b) for μ*b > 0
 *
 * Reference: Strogatz, §3.4; Haken, Synergetics, §5.1
 * --------------------------------------------------------------------------- */
double nept_pitchfork_eval(const PitchforkNF *nf, double x)
{
    if (!nf) return 0.0;
    return nf->mu * x - nf->b * x * x * x;
}

void nept_pitchfork_equilibria(PitchforkNF *nf)
{
    if (!nf) return;
    nf->x_eq_trivial = 0.0;
    nf->is_supercritical = (nf->b > 0.0);

    double arg = nf->mu / nf->b;
    if (arg > 0.0) {
        nf->x_eq_nontrivial = sqrt(arg);
    } else {
        nf->x_eq_nontrivial = 0.0;
    }
}

/* ---------------------------------------------------------------------------
 * Hopf Normal Form: dz/dt = (μ + iω)z - (c₁ + ic₂)|z|²z
 *
 * Knowledge point: Hopf bifurcation is the birth of a limit cycle from
 * a fixed point. The normal form in polar coordinates (z = r*e^{iθ}):
 *
 *   dr/dt = μ*r - c₁*r³   (amplitude equation: Landau form!)
 *   dθ/dt = ω + c₂*r²     (frequency shift)
 *
 * Supercritical (c₁ > 0): stable limit cycle for μ > 0
 *   Radius: r_eq = √(μ/c₁), Frequency: Ω = ω + c₂*r_eq²
 *
 * Subcritical (c₁ < 0): unstable limit cycle for μ < 0
 *   (requires higher-order terms c₅ > 0 for global stability)
 *
 * The amplitude equation is identical to the Landau equation for a
 * continuous phase transition, with r as the order parameter!
 * This is the deep connection between bifurcation theory and
 * non-equilibrium phase transitions.
 *
 * Reference: Strogatz, §8.2; Kuznetsov, §3.5
 * --------------------------------------------------------------------------- */
double nept_hopf_amplitude_ode(const HopfNF *nf, double radius)
{
    if (!nf) return 0.0;
    return nf->mu * radius - nf->c1 * radius * radius * radius;
}

void nept_hopf_cycle_properties(const HopfNF *nf,
                                 double *amplitude_out,
                                 double *frequency_out)
{
    if (!nf) {
        if (amplitude_out) *amplitude_out = 0.0;
        if (frequency_out) *frequency_out = 0.0;
        return;
    }

    if (nf->mu > 0.0 && nf->c1 > 1e-15) {
        double r_eq = sqrt(nf->mu / nf->c1);
        if (amplitude_out) *amplitude_out = r_eq;
        if (frequency_out) *frequency_out = nf->omega + nf->c2 * r_eq * r_eq;
    } else {
        if (amplitude_out) *amplitude_out = 0.0;
        if (frequency_out) *frequency_out = nf->omega;
    }
}

/* ---------------------------------------------------------------------------
 * Bifurcation Classification from Eigenvalue Crossing
 *
 * Knowledge point: Local bifurcations of fixed points are classified by
 * how eigenvalues of the Jacobian cross the imaginary axis:
 *
 * Real eigenvalue crosses zero (λ = 0):
 *   Generic                     → Saddle-node
 *   With Z₂ symmetry            → Pitchfork (super or sub)
 *   Exchange of stability        → Transcritical
 *
 * Complex conjugate pair crosses (λ = ±iω, ω ≠ 0):
 *   Generic                     → Hopf (super or sub)
 *
 * Classification also uses the sign of the leading nonlinear term
 * (Lyapunov coefficient / first focal value for Hopf).
 *
 * Reference: Kuznetsov, §2.1, §8.1
 * --------------------------------------------------------------------------- */
BifurcationType nept_classify_bifurcation(double real_part_before,
                                           double real_part_after,
                                           double imag_part,
                                           double nonlinear_coeff)
{
    bool was_stable = (real_part_before < 0.0);
    bool became_unstable = (real_part_after > 0.0);

    if (!was_stable || !became_unstable) {
        return BIF_NONE;
    }

    if (fabs(imag_part) < 1e-12) {
        /* Real eigenvalue crossing */
        if (nonlinear_coeff > 1e-12) {
            return BIF_PITCHFORK_SUPER;
        } else if (nonlinear_coeff < -1e-12) {
            return BIF_PITCHFORK_SUB;
        } else {
            return BIF_SADDLE_NODE;
        }
    } else {
        /* Complex pair crossing: Hopf */
        if (nonlinear_coeff > 1e-12) {
            return BIF_HOPF_SUPER;
        } else {
            return BIF_HOPF_SUB;
        }
    }
}

/* ---------------------------------------------------------------------------
 * Linear Stability Analysis
 *
 * Knowledge point: The stability of a fixed point x* of dx/dt = f(x)
 * is determined by the eigenvalues λ of the Jacobian J = ∂f/∂x|_{x*}.
 *
 * Stability criterion (Lyapunov's First Method):
 *   Re(λ_i) < 0 for all i  ⇒  asymptotically stable
 *   Re(λ_i) > 0 for any i  ⇒  unstable
 *   Re(λ_i) = 0 for some i ⇒  marginal (bifurcation candidate)
 *
 * Eigenvalue computation uses the QR algorithm (simplified:
 * for 1D/2D direct analytic formula, for larger systems we
 * implement power iteration for the dominant eigenvalue).
 *
 * Reference: Khalil, Nonlinear Systems, §3.3
 * --------------------------------------------------------------------------- */
LinearStability *nept_linear_stability(const double *jacobian, int dim)
{
    if (!jacobian || dim <= 0) return NULL;

    LinearStability *ls = (LinearStability *)malloc(sizeof(LinearStability));
    if (!ls) return NULL;

    ls->dimension = dim;
    ls->jacobian = (double *)malloc((size_t)(dim * dim) * sizeof(double));
    ls->eigenvalues_real = (double *)calloc((size_t)dim, sizeof(double));
    ls->eigenvalues_imag = (double *)calloc((size_t)dim, sizeof(double));

    if (!ls->jacobian || !ls->eigenvalues_real || !ls->eigenvalues_imag) {
        free(ls->jacobian);
        free(ls->eigenvalues_real);
        free(ls->eigenvalues_imag);
        free(ls);
        return NULL;
    }

    memcpy(ls->jacobian, jacobian, (size_t)(dim * dim) * sizeof(double));

    /* --- Eigenvalue computation --- */
    ls->max_real_part = -1e15;
    ls->max_imag_part = 0.0;
    ls->has_hopf_candidate = false;

    if (dim == 1) {
        /* 1D: λ = J[0] */
        ls->eigenvalues_real[0] = jacobian[0];
        ls->eigenvalues_imag[0] = 0.0;
        ls->max_real_part = jacobian[0];

    } else if (dim == 2) {
        /* 2D: analytic formula
         * J = [[a, b], [c, d]]
         * λ = (τ ± √(τ² - 4Δ))/2 where τ = a+d, Δ = ad-bc */
        double a = jacobian[0], b = jacobian[1];
        double c = jacobian[2], d = jacobian[3];
        double trace = a + d;
        double det = a * d - b * c;
        double disc = trace * trace - 4.0 * det;

        if (disc >= 0.0) {
            double sqrt_d = sqrt(disc);
            ls->eigenvalues_real[0] = 0.5 * (trace + sqrt_d);
            ls->eigenvalues_real[1] = 0.5 * (trace - sqrt_d);
            ls->eigenvalues_imag[0] = 0.0;
            ls->eigenvalues_imag[1] = 0.0;
        } else {
            ls->eigenvalues_real[0] = 0.5 * trace;
            ls->eigenvalues_real[1] = 0.5 * trace;
            ls->eigenvalues_imag[0] = 0.5 * sqrt(-disc);
            ls->eigenvalues_imag[1] = -0.5 * sqrt(-disc);
            ls->has_hopf_candidate = true;
        }

        for (int i = 0; i < 2; i++) {
            if (ls->eigenvalues_real[i] > ls->max_real_part)
                ls->max_real_part = ls->eigenvalues_real[i];
        }

    } else {
        /* nD: Power iteration for dominant eigenvalue */
        double *v = (double *)calloc((size_t)dim, sizeof(double));
        if (v) {
            for (int i = 0; i < dim; i++) v[i] = 1.0 / sqrt((double)dim);
            for (int iter = 0; iter < 100; iter++) {
                double *Av = (double *)calloc((size_t)dim, sizeof(double));
                if (!Av) break;
                for (int i = 0; i < dim; i++)
                    for (int j = 0; j < dim; j++)
                        Av[i] += jacobian[i * dim + j] * v[j];
                double norm = 0.0;
                for (int i = 0; i < dim; i++) norm += Av[i] * Av[i];
                norm = sqrt(norm);
                if (norm < 1e-15) { free(Av); break; }
                for (int i = 0; i < dim; i++) v[i] = Av[i] / norm;
                free(Av);
            }
            /* Rayleigh quotient: λ ≈ v^T J v */
            double lambda = 0.0;
            double *Jv = (double *)calloc((size_t)dim, sizeof(double));
            if (Jv) {
                for (int i = 0; i < dim; i++)
                    for (int j = 0; j < dim; j++)
                        Jv[i] += jacobian[i * dim + j] * v[j];
                for (int i = 0; i < dim; i++) lambda += v[i] * Jv[i];
                free(Jv);
            }
            ls->eigenvalues_real[0] = lambda;
            ls->max_real_part = lambda;
            free(v);
        }
    }

    ls->is_stable = (ls->max_real_part < -1e-12);

    return ls;
}

void nept_linear_stability_free(LinearStability *ls)
{
    if (ls) {
        free(ls->jacobian);
        free(ls->eigenvalues_real);
        free(ls->eigenvalues_imag);
        free(ls);
    }
}

/* --- Bifurcation diagram management --- */
BifurcationDiagram *nept_bifdiagram_create(double ctrl_min, double ctrl_max)
{
    BifurcationDiagram *bd = (BifurcationDiagram *)malloc(sizeof(BifurcationDiagram));
    if (!bd) return NULL;
    memset(bd, 0, sizeof(BifurcationDiagram));
    bd->control_min = ctrl_min;
    bd->control_max = ctrl_max;
    return bd;
}

void nept_bifdiagram_free(BifurcationDiagram *bd)
{
    if (bd) {
        for (int i = 0; i < bd->n_branches; i++) {
            free(bd->branches[i].control_values);
            free(bd->branches[i].state_values);
        }
        free(bd);
    }
}

int nept_bifdiagram_add_branch(BifurcationDiagram *bd,
                                const double *ctrl, const double *state,
                                int n_points, bool is_stable)
{
    if (!bd || bd->n_branches >= BIF_MAX_BRANCHES || n_points <= 0) return -1;
    BifurcationBranch *br = &bd->branches[bd->n_branches];
    br->control_values = (double *)malloc((size_t)n_points * sizeof(double));
    br->state_values = (double *)malloc((size_t)n_points * sizeof(double));
    if (!br->control_values || !br->state_values) {
        free(br->control_values);
        free(br->state_values);
        return -1;
    }
    memcpy(br->control_values, ctrl, (size_t)n_points * sizeof(double));
    memcpy(br->state_values, state, (size_t)n_points * sizeof(double));
    br->n_points = n_points;
    br->is_stable = is_stable;
    br->branch_id = bd->n_branches;
    bd->n_branches++;
    return br->branch_id;
}

/* --- Centre manifold reduction (simplified: return allocated structure for centre manifold reduction) --- */
CentreManifold *nept_centre_manifold_reduce(const double *jacobian, int dim,
                                              const double *nonlinear_terms,
                                              int n_nonlinear)
{
    if (!jacobian || dim <= 0) return NULL;
    CentreManifold *cm = (CentreManifold *)malloc(sizeof(CentreManifold));
    if (!cm) return NULL;
    cm->dim_centre = 1; /* assume codimension-1 bifurcation */
    cm->reduced_jacobian = (double *)malloc(sizeof(double));
    cm->normal_form_coeffs = (double *)calloc(4, sizeof(double));
    cm->centre_eigenvectors = NULL;
    if (!cm->reduced_jacobian || !cm->normal_form_coeffs) {
        free(cm->reduced_jacobian);
        free(cm->normal_form_coeffs);
        free(cm);
        return NULL;
    }
    /* Placeholder: set reduced Jacobian to dominant eigenvalue */
    cm->reduced_jacobian[0] = jacobian[0];
    (void)nonlinear_terms;
    (void)n_nonlinear;
    return cm;
}

void nept_centre_manifold_free(CentreManifold *cm)
{
    if (cm) {
        free(cm->reduced_jacobian);
        free(cm->normal_form_coeffs);
        free(cm->centre_eigenvectors);
        free(cm);
    }
}

/* --- Supercritical/subcritical check --- */
bool nept_is_supercritical(double cubic_coefficient)
{
    return (cubic_coefficient > 0.0);
}

/* ---------------------------------------------------------------------------
 * Lyapunov Exponent (Maximum)
 *
 * Knowledge point: The maximum Lyapunov exponent λ₁ characterizes the
 * average exponential rate of divergence of nearby trajectories:
 *
 *   |δx(t)| ≈ |δx(0)| * exp(λ₁*t)
 *
 * λ₁ > 0: chaos (sensitive dependence on initial conditions)
 * λ₁ = 0: marginal (periodic or quasiperiodic)
 * λ₁ < 0: stable fixed point
 *
 * Computation via Benettin's algorithm: evolve the system and a
 * tangent vector simultaneously, periodically renormalizing.
 *
 * Relationship to phase transitions: the Lyapunov time τ_L = 1/λ₁
 * diverges at a bifurcation point (critical slowing down).
 *
 * Reference: Benettin et al., Meccanica 15, 9 (1980)
 * --------------------------------------------------------------------------- */
double nept_max_lyapunov_exponent(
    void (*ode)(double t, const double *x, double *dxdt, void *params),
    const double *x0, int dim, void *params,
    double dt, int n_steps)
{
    if (!ode || !x0 || dim <= 0 || n_steps <= 0) return 0.0;

    /* Allocate and initialize state and tangent vectors */
    double *x = (double *)malloc((size_t)dim * sizeof(double));
    double *v = (double *)malloc((size_t)dim * sizeof(double));
    double *dx = (double *)malloc((size_t)dim * sizeof(double));
    double *dv = (double *)malloc((size_t)dim * sizeof(double));
    if (!x || !v || !dx || !dv) {
        free(x); free(v); free(dx); free(dv);
        return 0.0;
    }

    for (int i = 0; i < dim; i++) {
        x[i] = x0[i];
        v[i] = (i == 0) ? 1.0 : 0.0; /* initial unit tangent */
    }

    /* Normalize initial tangent vector */
    double norm = 0.0;
    for (int i = 0; i < dim; i++) norm += v[i] * v[i];
    norm = sqrt(norm);
    if (norm > 1e-15)
        for (int i = 0; i < dim; i++) v[i] /= norm;

    double sum_log = 0.0;
    int n_renorm = 0;

    for (int step = 0; step < n_steps; step++) {
        double t = (double)step * dt;

        /* Compute dx/dt = f(x) */
        ode(t, x, dx, params);

        /* Linearized dynamics: dv/dt = J(x)*v
         * Approximate J(x)*v ≈ [f(x + ε*v) - f(x - ε*v)] / (2ε) */
        double eps = 1e-6;
        double *x_plus = (double *)malloc((size_t)dim * sizeof(double));
        double *x_minus = (double *)malloc((size_t)dim * sizeof(double));
        double *f_plus = (double *)malloc((size_t)dim * sizeof(double));
        double *f_minus = (double *)malloc((size_t)dim * sizeof(double));

        if (x_plus && x_minus && f_plus && f_minus) {
            for (int i = 0; i < dim; i++) {
                x_plus[i] = x[i] + eps * v[i];
                x_minus[i] = x[i] - eps * v[i];
            }
            ode(t, x_plus, f_plus, params);
            ode(t, x_minus, f_minus, params);
            for (int i = 0; i < dim; i++) {
                dv[i] = (f_plus[i] - f_minus[i]) / (2.0 * eps);
            }
        }

        /* Euler integration */
        for (int i = 0; i < dim; i++) {
            x[i] += dx[i] * dt;
            v[i] += dv[i] * dt;
        }

        free(x_plus); free(x_minus); free(f_plus); free(f_minus);

        /* Renormalize every 10 steps */
        if ((step + 1) % 10 == 0 || step == n_steps - 1) {
            norm = 0.0;
            for (int i = 0; i < dim; i++) norm += v[i] * v[i];
            norm = sqrt(norm);
            if (norm > 1e-15) {
                sum_log += log(norm);
                for (int i = 0; i < dim; i++) v[i] /= norm;
                n_renorm++;
            }
        }
    }

    free(x); free(v); free(dx); free(dv);

    if (n_renorm == 0) return 0.0;
    double total_time = (double)n_steps * dt;
    return sum_log / total_time;
}

/* ---------------------------------------------------------------------------
 * Simple Numerical Continuation (1D)
 *
 * Knowledge point: Parameter continuation follows a branch of equilibria
 * as a control parameter varies. For 1D: f(x, μ) = 0.
 *
 * At each μ step, Newton's method is used to find x such that f(x,μ)=0:
 *   x_{n+1} = x_n - f(x_n, μ) / (∂f/∂x(x_n, μ))
 *
 * The initial guess for the next μ is the converged solution from the
 * previous step (natural continuation).
 *
 * Reference: Kuznetsov, §10.2; Seydel, Practical Bifurcation... (2010)
 * --------------------------------------------------------------------------- */
int nept_continue_branch(
    double (*f)(double x, double mu),
    double mu_start, double mu_end, int n_steps,
    double x_guess,
    double *ctrl_out, double *state_out)
{
    if (!f || !ctrl_out || !state_out || n_steps <= 0) return -1;

    double dmu = (mu_end - mu_start) / (double)n_steps;
    double x = x_guess;
    int count = 0;

    for (int i = 0; i <= n_steps; i++) {
        double mu = mu_start + (double)i * dmu;

        /* Newton iteration for f(x, mu) = 0 */
        for (int iter = 0; iter < 50; iter++) {
            double fx = f(x, mu);
            if (fabs(fx) < 1e-12) break;

            /* Numerical derivative: df/dx ≈ (f(x+h) - f(x-h))/(2h) */
            double h = 1e-6 * (fabs(x) + 1e-6);
            double dfdx = (f(x + h, mu) - f(x - h, mu)) / (2.0 * h);

            if (fabs(dfdx) < 1e-15) break; /* singular: bifurcation point */

            double dx = -fx / dfdx;
            x += dx;
            if (fabs(dx) < 1e-12) break;
        }

        ctrl_out[i] = mu;
        state_out[i] = x;
        count++;
    }

    return count;
}

/* --- Eigenvalue crossing check --- */
bool nept_eigenvalue_crossing(double lambda_real, double tolerance)
{
    return (fabs(lambda_real) < tolerance);
}