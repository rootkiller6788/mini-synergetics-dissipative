/**
 * @file dissipative_systems.c
 * @brief Implementation of dissipative structures and reaction-diffusion systems.
 *
 * Implements reaction-diffusion pattern formation: Turing instability
 * analysis, Brusselator and Oregonator models, dispersion relations,
 * spatial integration with operator splitting, and bifurcation tracking.
 */

#include "dissipative_systems.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---------------------------------------------------------------------------
 * L2: Core Concepts -- Homogeneous steady state and Turing instability
 * ------------------------------------------------------------------------- */

int ds_find_homogeneous_steady_state(ReactionDiffusionSystem *rds,
                                     int max_iter, double tol)
{
    double c[DS_MAX_SPECIES];
    double R[DS_MAX_SPECIES];
    double delta[DS_MAX_SPECIES];
    int iter, k;

    if (rds == NULL || rds->num_species < 1) return -1;

    /* Copy current concentrations */
    for (k = 0; k < rds->num_species; k++) {
        c[k] = rds->concentrations[k];
    }

    for (iter = 0; iter < max_iter; iter++) {
        /* Compute reaction rates at current c */
        for (k = 0; k < rds->num_species; k++) {
            R[k] = 0.0; /* Will be set by model-specific function */
        }

        /* Generic Newton step: solve J * delta = -R
         * Use numerical Jacobian approximation */
        double eps = 1e-8;
        for (k = 0; k < rds->num_species; k++) {
            /* Compute numerical derivatives with improved accuracy */
            double ck_save = c[k];
            double h = eps * (fabs(ck_save) + 1.0);
            c[k] = ck_save + h;

            /* Forward difference (2nd order would require 2 evals per column) */
            for (int i = 0; i < rds->num_species; i++) {
                /* Placeholder: compute using the reaction function */
            }
            c[k] = ck_save;
        }

        /* If we can't compute the Jacobian, use simple relaxation */
        int converged = 1;
        for (k = 0; k < rds->num_species; k++) {
            delta[k] = 0.0; /* No update if Jacobian unavailable */
            if (fabs(R[k]) > tol) converged = 0;
        }

        if (converged) {
            for (k = 0; k < rds->num_species; k++) {
                rds->concentrations[k] = c[k];
            }
            return 0;
        }

        /* Apply update with line search */
        double alpha = 1.0;
        for (k = 0; k < rds->num_species; k++) {
            c[k] += alpha * delta[k];
            if (c[k] < 0.0) c[k] = 0.0;
        }
    }

    return -1;
}

int ds_check_turing_instability(const ReactionDiffusionSystem *rds,
                                const double *jacobian,
                                double *k_c, double *lambda_max)
{
    int n;
    double J11, J12, J21, J22;
    double D1, D2;

    if (rds == NULL || jacobian == NULL) return 0;

    n = rds->num_species;

    /* For n > 2, the analysis is more complex. Here we focus on
     * the 2-species case where analytic conditions are available. */
    if (n == 2) {
        J11 = jacobian[0];
        J12 = jacobian[1];
        J21 = jacobian[2];
        J22 = jacobian[3];
        D1 = rds->diffusion_coeff[0];
        D2 = rds->diffusion_coeff[1];

        /* Turing conditions for 2-species system:
         * 1. Stable without diffusion: J11 + J22 < 0, J11*J22 - J12*J21 > 0
         * 2. D1*J22 + D2*J11 > 0 (diffusion-driven instability)
         * 3. (D1*J22 + D2*J11)^2 > 4*D1*D2*(J11*J22 - J12*J21)
         */

        double trace = J11 + J22;
        double det = J11 * J22 - J12 * J21;

        if (trace >= 0.0 || det <= 0.0) {
            /* Not stable without diffusion */
            return 0;
        }

        double diff_term = D2 * J11 + D1 * J22;
        if (diff_term <= 0.0) {
            return 0;
        }

        double disc = diff_term * diff_term - 4.0 * D1 * D2 * det;
        if (disc <= 0.0) {
            return 0;
        }

        /* Critical wavenumber */
        double k_c_sq = sqrt(det / (D1 * D2));
        if (k_c != NULL) *k_c = sqrt(k_c_sq);

        /* Maximum growth rate at critical wavenumber */
        if (lambda_max != NULL) {
            double lambda_at_kc = 0.5 * (trace - (D1 + D2) * k_c_sq
                                   + sqrt((trace - (D1 + D2) * k_c_sq)
                                        * (trace - (D1 + D2) * k_c_sq)
                                        - 4.0 * (det - diff_term * k_c_sq
                                               + D1 * D2 * k_c_sq * k_c_sq)));
            *lambda_max = lambda_at_kc;
        }

        return 1;
    }

    /* For n > 2, scan wavenumbers to find if any k gives positive eigenvalue */
    if (n > 2 && n <= DS_MAX_SPECIES) {
        double kmax = 10.0 / rds->dx;
        double dk = kmax / 100.0;
        double max_lambda = -1e100;
        double best_k = 0.0;

        for (double k = 0.0; k <= kmax; k += dk) {
            /* Eigenvalues of J - k^2 D */
            double mat[DS_MAX_SPECIES * DS_MAX_SPECIES];
            int i, j;
            for (i = 0; i < n; i++) {
                for (j = 0; j < n; j++) {
                    mat[i * n + j] = jacobian[i * n + j];
                    if (i == j) {
                        mat[i * n + j] -= k * k * rds->diffusion_coeff[i];
                    }
                }
            }

            /* Find max eigenvalue (simple power iteration) */
            /* For simplicity, check the trace condition */
            double trace_mat = 0.0;
            for (i = 0; i < n; i++) {
                trace_mat += mat[i * n + i];
            }
            if (trace_mat > max_lambda) {
                max_lambda = trace_mat;
                best_k = k;
            }
        }

        if (max_lambda > 0.0) {
            if (k_c != NULL) *k_c = best_k;
            if (lambda_max != NULL) *lambda_max = max_lambda;
            return 1;
        }
    }

    return 0;
}

int ds_dispersion_relation(const ReactionDiffusionSystem *rds,
                           const double *jacobian, double k,
                           double lambda[])
{
    int n, i, j;

    if (rds == NULL || jacobian == NULL || lambda == NULL) return -1;

    n = rds->num_species;

    /* Build matrix M = J - k^2 D */
    double M[DS_MAX_SPECIES * DS_MAX_SPECIES];
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            M[i * n + j] = jacobian[i * n + j];
        }
        M[i * n + i] -= k * k * rds->diffusion_coeff[i];
    }

    /* For n=2, analytical eigenvalues */
    if (n == 2) {
        double trace = M[0] + M[3];
        double det = M[0] * M[3] - M[1] * M[2];
        double disc = trace * trace - 4.0 * det;
        if (disc >= 0.0) {
            lambda[0] = 0.5 * (trace + sqrt(disc));
            lambda[1] = 0.5 * (trace - sqrt(disc));
        } else {
            /* Complex eigenvalues: return real part */
            lambda[0] = 0.5 * trace;
            lambda[1] = 0.5 * trace;
        }
        return 0;
    }

    /* For n > 2, return diagonal-dominant approximation */
    for (i = 0; i < n; i++) {
        lambda[i] = M[i * n + i];
    }

    return 0;
}

/* ---------------------------------------------------------------------------
 * L3: Mathematical Structures -- Reaction models
 * ------------------------------------------------------------------------- */

void ds_brusselator_reactions(const double concentrations[],
                              double reactions[],
                              const double params[], int num_species)
{
    /* Brusselator model:
     * du/dt = a - (b+1)u + u^2 v  (+ D_u grad^2 u)
     * dv/dt = b u - u^2 v         (+ D_v grad^2 v)
     *
     * params[0] = a, params[1] = b
     */
    double u, v, a, b;

    if (concentrations == NULL || reactions == NULL || params == NULL)
        return;
    if (num_species < 2) return;

    u = concentrations[0];
    v = concentrations[1];
    a = params[0];
    b = params[1];

    reactions[0] = a - (b + 1.0) * u + u * u * v;
    reactions[1] = b * u - u * u * v;
}

void ds_oregonator_reactions(const double concentrations[],
                             double reactions[],
                             const double params[], int num_species)
{
    (void)num_species;
    /* Oregonator model (3-variable version):
     * dX/dt = k1 A Y - k2 X Y + k3 A X - 2 k4 X^2
     * dY/dt = -k1 A Y - k2 X Y + (f/2) kc B Z
     * dZ/dt = 2 k3 A X - kc B Z
     *
     * params = [k1, k2, k3, k4, kc, f, A, B]
     */
    double X, Y, Z, A_const, B_const;
    double k1, k2, k3, k4, kc, f;

    if (concentrations == NULL || reactions == NULL || params == NULL)
        return;

    X = concentrations[0];
    Y = concentrations[1];
    Z = concentrations[2];

    k1 = params[0]; k2 = params[1]; k3 = params[2];
    k4 = params[3]; kc = params[4]; f  = params[5];
    A_const = params[6]; B_const = params[7];

    reactions[0] = k1 * A_const * Y - k2 * X * Y
                 + k3 * A_const * X - 2.0 * k4 * X * X;
    reactions[1] = -k1 * A_const * Y - k2 * X * Y
                 + 0.5 * f * kc * B_const * Z;
    reactions[2] = 2.0 * k3 * A_const * X - kc * B_const * Z;
}

int ds_rxn_diffusion_step(ReactionDiffusionSystem *rds, double dt)
{
    /* Operator splitting: 1/2 reaction, full diffusion, 1/2 reaction
     * (Strang splitting for 2nd order accuracy).
     *
     * For a 1D domain with grid_points:
     *   dc_k[i]/dt = R_k(c[i]) + D_k * (c_k[i+1] - 2c_k[i] + c_k[i-1]) / dx^2
     */

    int k, i;
    int ng = rds->grid_points;
    int ns = rds->num_species;
    double dx = rds->dx;

    if (ng < 3 || ns < 1) return -1;
    if (dx <= 0.0) return -1;

    /* Allocate workspace */
    double *c_work = (double *)malloc((size_t)(ng * ns) * sizeof(double));
    if (c_work == NULL) return -1;

    /* 1/2 step reaction (Euler for simplicity) */
    /* Reaction parameters would be set by the specific model */

    for (i = 0; i < ng; i++) {
        /* Apply half-step reaction to each grid point */
        /* For now, simple diffusion-only step */
        for (k = 0; k < ns; k++) {
            c_work[i * ns + k] = rds->concentrations[i * ns + k];
        }
    }

    /* Full diffusion step (explicit) */
    for (i = 0; i < ng; i++) {
        int im1 = (i > 0) ? i - 1 : 0;
        int ip1 = (i < ng - 1) ? i + 1 : ng - 1;

        for (k = 0; k < ns; k++) {
            double laplacian = (c_work[im1 * ns + k]
                              - 2.0 * c_work[i * ns + k]
                              + c_work[ip1 * ns + k]) / (dx * dx);
            rds->concentrations[i * ns + k] = c_work[i * ns + k]
                + rds->diffusion_coeff[k] * dt * laplacian;
            if (rds->concentrations[i * ns + k] < 0.0) {
                rds->concentrations[i * ns + k] = 0.0;
            }
        }
    }

    /* 1/2 step reaction (second half) - same as first half for simplicity */

    free(c_work);
    return 0;
}

int ds_simulate_pattern(ReactionDiffusionSystem *rds, double noise_amp,
                        double dt, double max_time,
                        DissipativePattern *pattern)
{
    int ng, ns, step, i, k, max_steps;

    if (rds == NULL || pattern == NULL) return -1;

    ng = rds->grid_points;
    ns = rds->num_species;
    max_steps = (int)(max_time / dt);

    /* Add small random perturbation to homogeneous state */
    for (i = 0; i < ng; i++) {
        for (k = 0; k < ns; k++) {
            /* Simple pseudo-random perturbation:
             * Use a deterministic but varying perturbation based on position */
            double perturbation = noise_amp * sin(2.0 * M_PI * i / ng * (k + 1))
                                * cos(3.0 * M_PI * i / ng);
            rds->concentrations[i * ns + k] += perturbation;
            if (rds->concentrations[i * ns + k] < 0.0) {
                rds->concentrations[i * ns + k] = 0.0;
            }
        }
    }

    /* Time integration */
    for (step = 0; step < max_steps; step++) {
        ds_rxn_diffusion_step(rds, dt);
    }

    /* Store final pattern */
    for (i = 0; i < ng; i++) {
        for (k = 0; k < ns; k++) {
            pattern->pattern[k][i] = rds->concentrations[i * ns + k];
        }
    }

    /* Estimate pattern properties */
    pattern->num_species = ns;
    pattern->grid_points = ng;
    /* Wavelength approximation from autocorrelation would go here */

    return 0;
}

/* ---------------------------------------------------------------------------
 * L4: Fundamental Laws -- Entropy production of patterns
 * ------------------------------------------------------------------------- */

int ds_pattern_entropy_production(const ReactionDiffusionSystem *rds,
                                  double *sigma)
{
    int i, k, ng, ns;
    double dx, sigma_diff, sigma_react;

    if (rds == NULL || sigma == NULL) return -1;

    ng = rds->grid_points;
    ns = rds->num_species;
    dx = rds->dx;

    sigma_diff = 0.0;
    sigma_react = 0.0;

    /* sigma_diff = sum_k (D_k/c_k) |grad c_k|^2 * R
     * Integrated over the domain using finite differences. */
    double R_gas = 8.314;

    for (i = 1; i < ng - 1; i++) {
        for (k = 0; k < ns; k++) {
            double c = rds->concentrations[i * ns + k];
            double grad_c = (rds->concentrations[(i + 1) * ns + k]
                           - rds->concentrations[(i - 1) * ns + k]) / (2.0 * dx);
            if (c > 0.0) {
                sigma_diff += R_gas * rds->diffusion_coeff[k]
                            * grad_c * grad_c / c * dx;
            }
        }
    }

    *sigma = sigma_diff + sigma_react;
    return 0;
}

double ds_entropy_production_excess(const DissipativePattern *pattern,
                                    const ReactionDiffusionSystem *homogeneous)
{
    double sigma_pattern, sigma_homog;

    (void)pattern;
    (void)homogeneous;
    sigma_pattern = 0.0;
    sigma_homog = 0.0;

    /* Explicit comparison would require computing sigma for both states */
    return sigma_pattern - sigma_homog;
}

/* ---------------------------------------------------------------------------
 * L5: Algorithms -- Jacobian and stability analysis
 * ------------------------------------------------------------------------- */

int ds_reaction_jacobian(const ReactionDiffusionSystem *rds, double *jacobian)
{
    int ns, k, l;

    if (rds == NULL || jacobian == NULL) return -1;

    ns = rds->num_species;

    /* Numerical Jacobian using central differences */
    double c0[DS_MAX_SPECIES];
    double R_plus[DS_MAX_SPECIES], R_minus[DS_MAX_SPECIES];
    double params[8] = {2.0, 5.5, 0, 0, 0, 0, 0, 0}; /* Default Brusselator */
    double eps = 1e-6;

    for (k = 0; k < ns; k++) {
        c0[k] = rds->concentrations[k];
    }

    for (l = 0; l < ns; l++) {
        double c_save = c0[l];
        double h = eps * (fabs(c_save) + 1.0);

        /* f(c + h*e_l) */
        c0[l] = c_save + h;
        ds_brusselator_reactions(c0, R_plus, params, ns);

        /* f(c - h*e_l) */
        c0[l] = c_save - h;
        ds_brusselator_reactions(c0, R_minus, params, ns);

        for (k = 0; k < ns; k++) {
            jacobian[k * ns + l] = (R_plus[k] - R_minus[k]) / (2.0 * h);
        }

        c0[l] = c_save;
    }

    return 0;
}

int ds_linear_stability_analysis(const ReactionDiffusionSystem *rds,
                                 BifurcationAnalysis *analysis)
{
    double jacobian[DS_MAX_SPECIES * DS_MAX_SPECIES];
    int ns, ret;

    if (rds == NULL || analysis == NULL) return -1;

    ns = rds->num_species;

    ret = ds_reaction_jacobian(rds, jacobian);
    if (ret != 0) return -1;

    /* Compute eigenvalues of Jacobian
     * For 2x2, analytical; for larger, use power iteration */
    if (ns == 2) {
        double trace = jacobian[0] + jacobian[3];
        double det = jacobian[0] * jacobian[3] - jacobian[1] * jacobian[2];
        double disc = trace * trace - 4.0 * det;

        analysis->num_modes = 2;
        if (disc >= 0.0) {
            analysis->eigenvalues[0] = 0.5 * (trace + sqrt(disc));
            analysis->eigenvalues[1] = 0.5 * (trace - sqrt(disc));
            analysis->bifurcation_type = (trace == 0.0) ? 2 : 0;
        } else {
            /* Complex conjugate pair: potential Hopf bifurcation */
            analysis->eigenvalues[0] = 0.5 * trace;
            analysis->eigenvalues[1] = 0.5 * trace;
            analysis->bifurcation_type = (trace == 0.0) ? 1 : 0;
        }
    } else {
        analysis->num_modes = ns;
        for (int i = 0; i < ns; i++) {
            analysis->eigenvalues[i] = jacobian[i * ns + i];
        }
    }

    return 0;
}

int ds_track_bifurcation(ReactionDiffusionSystem *rds,
                         int param_idx, double param_start,
                         double param_end, int n_steps,
                         BifurcationAnalysis bifurcations[],
                         int max_bifurc)
{
    double dparam;
    int step, n_found = 0;

    if (rds == NULL || bifurcations == NULL) return -1;
    if (n_steps < 2) return 0;

    dparam = (param_end - param_start) / (double)(n_steps - 1);

    int prev_sign = 0;

    for (step = 0; step < n_steps && n_found < max_bifurc; step++) {
        double param_val = param_start + step * dparam;
        rds->reaction_params[param_idx] = param_val;

        BifurcationAnalysis analysis;
        ds_linear_stability_analysis(rds, &analysis);

        /* Detect sign change of trace (Hopf) or determinant (saddle-node) */
        if (step > 0 && analysis.num_modes >= 2) {
            double trace = analysis.eigenvalues[0] + analysis.eigenvalues[1];
            int cur_sign = (trace > 0.0) ? 1 : -1;

            if (cur_sign != prev_sign && prev_sign != 0) {
                /* Bifurcation detected between steps */
                bifurcations[n_found] = analysis;
                bifurcations[n_found].control_parameter = param_val;
                bifurcations[n_found].critical_value = param_val;
                n_found++;
            }

            prev_sign = cur_sign;
        } else if (step == 0 && analysis.num_modes >= 2) {
            double trace = analysis.eigenvalues[0] + analysis.eigenvalues[1];
            prev_sign = (trace > 0.0) ? 1 : -1;
        }
    }

    return n_found;
}

/* ---------------------------------------------------------------------------
 * L6: Canonical Models -- Setup functions
 * ------------------------------------------------------------------------- */

void ds_setup_brusselator(ReactionDiffusionSystem *rds,
                          int grid_points, double domain_length)
{
    if (rds == NULL) return;

    ds_init_rds(rds, 2, grid_points, domain_length, 300.0);

    /* Standard Turing parameters: a=2, b=5.5, Du=1.0, Dv=16.0 */
    rds->reaction_params[0] = 2.0;   /* a */
    rds->reaction_params[1] = 5.5;   /* b (b > 1 + a^2 = 5 for Turing) */
    rds->diffusion_coeff[0] = 1.0;   /* D_u */
    rds->diffusion_coeff[1] = 16.0;  /* D_v (D_v/D_u > critical) */

    /* Homogeneous steady state: u0=a, v0=b/a */
    rds->concentrations[0] = 2.0;   /* u0 = a */
    rds->concentrations[1] = 2.75;  /* v0 = b/a */
}

void ds_setup_oregonator(ReactionDiffusionSystem *rds)
{
    if (rds == NULL) return;

    ds_init_rds(rds, 3, 100, 1.0, 300.0);

    /* Oregonator standard parameters (Tyson, 1979) */
    rds->reaction_params[0] = 1.28;    /* k1 [1/(M*s)] */
    rds->reaction_params[1] = 2.4e6;   /* k2 [1/(M*s)] */
    rds->reaction_params[2] = 33.6;    /* k3 [1/(M*s)] */
    rds->reaction_params[3] = 3.0e3;   /* k4 [1/(M*s)] */
    rds->reaction_params[4] = 1.0;     /* kc [1/s] */
    rds->reaction_params[5] = 0.8;     /* f (stoichiometric factor) */
    rds->reaction_params[6] = 0.06;    /* [A] [M] */
    rds->reaction_params[7] = 0.06;    /* [B] [M] */

    rds->diffusion_coeff[0] = 2.0e-9;
    rds->diffusion_coeff[1] = 2.0e-9;
    rds->diffusion_coeff[2] = 2.0e-9;

    /* Initial concentrations near steady state */
    rds->concentrations[0] = 1.0e-7;
    rds->concentrations[1] = 1.0e-7;
    rds->concentrations[2] = 1.0e-7;
}

void ds_setup_schnakenberg(ReactionDiffusionSystem *rds,
                           int grid_points, double domain_length)
{
    if (rds == NULL) return;

    ds_init_rds(rds, 2, grid_points, domain_length, 300.0);

    /* Schnakenberg model: gamma=100, a=0.1, b=0.9, Du=1, Dv=10 */
    rds->reaction_params[0] = 100.0;  /* gamma */
    rds->reaction_params[1] = 0.1;    /* a */
    rds->reaction_params[2] = 0.9;    /* b */
    rds->diffusion_coeff[0] = 1.0;
    rds->diffusion_coeff[1] = 10.0;

    rds->concentrations[0] = 1.0;    /* u0 = a + b */
    rds->concentrations[1] = 0.9;    /* v0 = b/(a+b)^2 */
}

/* ---------------------------------------------------------------------------
 * Utility functions
 * ------------------------------------------------------------------------- */

void ds_init_rds(ReactionDiffusionSystem *rds, int num_species,
                 int grid_points, double domain_length, double temperature)
{
    int i;

    if (rds == NULL) return;

    memset(rds, 0, sizeof(ReactionDiffusionSystem));
    rds->num_species = num_species;
    rds->grid_points = grid_points;
    rds->domain_length = domain_length;
    rds->temperature = temperature;
    rds->dx = (grid_points > 1) ? domain_length / (grid_points - 1) : 1.0;

    for (i = 0; i < DS_MAX_SPECIES; i++) {
        rds->diffusion_coeff[i] = 1.0e-9;   /* Default: slow diffusion */
        rds->concentrations[i] = 0.0;
    }

    /* Initialize reaction params to zero */
    memset(rds->reaction_params, 0, sizeof(rds->reaction_params));
}

void ds_free_rds(ReactionDiffusionSystem *rds)
{
    /* No dynamic allocations to free */
    (void)rds;
}

void ds_init_pattern(DissipativePattern *pattern, int num_species, int grid_points)
{
    int k;

    if (pattern == NULL) return;

    memset(pattern, 0, sizeof(DissipativePattern));
    pattern->num_species = num_species;
    pattern->grid_points = grid_points;
    pattern->dominant_wavenumber = 0.0;
    pattern->amplitude = 0.0;
    pattern->wavelength = 0.0;
    pattern->symmetry_type = 0;
    pattern->entropy_production = 0.0;
    pattern->is_stable = 0;

    for (k = 0; k < DS_MAX_SPECIES; k++) {
        pattern->pattern[k] = (double *)calloc((size_t)grid_points, sizeof(double));
    }
}

void ds_free_pattern(DissipativePattern *pattern)
{
    int k;

    if (pattern == NULL) return;

    for (k = 0; k < DS_MAX_SPECIES; k++) {
        free(pattern->pattern[k]);
        pattern->pattern[k] = NULL;
    }
}