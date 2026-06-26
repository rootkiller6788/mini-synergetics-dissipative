/**
 * turing_models.c — Reaction Kinetics Models for Turing Pattern Formation
 *
 * Implements 7 classical reaction-diffusion kinetics models,
 * each with: reaction rates, Jacobian, HSS finder, and parameter validation.
 *
 * Reference: Turing, A.M. (1952) "The Chemical Basis of Morphogenesis"
 *            Murray, J.D. (2003) Mathematical Biology II: Spatial Models
 *            Pearson, J.E. (1993) Science 261:189-192
 *            Gierer & Meinhardt (1972) Kybernetik 12(1):30-39
 *            Prigogine & Lefever (1968) J. Chem. Phys. 48(4):1695-1700
 *
 * Knowledge Coverage:
 *   L1: Model definitions (Gray-Scott, Gierer-Meinhardt, Brusselator, etc.)
 *   L2: Activator-inhibitor concept, short-range activation, long-range inhibition
 *   L3: Reaction rate functions, Jacobian matrices, nonlinear kinetics
 *   L4: Turing condition verification for each model
 */

#include "turing_models.h"
#include "turing_analysis.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

/* ==========================================================================
 * Model 0: Gray-Scott
 *
 * f(u,v) = F*(1-u) - u*v^2
 * g(u,v) = u*v^2 - (F+k)*v
 *
 * HSS: Solve F*(1-u) = u*v^2, u*v^2 = (F+k)*v
 *      If v ≠ 0: u = (F+k)*v / v^2 = (F+k)/v → no, wrong
 *      Properly: u*v^2 = (F+k)*v → v*(u*v - (F+k)) = 0
 *      So either v=0 → u=1 (trivial HSS)
 *      Or v = (u*v^2)/(F+k) ... better: F*(1-u) = (F+k)*v
 *      From g=0: v*(u*v - (F+k)) = 0
 *      If v ≠ 0: u*v = F+k → v = (F+k)/u
 *      From f=0: F*(1-u) = u*v^2 = u*((F+k)/u)² = (F+k)²/u
 *      → F*u*(1-u) = (F+k)²
 *      → F*u - F*u² = (F+k)²
 *      → F*u² - F*u + (F+k)² = 0
 *      → u = [F ± sqrt(F² - 4F*(F+k)²)] / (2F)
 *      → u = [1 ± sqrt(1 - 4*(F+k)²/F)] / 2
 *
 *      Requires: 4*(F+k)²/F ≤ 1 for real roots
 * ==========================================================================*/

int grayscott_reaction(double u, double v, const ModelParams *p,
                       double *f_out, double *g_out) {
    if (!p || !f_out || !g_out) return -1;
    if (!is_finite_val(u) || !is_finite_val(v)) return -1;

    double u_use = (u < 0.0) ? 0.0 : u;
    double v_use = (v < 0.0) ? 0.0 : v;

    *f_out = p->F * (1.0 - u_use) - u_use * v_use * v_use;
    *g_out = u_use * v_use * v_use - (p->F + p->k) * v_use;

    return 0;
}

int grayscott_jacobian(double u, double v, const ModelParams *p,
                       double J[2][2]) {
    if (!p || !J) return -1;

    double u_use = (u < 0.0) ? 0.0 : u;
    double v_use = (v < 0.0) ? 0.0 : v;

    J[0][0] = -p->F - v_use * v_use;     /* ∂f/∂u */
    J[0][1] = -2.0 * u_use * v_use;       /* ∂f/∂v */
    J[1][0] = v_use * v_use;              /* ∂g/∂u */
    J[1][1] = 2.0 * u_use * v_use - (p->F + p->k); /* ∂g/∂v */

    return 0;
}

int grayscott_find_hss(const ModelParams *p, HomogeneousSteadyState *hss) {
    if (!p || !hss) return -1;

    double F = p->F, kap = p->k;
    double Fpk2 = (F + kap) * (F + kap);

    /* Discriminant: 1 - 4*(F+k)²/F */
    double disc = 1.0 - 4.0 * Fpk2 / F;

    if (disc < 0.0) {
        /* No real non-trivial HSS — return trivial */
        hss->u_star = 1.0;
        hss->v_star = 0.0;
        hss->converged = 1;
        hss->iterations = 0;
        hss->residual = 0.0;
        return 0;
    }

    double sqrt_disc = sqrt(disc);
    double u1 = (1.0 + sqrt_disc) / 2.0;
    double u2 = (1.0 - sqrt_disc) / 2.0;

    /* Choose the larger u (more physically relevant) */
    double u_star = (u1 > u2) ? u1 : u2;
    double v_star = (F + kap) / u_star;

    hss->u_star = u_star;
    hss->v_star = v_star;
    hss->converged = 1;
    hss->iterations = 1;

    /* Verify */
    double f_val, g_val;
    grayscott_reaction(u_star, v_star, p, &f_val, &g_val);
    hss->residual = sqrt(f_val*f_val + g_val*g_val);

    return 0;
}

int grayscott_validate_params(const ModelParams *p) {
    return (p && p->F > 0.0 && p->k > 0.0);
}

/* ==========================================================================
 * Model 1: Gierer-Meinhardt
 *
 * Standard simplified form:
 *   f(u,v) = a - mu_a*u + rho_a * u² / ((1 + kappa_u*u²) * v)
 *   g(u,v) = b - mu_h*v + rho_h * u²
 *
 * For simplicity, use the common form without saturation (kappa=0):
 *   f(u,v) = a - mu_a*u + rho_a * u² / v
 *   g(u,v) = b - mu_h*v + rho_h * u²
 *
 * HSS: a - mu_a*u + rho_a*u²/v = 0, b - mu_h*v + rho_h*u² = 0
 *      From g: v = (b + rho_h*u²) / mu_h
 *      Substitute into f: a - mu_a*u + rho_a*u² * mu_h/(b + rho_h*u²) = 0
 *      Solve via Newton iteration.
 * ==========================================================================*/

int giere_meinhardt_reaction(double u, double v, const ModelParams *p,
                              double *f_out, double *g_out) {
    if (!p || !f_out || !g_out) return -1;
    if (!is_finite_val(u) || !is_finite_val(v)) return -1;

    double u_use = (u < 0.0) ? 0.0 : u;
    double v_use = (v < 1e-12) ? 1e-12 : v;  /* avoid division by zero */

    double u2 = u_use * u_use;
    double denom = v_use * (1.0 + p->rho_a * u2);  /* saturation in denominator */

    *f_out = p->c + p->rho_a * u2 / v_use - p->mu_a * u_use;
    *g_out = p->c + p->rho_h * u2 - p->mu_h * v_use;

    return 0;
}

int giere_meinhardt_jacobian(double u, double v, const ModelParams *p,
                              double J[2][2]) {
    if (!p || !J) return -1;

    double u_use = (u < 0.0) ? 0.0 : u;
    double v_use = (v < 1e-12) ? 1e-12 : v;

    J[0][0] = 2.0 * p->rho_a * u_use / v_use - p->mu_a;   /* ∂f/∂u */
    J[0][1] = -p->rho_a * u_use * u_use / (v_use * v_use); /* ∂f/∂v */
    J[1][0] = 2.0 * p->rho_h * u_use;                       /* ∂g/∂u */
    J[1][1] = -p->mu_h;                                     /* ∂g/∂v */

    return 0;
}

int giere_meinhardt_find_hss(const ModelParams *p, HomogeneousSteadyState *hss) {
    if (!p || !hss) return -1;

    /* Use Newton's method: solve f(u,v)=0, g(u,v)=0
     * Initial guess from simplified analysis */
    double u = 1.0, v = 1.0;
    int max_iter = 50;
    double tol = 1e-10;

    for (int iter = 0; iter < max_iter; iter++) {
        double f, g, J[2][2];
        giere_meinhardt_reaction(u, v, p, &f, &g);
        giere_meinhardt_jacobian(u, v, p, J);

        double det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
        if (fabs(det) < 1e-15) {
            hss->converged = 0;
            hss->iterations = iter;
            return -1;
        }

        /* Newton step: [Δu,Δv]ᵀ = -J⁻¹[f,g]ᵀ */
        double du = -( J[1][1]*f - J[0][1]*g) / det;
        double dv = -(-J[1][0]*f + J[0][0]*g) / det;

        u += du;
        v += dv;

        /* Enforce positivity */
        if (u < 0.0) u = 1e-6;
        if (v < 0.0) v = 1e-6;

        if (fabs(du) < tol && fabs(dv) < tol) {
            hss->u_star = u;
            hss->v_star = v;
            hss->converged = 1;
            hss->iterations = iter + 1;
            hss->residual = sqrt(f*f + g*g);
            return 0;
        }
    }

    hss->converged = 0;
    hss->iterations = max_iter;
    return -1;
}

int giere_meinhardt_validate_params(const ModelParams *p) {
    return (p && p->rho_a > 0.0 && p->rho_h > 0.0 && p->mu_a > 0.0 && p->mu_h > 0.0);
}

/* ==========================================================================
 * Model 2: Brusselator
 *
 * f(u,v) = A - (B+1)*u + u²*v
 * g(u,v) = B*u - u²*v
 *
 * HSS: u* = A, v* = B/A
 *      f(A, B/A) = A - (B+1)*A + A²*(B/A) = A - (B+1)A + A*B = A - A*B - A + A*B = 0 ✓
 *      g(A, B/A) = B*A - A²*(B/A) = A*B - A*B = 0 ✓
 * ==========================================================================*/

int brusselator_reaction(double u, double v, const ModelParams *p,
                         double *f_out, double *g_out) {
    if (!p || !f_out || !g_out) return -1;

    double u_use = (u < 0.0) ? 0.0 : u;
    double v_use = (v < 0.0) ? 0.0 : v;

    *f_out = p->A - (p->B + 1.0) * u_use + u_use * u_use * v_use;
    *g_out = p->B * u_use - u_use * u_use * v_use;

    return 0;
}

int brusselator_jacobian(double u, double v, const ModelParams *p,
                         double J[2][2]) {
    if (!p || !J) return -1;

    J[0][0] = -(p->B + 1.0) + 2.0 * u * v;  /* ∂f/∂u */
    J[0][1] = u * u;                          /* ∂f/∂v */
    J[1][0] = p->B - 2.0 * u * v;            /* ∂g/∂u */
    J[1][1] = -u * u;                         /* ∂g/∂v */

    return 0;
}

int brusselator_find_hss(const ModelParams *p, HomogeneousSteadyState *hss) {
    if (!p || !hss) return -1;

    if (p->A <= 0.0) return -1;

    hss->u_star = p->A;
    hss->v_star = p->B / p->A;
    hss->converged = 1;
    hss->iterations = 0;
    hss->residual = 0.0;

    return 0;
}

int brusselator_validate_params(const ModelParams *p) {
    return (p && p->A > 0.0 && p->B > 0.0);
}

/* ==========================================================================
 * Model 3: Schnakenberg
 *
 * f(u,v) = gamma * (a - u + u²*v)
 * g(u,v) = gamma * (b - u²*v)
 *
 * In many formulations gamma=1 (absorbed into diffusion ratio).
 *
 * HSS: u* = a + b, v* = b / ((a+b)²)
 *      f(a+b, v*) = γ*(a - (a+b) + (a+b)²*v*) = γ*(-b + (a+b)²*b/(a+b)²) = 0 ✓
 *      g(a+b, v*) = γ*(b - (a+b)²*v*) = γ*(b - (a+b)²*b/(a+b)²) = 0 ✓
 * ==========================================================================*/

int schnakenberg_reaction(double u, double v, const ModelParams *p,
                          double *f_out, double *g_out) {
    if (!p || !f_out || !g_out) return -1;

    double u_use = (u < 0.0) ? 0.0 : u;
    double u2 = u_use * u_use;

    *f_out = p->a - u_use + u2 * v;
    *g_out = p->b - u2 * v;

    return 0;
}

int schnakenberg_jacobian(double u, double v, const ModelParams *p,
                          double J[2][2]) {
    if (!p || !J) return -1;

    (void)p; /* a,b not needed in Jacobian */

    J[0][0] = -1.0 + 2.0 * u * v;   /* ∂f/∂u */
    J[0][1] = u * u;                  /* ∂f/∂v */
    J[1][0] = -2.0 * u * v;          /* ∂g/∂u */
    J[1][1] = -u * u;                /* ∂g/∂v */

    return 0;
}

int schnakenberg_find_hss(const ModelParams *p, HomogeneousSteadyState *hss) {
    if (!p || !hss) return -1;

    double u_star = p->a + p->b;
    if (u_star < 1e-15) return -1;

    hss->u_star = u_star;
    hss->v_star = p->b / (u_star * u_star);
    hss->converged = 1;
    hss->iterations = 0;
    hss->residual = 0.0;

    return 0;
}

int schnakenberg_validate_params(const ModelParams *p) {
    return (p && p->a >= 0.0 && p->b >= 0.0 && (p->a + p->b) > 0.0);
}

/* ==========================================================================
 * Model 4: FitzHugh-Nagumo (excitable media variant)
 *
 * f(u,v) = u - u³/3 - v + I_ext
 * g(u,v) = ε*(u + β - γ*v)
 *
 * The cubic nonlinearity allows bistability and excitability.
 *
 * HSS: From g=0: v = (u + β)/γ
 *      Substitute into f: u - u³/3 - (u+β)/γ + I_ext = 0
 *      Solve cubic via Newton.
 * ==========================================================================*/

int fitzhugh_nagumo_reaction(double u, double v, const ModelParams *p,
                             double *f_out, double *g_out) {
    if (!p || !f_out || !g_out) return -1;

    double u3 = u * u * u;
    *f_out = u - u3 / 3.0 - v + p->I_ext;
    *g_out = p->epsilon * (u + p->beta - p->gamma * v);

    return 0;
}

int fitzhugh_nagumo_jacobian(double u, double v, const ModelParams *p,
                             double J[2][2]) {
    if (!p || !J) return -1;

    (void)v;
    J[0][0] = 1.0 - u * u;                /* ∂f/∂u = 1 - u² */
    J[0][1] = -1.0;                        /* ∂f/∂v */
    J[1][0] = p->epsilon;                  /* ∂g/∂u */
    J[1][1] = -p->epsilon * p->gamma;      /* ∂g/∂v */

    return 0;
}

int fitzhugh_nagumo_find_hss(const ModelParams *p, HomogeneousSteadyState *hss) {
    if (!p || !hss) return -1;

    /* Newton solve for u (v eliminated): u - u³/3 - (u+β)/γ + I_ext = 0 */
    double u = 0.0;  /* initial guess */
    int max_iter = 100;
    double tol = 1e-12;

    for (int iter = 0; iter < max_iter; iter++) {
        double f = u - u*u*u/3.0 - (u + p->beta)/p->gamma + p->I_ext;
        double df = 1.0 - u*u - 1.0/p->gamma;

        if (fabs(df) < 1e-15) {
            hss->converged = 0;
            hss->iterations = iter;
            return -1;
        }

        double du = -f / df;
        u += du;

        if (fabs(du) < tol) {
            hss->u_star = u;
            hss->v_star = (u + p->beta) / p->gamma;
            hss->converged = 1;
            hss->iterations = iter + 1;
            hss->residual = fabs(f);
            return 0;
        }
    }

    hss->converged = 0;
    hss->iterations = max_iter;
    return -1;
}

int fitzhugh_nagumo_validate_params(const ModelParams *p) {
    return (p && p->epsilon > 0.0 && p->gamma > 0.0);
}

/* ==========================================================================
 * Model 5: Lengyel-Epstein (CIMA reaction)
 *
 * f(u,v) = a - u - 4uv/(1+u²)
 * g(u,v) = σ * [b * (u - uv/(1+u²))]
 *
 * where σ = c (starch complexation factor).
 *
 * HSS: From g=0: u - uv/(1+u²) = 0 → v = 1 + u² (if u ≠ 0)
 *      From f=0: a - u - 4u(1+u²)/(1+u²) = a - u - 4u = a - 5u → u* = a/5
 *      Then: v* = 1 + (a/5)²
 * ==========================================================================*/

int lengyel_epstein_reaction(double u, double v, const ModelParams *p,
                             double *f_out, double *g_out) {
    if (!p || !f_out || !g_out) return -1;

    double u_use = (u < 0.0) ? 0.0 : u;
    double denom = 1.0 + u_use * u_use;
    double uv_over_denom = u_use * v / denom;

    *f_out = p->a - u_use - 4.0 * uv_over_denom;
    *g_out = p->c * p->b * (u_use - uv_over_denom);  /* σ = c, factor b included */

    return 0;
}

int lengyel_epstein_jacobian(double u, double v, const ModelParams *p,
                             double J[2][2]) {
    if (!p || !J) return -1;

    double u2 = u * u;
    double denom = 1.0 + u2;
    double denom2 = denom * denom;

    /* ∂/∂u [uv/(1+u²)] = v * (1+u² - 2u²) / (1+u²)² = v*(1-u²)/(1+u²)² */
    double d_uv_d_u = v * (1.0 - u2) / denom2;
    double d_uv_d_v = u / denom;

    J[0][0] = -1.0 - 4.0 * d_uv_d_u;   /* ∂f/∂u */
    J[0][1] = -4.0 * d_uv_d_v;          /* ∂f/∂v */
    J[1][0] = p->c * p->b * (1.0 - d_uv_d_u); /* ∂g/∂u */
    J[1][1] = p->c * p->b * (-d_uv_d_v);       /* ∂g/∂v */

    return 0;
}

int lengyel_epstein_find_hss(const ModelParams *p, HomogeneousSteadyState *hss) {
    if (!p || !hss) return -1;

    double u_star = p->a / 5.0;
    double v_star = 1.0 + u_star * u_star;

    hss->u_star = u_star;
    hss->v_star = v_star;
    hss->converged = 1;
    hss->iterations = 0;
    hss->residual = 0.0;

    return 0;
}

int lengyel_epstein_validate_params(const ModelParams *p) {
    return (p && p->a > 0.0 && p->b > 0.0);
}

/* ==========================================================================
 * Model 6: Thomas (Substrate-Inhibition Kinetics)
 *
 * f(u,v) = a - u - ρ * R(u,v)
 * g(u,v) = α*(b - v) - ρ * R(u,v)
 *
 * where R(u,v) = uv / (1 + u + K*u²)
 *
 * This is a realistic model for the uricase-urate-O₂ reaction.
 * ==========================================================================*/

static double thomas_R(double u, double v, double K) {
    double denom = 1.0 + u + K * u * u;
    if (fabs(denom) < 1e-15) return 0.0;
    return u * v / denom;
}

int thomas_reaction(double u, double v, const ModelParams *p,
                    double *f_out, double *g_out) {
    if (!p || !f_out || !g_out) return -1;

    double u_use = (u < 0.0) ? 0.0 : u;
    double v_use = (v < 0.0) ? 0.0 : v;
    double R = thomas_R(u_use, v_use, p->k);  /* K = p->k for Thomas */

    *f_out = p->a - u_use - p->rho * R;
    *g_out = p->alpha * (p->b - v_use) - p->rho * R;

    return 0;
}

int thomas_jacobian(double u, double v, const ModelParams *p,
                    double J[2][2]) {
    if (!p || !J) return -1;

    double u_use = (u < 0.0) ? 0.0 : u;
    double denom = 1.0 + u_use + p->k * u_use * u_use;
    if (fabs(denom) < 1e-15) {
        J[0][0] = J[0][1] = J[1][0] = J[1][1] = 0.0;
        return 0;
    }

    double denom2 = denom * denom;
    /* R = uv/D, D=1+u+Ku²
     * ∂R/∂u = v*(D - u*D')/D² = v*(1+u+Ku² - u*(1+2Ku))/D² = v*(1-Ku²)/D² */
    double dR_du = v * (1.0 - p->k * u_use * u_use) / denom2;
    double dR_dv = u_use / denom;

    J[0][0] = -1.0 - p->rho * dR_du;                /* ∂f/∂u */
    J[0][1] = -p->rho * dR_dv;                       /* ∂f/∂v */
    J[1][0] = -p->rho * dR_du;                       /* ∂g/∂u */
    J[1][1] = -p->alpha - p->rho * dR_dv;            /* ∂g/∂v */

    return 0;
}

int thomas_find_hss(const ModelParams *p, HomogeneousSteadyState *hss) {
    if (!p || !hss) return -1;

    /* Newton solve for (u,v) */
    double u = p->a / 2.0;  /* initial guess */
    double v = p->b / 2.0;
    int max_iter = 100;
    double tol = 1e-12;

    for (int iter = 0; iter < max_iter; iter++) {
        double f, g, J[2][2];
        thomas_reaction(u, v, p, &f, &g);
        thomas_jacobian(u, v, p, J);

        double det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
        if (fabs(det) < 1e-15) { u += 0.1; v += 0.1; continue; }

        double du = -( J[1][1]*f - J[0][1]*g) / det;
        double dv = -(-J[1][0]*f + J[0][0]*g) / det;

        u += du; v += dv;
        if (u < 0.0) u = 1e-6;
        if (v < 0.0) v = 1e-6;

        if (fabs(du) < tol && fabs(dv) < tol) {
            hss->u_star = u;
            hss->v_star = v;
            hss->converged = 1;
            hss->iterations = iter + 1;
            hss->residual = sqrt(f*f + g*g);
            return 0;
        }
    }

    hss->converged = 0;
    hss->iterations = max_iter;
    return -1;
}

int thomas_validate_params(const ModelParams *p) {
    return (p && p->a > 0.0 && p->rho > 0.0 && p->alpha >= 0.0);
}

/* ==========================================================================
 * Universal Dispatchers
 * ==========================================================================*/

int model_reaction(ModelType type, double u, double v, const ModelParams *p,
                   double *f_out, double *g_out) {
    switch (type) {
        case MODEL_GRAY_SCOTT:       return grayscott_reaction(u, v, p, f_out, g_out);
        case MODEL_GIERER_MEINHARDT: return giere_meinhardt_reaction(u, v, p, f_out, g_out);
        case MODEL_BRUSSELATOR:      return brusselator_reaction(u, v, p, f_out, g_out);
        case MODEL_SCHNAKENBERG:     return schnakenberg_reaction(u, v, p, f_out, g_out);
        case MODEL_FITZHUGH_NAGUMO:  return fitzhugh_nagumo_reaction(u, v, p, f_out, g_out);
        case MODEL_LENGYEL_EPSTEIN:  return lengyel_epstein_reaction(u, v, p, f_out, g_out);
        case MODEL_THOMAS:           return thomas_reaction(u, v, p, f_out, g_out);
        default: return -1;
    }
}

int model_jacobian(ModelType type, double u, double v, const ModelParams *p,
                   double J[2][2]) {
    switch (type) {
        case MODEL_GRAY_SCOTT:       return grayscott_jacobian(u, v, p, J);
        case MODEL_GIERER_MEINHARDT: return giere_meinhardt_jacobian(u, v, p, J);
        case MODEL_BRUSSELATOR:      return brusselator_jacobian(u, v, p, J);
        case MODEL_SCHNAKENBERG:     return schnakenberg_jacobian(u, v, p, J);
        case MODEL_FITZHUGH_NAGUMO:  return fitzhugh_nagumo_jacobian(u, v, p, J);
        case MODEL_LENGYEL_EPSTEIN:  return lengyel_epstein_jacobian(u, v, p, J);
        case MODEL_THOMAS:           return thomas_jacobian(u, v, p, J);
        default: return -1;
    }
}

int model_find_hss(ModelType type, const ModelParams *p,
                   HomogeneousSteadyState *hss) {
    switch (type) {
        case MODEL_GRAY_SCOTT:       return grayscott_find_hss(p, hss);
        case MODEL_GIERER_MEINHARDT: return giere_meinhardt_find_hss(p, hss);
        case MODEL_BRUSSELATOR:      return brusselator_find_hss(p, hss);
        case MODEL_SCHNAKENBERG:     return schnakenberg_find_hss(p, hss);
        case MODEL_FITZHUGH_NAGUMO:  return fitzhugh_nagumo_find_hss(p, hss);
        case MODEL_LENGYEL_EPSTEIN:  return lengyel_epstein_find_hss(p, hss);
        case MODEL_THOMAS:           return thomas_find_hss(p, hss);
        default: return -1;
    }
}

int model_validate_params(ModelType type, const ModelParams *p) {
    switch (type) {
        case MODEL_GRAY_SCOTT:       return grayscott_validate_params(p);
        case MODEL_GIERER_MEINHARDT: return giere_meinhardt_validate_params(p);
        case MODEL_BRUSSELATOR:      return brusselator_validate_params(p);
        case MODEL_SCHNAKENBERG:     return schnakenberg_validate_params(p);
        case MODEL_FITZHUGH_NAGUMO:  return fitzhugh_nagumo_validate_params(p);
        case MODEL_LENGYEL_EPSTEIN:  return lengyel_epstein_validate_params(p);
        case MODEL_THOMAS:           return thomas_validate_params(p);
        default: return 0;
    }
}

const char *model_name(ModelType type) {
    switch (type) {
        case MODEL_GRAY_SCOTT:       return "Gray-Scott";
        case MODEL_GIERER_MEINHARDT: return "Gierer-Meinhardt";
        case MODEL_BRUSSELATOR:      return "Brusselator";
        case MODEL_SCHNAKENBERG:     return "Schnakenberg";
        case MODEL_FITZHUGH_NAGUMO:  return "FitzHugh-Nagumo";
        case MODEL_LENGYEL_EPSTEIN:  return "Lengyel-Epstein";
        case MODEL_THOMAS:           return "Thomas";
        case MODEL_USER_DEFINED:     return "User-Defined";
        default:                     return "Unknown";
    }
}

const char *model_description(ModelType type) {
    switch (type) {
        case MODEL_GRAY_SCOTT:
            return "Isothermal autocatalytic system. Produces spots, stripes, "
                   "labyrinths, and spatiotemporal chaos. (Gray & Scott 1984)";
        case MODEL_GIERER_MEINHARDT:
            return "Classic activator-inhibitor model for biological morphogenesis. "
                   "Short-range activation, long-range inhibition. (Gierer & Meinhardt 1972)";
        case MODEL_BRUSSELATOR:
            return "Trimolecular reaction-diffusion model from Prigogine school. "
                   "Prototype for nonequilibrium self-organization. (Prigogine & Lefever 1968)";
        case MODEL_SCHNAKENBERG:
            return "Simplest realistic trimolecular reaction with clean time-scale "
                   "separation. Popular in mathematical biology. (Schnakenberg 1979)";
        case MODEL_FITZHUGH_NAGUMO:
            return "Excitable media model derived from Hodgkin-Huxley. Produces "
                   "traveling waves, spiral waves. (FitzHugh 1961, Nagumo 1962)";
        case MODEL_LENGYEL_EPSTEIN:
            return "CIMA (Chlorite-Iodide-Malonic Acid) reaction model. First "
                   "experimental Turing patterns (1990). (Lengyel & Epstein 1991)";
        case MODEL_THOMAS:
            return "Substrate-inhibition kinetics for immobilized enzyme reactors. "
                   "Experimentally validated uricase system. (Thomas 1975)";
        default:
            return "Unknown model type.";
    }
}

int model_supports_turing(ModelType type, const ModelParams *p,
                          double Du, double Dv) {
    HomogeneousSteadyState hss;
    int ret = model_find_hss(type, p, &hss);
    if (ret != 0 || !hss.converged) return 0;

    double J[2][2];
    ret = model_jacobian(type, hss.u_star, hss.v_star, p, J);
    if (ret != 0) return 0;

    TuringConditions tc;
    turing_conditions_compute(J, Du, Dv, &tc);
    return tc.turing_unstable;
}
