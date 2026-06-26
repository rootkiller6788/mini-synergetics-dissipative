#include <stdlib.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <stdio.h>

/* ============================================================================
 * Pattern Formation via Non-Equilibrium Phase Transitions
 *
 * References:
 *   A.M. Turing - The Chemical Basis of Morphogenesis (1952)
 *   M.C. Cross & P.C. Hohenberg - Pattern Formation Outside Equilibrium (1993)
 *   R. Hoyle - Pattern Formation: An Introduction to Methods (2006)
 * ============================================================================ */

typedef struct {
    double a11, a12, a21, a22;
    double D1, D2;
    double k_min, k_max;
    int n_k;
} TuringAnalysis;

typedef struct {
    double *k_values;
    double *lambda_max;
    int n_points;
    double k_critical;
    double lambda_critical;
    int is_unstable;
    double pattern_wavelength;
} DispersionRelation;

/*
 * Allocate dispersion relation structure
 */
static DispersionRelation *turing_dispersion_create(int n_k)
{
    if (n_k < 2) return NULL;
    DispersionRelation *dr = malloc(sizeof(DispersionRelation));
    if (!dr) return NULL;
    dr->k_values = calloc((size_t)n_k, sizeof(double));
    dr->lambda_max = calloc((size_t)n_k, sizeof(double));
    if (!dr->k_values || !dr->lambda_max) {
        free(dr->k_values); free(dr->lambda_max); free(dr);
        return NULL;
    }
    dr->n_points = n_k;
    dr->k_critical = 0.0;
    dr->lambda_critical = -1e15;
    dr->is_unstable = 0;
    dr->pattern_wavelength = 0.0;
    return dr;
}

static void turing_dispersion_free(DispersionRelation *dr)
{
    if (dr) { free(dr->k_values); free(dr->lambda_max); free(dr); }
}

/*
 * Turing Instability Analysis
 *
 * Linear stability of homogeneous steady state against perturbations
 * ~ exp(lambda*t + i*k*x) in a two-component reaction-diffusion system:
 *
 *   lambda(k) = (tau(k) +/- sqrt[tau(k)^2 - 4*Delta(k)]) / 2
 *
 *   tau(k) = a11 + a22 - (D1 + D2)*k^2
 *   Delta(k) = (a11 - D1*k^2)*(a22 - D2*k^2) - a12*a21
 *
 * Turing conditions:
 *   1. Homogeneous stable: a11+a22 < 0, a11*a22 - a12*a21 > 0
 *   2. Diffusion destabilizes: Re(lambda(k)) > 0 for some k != 0
 *   3. D2/D1 must exceed a critical threshold
 *
 * Reference: Turing (1952); Murray, Mathematical Biology (2003)
 */
static DispersionRelation *compute_turing_dispersion(
    double a11, double a12, double a21, double a22,
    double D1, double D2,
    double k_min, double k_max, int n_k)
{
    if (n_k < 2) return NULL;
    DispersionRelation *dr = turing_dispersion_create(n_k);
    if (!dr) return NULL;

    double dk = (k_max - k_min) / (double)(n_k - 1);

    for (int i = 0; i < n_k; i++) {
        double k = k_min + (double)i * dk;
        double k2 = k * k;
        dr->k_values[i] = k;

        double a11k = a11 - D1 * k2;
        double a22k = a22 - D2 * k2;
        double trace_k = a11k + a22k;
        double det_k = a11k * a22k - a12 * a21;
        double disc = trace_k * trace_k - 4.0 * det_k;

        if (disc >= 0.0) {
            double sqrt_disc = sqrt(disc);
            double lam1 = 0.5 * (trace_k + sqrt_disc);
            double lam2 = 0.5 * (trace_k - sqrt_disc);
            dr->lambda_max[i] = (lam1 > lam2) ? lam1 : lam2;
        } else {
            dr->lambda_max[i] = 0.5 * trace_k;
        }

        if (dr->lambda_max[i] > dr->lambda_critical) {
            dr->lambda_critical = dr->lambda_max[i];
            dr->k_critical = k;
        }
    }

    dr->is_unstable = (dr->lambda_critical > 1e-12) ? 1 : 0;

    if (dr->is_unstable && dr->k_critical > 1e-12) {
        dr->pattern_wavelength = 2.0 * M_PI / dr->k_critical;
    }

    return dr;
}

/*
 * Swift-Hohenberg Amplitude Equation
 *
 *   dA/dt = r*A - g2*A^2 - g3*A^3
 *
 * The SH equation is the normal form for a stationary pattern-forming
 * instability (type Is). The amplitude A measures pattern strength.
 *
 * r < 0: uniform state stable
 * r > 0: pattern emerges
 * g2 != 0: quadratic term -> hexagons via transcritical bifurcation
 * g3 > 0: cubic saturation (stripe/hexagon stabilization)
 *
 * Equilibrium: A*(r - g2*A - g3*A^2) = 0
 *
 * Reference: Swift & Hohenberg, Phys. Rev. A 15, 319 (1977)
 */
static double sh_amplitude_ode(double r, double g2, double g3, double A)
{
    return r * A - g2 * A * A - g3 * A * A * A;
}

static double sh_equilibrium(double r, double g2, double g3)
{
    if (r <= 0.0 && g2 >= 0.0) return 0.0;

    if (fabs(g3) < 1e-15) {
        if (fabs(g2) < 1e-15) return (r > 0.0) ? sqrt(r) : 0.0;
        return (r > 0.0) ? r / g2 : 0.0;
    }

    double disc = g2 * g2 + 4.0 * g3 * r;
    if (disc < 0.0) return 0.0;
    double sqrt_d = sqrt(disc);
    double A1 = (-g2 + sqrt_d) / (2.0 * g3);
    double A2 = (-g2 - sqrt_d) / (2.0 * g3);
    return (A1 > 0.0) ? A1 : ((A2 > 0.0) ? A2 : 0.0);
}

/*
 * Pattern Selection: Stripes vs Hexagons
 *
 * In 2D, three modes at 120-degree angles interact:
 * Stripes (one mode):    dA/dt = r*A - g3*A^3
 * Hexagons (three equal): dA/dt = r*A + g2*A^2 - 3*g3*A^3
 *
 * Lower effective Lyapunov functional V(A) determines selection.
 *
 * Reference: Cross & Hohenberg, Rev. Mod. Phys. 65, 851 (1993)
 */
static int pattern_sel(double r, double g2, double g3,
                       double *s_amp, double *h_amp, int *pref)
{
    double As = sh_equilibrium(r, 0.0, g3);
    if (s_amp) *s_amp = As;

    double Ah = sh_equilibrium(r, g2, 3.0 * g3);
    if (h_amp) *h_amp = Ah;

    if (pref) {
        if (Ah > 1e-12 && As < 1e-12) {
            *pref = 1;
        } else if (As > 1e-12 && Ah < 1e-12) {
            *pref = 0;
        } else if (Ah > 1e-12 && As > 1e-12) {
            double Vs = -0.5 * r * As * As + 0.25 * g3 * As * As * As * As;
            double Vh = -0.5 * r * Ah * Ah - (g2 / 3.0) * Ah * Ah * Ah
                        + 0.75 * g3 * Ah * Ah * Ah * Ah;
            *pref = (Vh < Vs) ? 1 : 0;
        } else {
            *pref = -1;
        }
    }
    return 0;
}

/* === Public Interface === */

void *pattern_turing_analyze(double a11, double a12, double a21, double a22,
                              double D1, double D2,
                              double k_min, double k_max, int n_k)
{
    return compute_turing_dispersion(a11, a12, a21, a22, D1, D2,
                                      k_min, k_max, n_k);
}

void pattern_disperse_free(void *dr_ptr)
{
    turing_dispersion_free((DispersionRelation *)dr_ptr);
}

double pattern_get_k_critical(void *dr_ptr)
{
    DispersionRelation *dr = (DispersionRelation *)dr_ptr;
    return dr ? dr->k_critical : 0.0;
}

int pattern_is_turing_unstable(void *dr_ptr)
{
    DispersionRelation *dr = (DispersionRelation *)dr_ptr;
    return (dr && dr->is_unstable) ? 1 : 0;
}

double pattern_get_wavelength(void *dr_ptr)
{
    DispersionRelation *dr = (DispersionRelation *)dr_ptr;
    return dr ? dr->pattern_wavelength : 0.0;
}

double pattern_sh_equilibrium(double r, double g2, double g3)
{
    return sh_equilibrium(r, g2, g3);
}

int pattern_selection(double r, double g2, double g3,
                      double *stripe_amp, double *hex_amp,
                      int *preferred)
{
    return pattern_sel(r, g2, g3, stripe_amp, hex_amp, preferred);
}

double pattern_amplitude_growth(double r, double g2, double g3,
                                 double amplitude)
{
    return sh_amplitude_ode(r, g2, g3, amplitude);
}
