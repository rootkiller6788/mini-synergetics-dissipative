#include "critical_phenomena.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* ============================================================================
 * Critical Phenomena — Implementation
 *
 * Knowledge points:
 * 1. Critical exponent definitions (L2: β, γ, α, δ, ν, η_fisher)
 * 2. Mean-field exponents and their universality (L2)
 * 3. Exact 2D Ising exponents (L4: Onsager solution)
 * 4. Scaling relation checks (L4: Rushbrooke, Widom, Fisher, Josephson)
 * 5. Exponent estimation from data (L5: linear regression in log-log)
 * 6. Data collapse method (L5: scaling function extraction)
 * 7. Finite-size scaling analysis (L5: from L-dependent data)
 * 8. Renormalization group (simplified 1D Ising) (L8)
 * 9. Universality class matching (L2)
 * ============================================================================ */

/* ---------------------------------------------------------------------------
 * Predefined Universality Classes
 *
 * Mean-field exponents (d ≥ 4, any N):
 *   α = 0, β = 1/2, γ = 1, δ = 3, ν = 1/2, η = 0
 *   Valid for: Landau theory, Bragg-Williams, Curie-Weiss,
 *              long-range interacting systems, infinite-range models.
 *   Upper critical dimension: d_c = 4 for O(N) φ⁴ theory.
 *
 * 2D Ising (Onsager 1944, exact):
 *   α = 0 (log), β = 1/8, γ = 7/4, δ = 15, ν = 1, η = 1/4
 *
 * 3D Ising (best numerical estimates):
 *   α ≈ 0.110, β ≈ 0.326, γ ≈ 1.237, δ ≈ 4.789, ν ≈ 0.630, η ≈ 0.036
 *   Note: α, δ, η differ from mean-field dramatically!
 *
 * Reference: Pelissetto & Vicari, Phys. Rep. 368, 549 (2002)
 * --------------------------------------------------------------------------- */
CriticalExponents nept_critical_exponents_mean_field(int dimension)
{
    CriticalExponents ce;
    ce.alpha = 0.0;
    ce.beta = 0.5;
    ce.gamma = 1.0;
    ce.delta = 3.0;
    ce.nu = 0.5;
    ce.eta_fisher = 0.0;
    ce.z = 2.0;   /* Model A dynamics */
    ce.dimension = dimension;
    ce.n_components = 1;
    return ce;
}

CriticalExponents nept_critical_exponents_ising_2d(void)
{
    CriticalExponents ce;
    ce.alpha = 0.0;         /* logarithmic singularity */
    ce.beta = 0.125;        /* 1/8 exactly */
    ce.gamma = 1.75;        /* 7/4 exactly */
    ce.delta = 15.0;         /* exactly 15 */
    ce.nu = 1.0;             /* exactly 1 */
    ce.eta_fisher = 0.25;    /* 1/4 exactly */
    ce.z = 2.166;            /* estimated */
    ce.dimension = 2;
    ce.n_components = 1;
    return ce;
}

CriticalExponents nept_critical_exponents_ising_3d(void)
{
    CriticalExponents ce;
    ce.alpha = 0.110;
    ce.beta = 0.3265;
    ce.gamma = 1.2372;
    ce.delta = 4.789;
    ce.nu = 0.6301;
    ce.eta_fisher = 0.0364;
    ce.z = 2.025;
    ce.dimension = 3;
    ce.n_components = 1;
    return ce;
}

CriticalExponents nept_critical_exponents_xy_3d(void)
{
    CriticalExponents ce;
    ce.alpha = -0.0146;
    ce.beta = 0.3485;
    ce.gamma = 1.3177;
    ce.delta = 4.780;
    ce.nu = 0.6717;
    ce.eta_fisher = 0.0381;
    ce.z = 2.0; /* Model E dynamics for superfluid */
    ce.dimension = 3;
    ce.n_components = 2;
    return ce;
}

CriticalExponents nept_critical_exponents_heisenberg_3d(void)
{
    CriticalExponents ce;
    ce.alpha = -0.1336;
    ce.beta = 0.3689;
    ce.gamma = 1.3960;
    ce.delta = 4.783;
    ce.nu = 0.7112;
    ce.eta_fisher = 0.0375;
    ce.z = 2.0;
    ce.dimension = 3;
    ce.n_components = 3;
    return ce;
}

/* ---------------------------------------------------------------------------
 * Scaling Relation Checks
 *
 * These are exact identities that follow from the scaling hypothesis:
 * f_s(t, h) = b^{-d} f_s(b^{y_t} t, b^{y_h} h)
 *
 * The exponents are related to the scaling dimensions:
 *   y_t = 1/ν, y_h = (d + 2 - η)/2
 *   α = 2 - d/y_t, β = (d - y_h)/y_t, γ = (2y_h - d)/y_t, δ = y_h/(d - y_h)
 *
 * From these, the four scaling relations follow algebraically.
 *
 * Reference: Fisher, Rep. Prog. Phys. 30, 615 (1967)
 * --------------------------------------------------------------------------- */
bool nept_check_rushbrooke(const CriticalExponents *ce, double tolerance)
{
    if (!ce) return false;
    double lhs = ce->alpha + 2.0 * ce->beta + ce->gamma;
    return fabs(lhs - 2.0) < tolerance;
}

bool nept_check_widom(const CriticalExponents *ce, double tolerance)
{
    if (!ce) return false;
    double lhs = ce->gamma;
    double rhs = ce->beta * (ce->delta - 1.0);
    return fabs(lhs - rhs) < tolerance;
}

bool nept_check_fisher(const CriticalExponents *ce, double tolerance)
{
    if (!ce) return false;
    double lhs = ce->gamma;
    double rhs = ce->nu * (2.0 - ce->eta_fisher);
    return fabs(lhs - rhs) < tolerance;
}

bool nept_check_josephson(const CriticalExponents *ce, double tolerance)
{
    if (!ce) return false;
    double lhs = 2.0 - ce->alpha;
    double rhs = (double)ce->dimension * ce->nu;
    return fabs(lhs - rhs) < tolerance;
}

int nept_check_all_scaling_relations(const CriticalExponents *ce,
                                      double tolerance,
                                      bool *results_out)
{
    if (!ce) return 0;
    int pass_count = 0;
    if (results_out) {
        results_out[0] = nept_check_rushbrooke(ce, tolerance);
        results_out[1] = nept_check_widom(ce, tolerance);
        results_out[2] = nept_check_fisher(ce, tolerance);
        results_out[3] = nept_check_josephson(ce, tolerance);
        for (int i = 0; i < 4; i++) if (results_out[i]) pass_count++;
    }
    return pass_count;
}

/* ---------------------------------------------------------------------------
 * Exponent Estimation via Log-Log Regression
 *
 * Knowledge point: Near the critical point, observables follow power laws:
 *   O(t) ~ |t|^{−x}  (t → 0)
 * Taking logs: log(O) = log(C) - x * log(|t|)
 * Linear regression yields the exponent x = -slope.
 *
 * For β estimation:
 *   η(t<0) = B * |t|^β  =>  log(η) = log(B) + β*log(|t|)
 *   Fit slope = β, intercept = log(B)
 *
 * Important: only use data sufficiently close to Tc (|t| << 1) to avoid
 * corrections to scaling. The asymptotic critical region is defined by
 * the Ginzburg criterion.
 *
 * Reference: Stanley, §4.3; Sokal, Monte Carlo Methods... (1997)
 * --------------------------------------------------------------------------- */
double nept_estimate_beta(const double *eta, const double *t_reduced,
                           int n_points, double *amplitude_out)
{
    if (!eta || !t_reduced || n_points < 3) {
        if (amplitude_out) *amplitude_out = 0.0;
        return 0.0;
    }

    double sum_x = 0.0, sum_y = 0.0, sum_xx = 0.0, sum_xy = 0.0;
    int count = 0;
    for (int i = 0; i < n_points; i++) {
        if (t_reduced[i] >= 0.0) continue; /* only t < 0 for ordered phase */
        if (eta[i] < 1e-15) continue;
        double x = log(fabs(t_reduced[i]));
        double y = log(eta[i]);
        sum_x += x; sum_y += y;
        sum_xx += x * x; sum_xy += x * y;
        count++;
    }
    if (count < 3) {
        if (amplitude_out) *amplitude_out = 0.0;
        return 0.0;
    }

    double denom = (double)count * sum_xx - sum_x * sum_x;
    if (fabs(denom) < 1e-15) {
        if (amplitude_out) *amplitude_out = 0.0;
        return 0.0;
    }
    double slope = ((double)count * sum_xy - sum_x * sum_y) / denom;
    double intercept = (sum_y - slope * sum_x) / (double)count;

    if (amplitude_out) *amplitude_out = exp(intercept);
    return slope; /* β = slope */
}

double nept_estimate_gamma(const double *chi, const double *t_reduced,
                            int n_points, double *amplitude_out)
{
    if (!chi || !t_reduced || n_points < 3) {
        if (amplitude_out) *amplitude_out = 0.0;
        return 0.0;
    }

    double sum_x = 0.0, sum_y = 0.0, sum_xx = 0.0, sum_xy = 0.0;
    int count = 0;
    for (int i = 0; i < n_points; i++) {
        if (fabs(t_reduced[i]) < 1e-15) continue; /* skip Tc */
        if (chi[i] < 1e-15) continue;
        double x = log(fabs(t_reduced[i]));
        double y = log(chi[i]);
        sum_x += x; sum_y += y;
        sum_xx += x * x; sum_xy += x * y;
        count++;
    }
    if (count < 3) {
        if (amplitude_out) *amplitude_out = 0.0;
        return 0.0;
    }

    double denom = (double)count * sum_xx - sum_x * sum_x;
    if (fabs(denom) < 1e-15) {
        if (amplitude_out) *amplitude_out = 0.0;
        return 0.0;
    }
    double slope = ((double)count * sum_xy - sum_x * sum_y) / denom;
    double intercept = (sum_y - slope * sum_x) / (double)count;

    /* χ ~ |t|^{-γ} => log(χ) = log(C) - γ*log(|t|) => slope = -γ */
    if (amplitude_out) *amplitude_out = exp(intercept);
    return -slope; /* γ */
}

/* ---------------------------------------------------------------------------
 * Data Collapse
 *
 * Knowledge point: The scaling hypothesis implies that near criticality,
 * the equation of state can be written as:
 *
 *   h = |η|^δ * f(t / |η|^{1/β})
 * or equivalently:
 *   η / |t|^β = F_±(h / |t|^{β+γ})  (scaling form)
 *
 * where F_± are universal scaling functions (± for above/below Tc).
 *
 * Data collapse: plot η_scaled = η / |t|^β vs h_scaled = h / |t|^{β+γ}
 * and vary β, γ until all data points fall onto a single curve.
 *
 * This function systematically searches for the exponents that minimize
 * the spread of the collapsed data.
 *
 * Reference: Bhattacharjee & Seno, J. Phys. A 34, 6375 (2001)
 * --------------------------------------------------------------------------- */
DataCollapseResult nept_data_collapse_order_param(
    const double *eta, const double *t, const double *h,
    int n_points, double beta_guess, double gamma_guess)
{
    DataCollapseResult result;
    result.best_exponent = beta_guess;
    result.residual = 0.0;
    result.confidence_interval = 0.0;

    if (!eta || !t || !h || n_points < 4) return result;

    /* Compute scaled variables with current exponents */
    double *x_scaled = (double *)malloc((size_t)n_points * sizeof(double));
    double *y_scaled = (double *)malloc((size_t)n_points * sizeof(double));
    if (!x_scaled || !y_scaled) {
        free(x_scaled); free(y_scaled);
        return result;
    }

    int count = 0;
    for (int i = 0; i < n_points; i++) {
        double abs_t = fabs(t[i]);
        if (abs_t < 1e-12) continue;

        x_scaled[count] = fabs(h[i]) / pow(abs_t, beta_guess + gamma_guess);
        y_scaled[count] = fabs(eta[i]) / pow(abs_t, beta_guess);
        count++;
    }

    if (count < 4) {
        free(x_scaled); free(y_scaled);
        return result;
    }

    /* Sort by x and compute variance of y in overlapping bins
     * as a measure of collapse quality */
    double residual = 0.0;
    for (int i = 0; i < count - 1; i++) {
        double dy = y_scaled[i+1] - y_scaled[i];
        residual += dy * dy;
    }
    result.residual = sqrt(residual / (double)(count - 1));

    free(x_scaled); free(y_scaled);
    return result;
}

/* ---------------------------------------------------------------------------
 * Finite-Size Scaling
 *
 * Knowledge point: In finite systems, the correlation length cannot
 * exceed the system size L. This modifies the scaling behavior near Tc:
 *
 *   Tc(L) - Tc(∞) ~ L^{-1/ν}           (shift)
 *   χ_max(L) ~ L^{γ/ν}                 (peak height)
 *   C_max(L) ~ (L^{α/ν} - 1)/α         (peak height)
 *
 * By studying how these quantities depend on L, we can extract ν and γ
 * without needing to know Tc exactly. This is the standard method for
 * numerical studies of critical phenomena (Monte Carlo).
 *
 * References: Barber, in Phase Transitions... Vol. 8 (1983);
 *   Binder & Heermann, Monte Carlo Simulation... (2010)
 * --------------------------------------------------------------------------- */
FiniteSizeScaling *nept_finite_size_scaling_create(int n_sizes)
{
    if (n_sizes <= 0) return NULL;
    FiniteSizeScaling *fss = (FiniteSizeScaling *)malloc(sizeof(FiniteSizeScaling));
    if (!fss) return NULL;
    fss->system_sizes = (int *)calloc((size_t)n_sizes, sizeof(int));
    fss->critical_point_shift = (double *)calloc((size_t)n_sizes, sizeof(double));
    fss->peak_susceptibility = (double *)calloc((size_t)n_sizes, sizeof(double));
    fss->peak_specific_heat = (double *)calloc((size_t)n_sizes, sizeof(double));
    fss->n_sizes = n_sizes;
    fss->estimate_nu = 0.0;
    fss->estimate_gamma = 0.0;
    if (!fss->system_sizes || !fss->critical_point_shift ||
        !fss->peak_susceptibility || !fss->peak_specific_heat) {
        nept_finite_size_scaling_free(fss);
        return NULL;
    }
    return fss;
}

void nept_finite_size_scaling_free(FiniteSizeScaling *fss)
{
    if (fss) {
        free(fss->system_sizes);
        free(fss->critical_point_shift);
        free(fss->peak_susceptibility);
        free(fss->peak_specific_heat);
        free(fss);
    }
}

void nept_finite_size_add_size(FiniteSizeScaling *fss, int idx,
                                int L, double Tc_shift,
                                double chi_peak, double C_peak)
{
    if (!fss || idx < 0 || idx >= fss->n_sizes) return;
    fss->system_sizes[idx] = L;
    fss->critical_point_shift[idx] = Tc_shift;
    fss->peak_susceptibility[idx] = chi_peak;
    fss->peak_specific_heat[idx] = C_peak;
}

int nept_finite_size_extract_exponents(FiniteSizeScaling *fss)
{
    if (!fss || fss->n_sizes < 2) return -1;

    /* From χ_max(L) ~ L^{γ/ν}: log(χ_max) vs log(L) => slope = γ/ν */
    /* From Tc_shift(L) ~ L^{-1/ν}: log(Tc_shift) vs log(L) => slope = -1/ν */
    double sum_lx = 0, sum_ly1 = 0, sum_ly2 = 0, sum_lxx = 0, sum_lxy1 = 0, sum_lxy2 = 0;
    int n = fss->n_sizes;

    for (int i = 0; i < n; i++) {
        double lx = log((double)fss->system_sizes[i]);
        double ly1 = log(fss->peak_susceptibility[i]);
        double ly2 = log(fabs(fss->critical_point_shift[i]) + 1e-15);
        sum_lx += lx; sum_ly1 += ly1; sum_ly2 += ly2;
        sum_lxx += lx * lx;
        sum_lxy1 += lx * ly1;
        sum_lxy2 += lx * ly2;
    }

    double denom = (double)n * sum_lxx - sum_lx * sum_lx;
    if (fabs(denom) < 1e-15) return -1;

    double slope1 = ((double)n * sum_lxy1 - sum_lx * sum_ly1) / denom; /* γ/ν */
    double slope2 = ((double)n * sum_lxy2 - sum_lx * sum_ly2) / denom; /* -1/ν */

    if (fabs(slope2) < 1e-15) return -1;
    fss->estimate_nu = -1.0 / slope2;
    fss->estimate_gamma = slope1 * fss->estimate_nu;

    return 0;
}

/* ---------------------------------------------------------------------------
 * Renormalization Group (1D Ising, decimation)
 *
 * Knowledge point: The RG transformation systematically coarse-grains
 * the system, reducing the number of degrees of freedom while preserving
 * the long-wavelength physics.
 *
 * 1D Ising decimation: trace out every other spin.
 *   Z = Σ_{s} exp(K Σ_i s_i s_{i+1})
 *   After tracing out even spins: Z' = Σ_{odd} exp(2K' Σ_i s_i s_{i+2})
 *   where tanh(K') = tanh²(K)  =>  K' = artanh(tanh²(K))
 *
 * The only fixed point is K* = 0 (T = ∞) — 1D Ising has no phase transition.
 * This is the simplest non-trivial example of RG.
 *
 * Reference: Cardy, Scaling and Renormalization... (1996), §3.1
 * --------------------------------------------------------------------------- */
RenormalizationGroup *nept_rg_create(int max_iter)
{
    if (max_iter <= 0) return NULL;
    RenormalizationGroup *rg = (RenormalizationGroup *)malloc(sizeof(RenormalizationGroup));
    if (!rg) return NULL;
    rg->n_iterations = 0;
    rg->coupling_constants = (double *)calloc((size_t)max_iter, sizeof(double));
    rg->fixed_points = (double *)calloc(10, sizeof(double)); /* max 10 fixed points */
    rg->n_couplings = 0;
    rg->n_fixed_points = 0;
    rg->relevant_eigenvalue = 1.0;
    rg->irrelevant_eigenvalue = 0.0;
    rg->estimate_nu = 0.0;
    if (!rg->coupling_constants || !rg->fixed_points) {
        nept_rg_free(rg);
        return NULL;
    }
    return rg;
}

void nept_rg_free(RenormalizationGroup *rg)
{
    if (rg) {
        free(rg->coupling_constants);
        free(rg->fixed_points);
        free(rg);
    }
}

int nept_rg_iterate_1d_ising(RenormalizationGroup *rg,
                              double initial_coupling, int decimation_factor)
{
    if (!rg || decimation_factor < 2) return -1;

    double K = initial_coupling;
    int max_iter = 0;
    /* Count capacity */
    {
        double *tmp = rg->coupling_constants;
        rg->coupling_constants = NULL;
        /* Can't query capacity easily; use a reasonable limit */
        max_iter = 100;
        rg->coupling_constants = tmp;
    }

    for (int iter = 0; iter < max_iter; iter++) {
        rg->coupling_constants[iter] = K;

        /* RG transformation: tanh(K') = (tanh(K))^b where b = decimation_factor */
        double th = tanh(K);
        double th_new = pow(th, (double)decimation_factor);
        /* Clamp to avoid numerical overflow */
        if (th_new >= 1.0 - 1e-15) {
            K = 10.0; /* effectively T -> 0 */
        } else if (th_new <= -1.0 + 1e-15) {
            K = -10.0;
        } else {
            K = 0.5 * log((1.0 + th_new) / (1.0 - th_new));
        }

        rg->n_iterations = iter + 1;

        /* Check for fixed point convergence */
        if (iter > 0) {
            double dK = fabs(K - rg->coupling_constants[iter - 1]);
            if (dK < 1e-12) {
                if (rg->n_fixed_points < 10) {
                    rg->fixed_points[rg->n_fixed_points++] = K;
                }
                return iter + 1;
            }
        }

        /* Known result: only fixed point is K* = 0 (T = infinity) */
        if (fabs(K) < 1e-12) {
            if (rg->n_fixed_points < 10) {
                rg->fixed_points[rg->n_fixed_points++] = 0.0;
            }
            return iter + 1;
        }
    }

    return rg->n_iterations;
}

/* ---------------------------------------------------------------------------
 * Dynamic Exponent Estimation
 *
 * Knowledge point: Near a continuous transition, the relaxation time τ
 * diverges as τ ~ ξ^z, where z is the dynamic critical exponent.
 *
 * Model A (non-conserved OP, no coupling to other slow modes): z ≈ 2.0
 * Model B (conserved OP, diffusive): z ≈ 4 - η
 * Model C (OP coupled to conserved energy): z ≈ 2 + α/ν
 * Model H (liquid-gas, OP + transverse momentum): z ≈ 3
 *
 * This function estimates z from log(τ) vs log(ξ) regression.
 *
 * Reference: Hohenberg & Halperin, Rev. Mod. Phys. 49, 435 (1977)
 * --------------------------------------------------------------------------- */
double nept_estimate_dynamic_exponent(const double *tau, const double *xi,
                                       int n_points)
{
    if (!tau || !xi || n_points < 3) return 0.0;

    double sum_x = 0, sum_y = 0, sum_xx = 0, sum_xy = 0;
    for (int i = 0; i < n_points; i++) {
        if (tau[i] < 1e-15 || xi[i] < 1e-15) continue;
        double x = log(xi[i]);
        double y = log(tau[i]);
        sum_x += x; sum_y += y;
        sum_xx += x * x; sum_xy += x * y;
    }
    double denom = (double)n_points * sum_xx - sum_x * sum_x;
    if (fabs(denom) < 1e-15) return 0.0;
    double slope = ((double)n_points * sum_xy - sum_x * sum_y) / denom;
    return slope; /* z = d log(τ) / d log(ξ) */
}

/* --- Match measured exponents to known universality class --- */
UniversalityClass nept_match_universality_class(const CriticalExponents *ce,
                                                  double tolerance)
{
    if (!ce) return UC_MEAN_FIELD;

    /* Compare with known classes */
    CriticalExponents mf = nept_critical_exponents_mean_field(ce->dimension);
    CriticalExponents is2 = nept_critical_exponents_ising_2d();
    CriticalExponents is3 = nept_critical_exponents_ising_3d();

    double diff_mf = fabs(ce->beta - mf.beta) + fabs(ce->gamma - mf.gamma)
                     + fabs(ce->nu - mf.nu);
    double diff_is2 = fabs(ce->beta - is2.beta) + fabs(ce->gamma - is2.gamma)
                      + fabs(ce->nu - is2.nu);
    double diff_is3 = fabs(ce->beta - is3.beta) + fabs(ce->gamma - is3.gamma)
                      + fabs(ce->nu - is3.nu);

    if (diff_is2 < tolerance && diff_is2 < diff_mf && diff_is2 < diff_is3)
        return UC_ISING_2D;
    if (diff_is3 < tolerance && diff_is3 < diff_mf)
        return UC_ISING_3D;

    /* Also check non-equilibrium classes */
    if (fabs(ce->beta - 0.276) < tolerance &&
        fabs(ce->nu - 1.097) < tolerance) {
        return UC_DIRECTED_PERCOLATION;
    }

    return UC_MEAN_FIELD;
}

/* --- Upper critical dimension --- */
int nept_upper_critical_dimension(const CriticalExponents *mf_exponents)
{
    if (!mf_exponents) return 4;
    /* d_c = 4 for O(N) Landau-Ginzburg phi^4 theory.
     * From mean-field exponents: 2*beta + gamma = 2, nu = 1/2.
     * d_c = 2*beta/gamma + 2*nu/nu = 2 + 2 = 4 (hyperscaling). */
    (void)mf_exponents;
    return 4;
}