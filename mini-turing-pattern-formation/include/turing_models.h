/**
 * turing_models.h — Reaction Kinetics Models for Turing Pattern Formation
 *
 * Reference: Turing, A.M. (1952) "The Chemical Basis of Morphogenesis"
 *            Gray & Scott (1983) Chem. Eng. Sci. 38(1):29-43
 *            Gierer & Meinhardt (1972) Kybernetik 12(1):30-39
 *            Prigogine & Lefever (1968) J. Chem. Phys. 48(4):1695-1700
 *            Schnakenberg (1979) J. Theor. Biol. 81(3):389-400
 *            FitzHugh (1961) Biophys. J. 1(6):445-466
 *            Lengyel & Epstein (1991) Science 251(4994):650-652
 *
 * This header declares the reaction kinetics functions for 7 classical
 * reaction-diffusion models that exhibit Turing pattern formation.
 *
 * Each model provides:
 *   - Reaction rates: f(u,v), g(u,v)
 *   - Jacobian matrix: J = [[∂f/∂u, ∂f/∂v], [∂g/∂u, ∂g/∂v]]
 *   - Homogeneous steady state finder
 *   - Parameter validity checker
 *
 * Knowledge Coverage:
 *   L1: Gray-Scott, Gierer-Meinhardt, Brusselator, etc. model definitions
 *   L2: Activator-inhibitor concept, short-range activation / long-range inhibition
 *   L3: Reaction rate functions, Jacobian matrices, PDE coupling
 */

#ifndef TURING_MODELS_H
#define TURING_MODELS_H

#include "turing_common.h"

/* ===========================================================================
 * Model 0: GRAY-SCOTT (Pearson 1993 classification)
 *
 * Chemical scheme:  U + 2V → 3V    (autocatalytic)
 *                   V → P           (decay)
 *                   → U             (feed)
 *
 * Dimensionless PDEs:
 *   ∂u/∂t = F(1 - u) - u·v² + D_u ∇²u
 *   ∂v/∂t = u·v² - (F + k)v + D_v ∇²v
 *
 * Parameters: F (feed rate), k (kill rate)
 * Diffusion ratio: σ = D_v/D_u (typically σ ≫ 1 for patterns)
 *
 * Pattern types (Pearson 1993, Science 261:189-192):
 *   α: spots (F small, k small)
 *   β: stripes (F moderate)
 *   γ: labyrinthine (F larger)
 *   δ: spots replication (intermediate)
 *   ε: chaos / spatiotemporal intermittency
 *
 * Reference: Gray, P. & Scott, S.K. (1984) Chem. Eng. Sci. 39(6):1087-1097
 * ===========================================================================*/

/** Gray-Scott reaction rates: f(u,v) = F*(1-u) - u*v^2, g(u,v) = u*v^2 - (F+k)*v */
int grayscott_reaction(double u, double v, const ModelParams *p,
                       double *f_out, double *g_out);

/** Gray-Scott Jacobian at (u,v):
 *  ∂f/∂u = -F - v²,  ∂f/∂v = -2uv
 *  ∂g/∂u = v²,       ∂g/∂v = 2uv - (F+k)
 */
int grayscott_jacobian(double u, double v, const ModelParams *p,
                       double J[2][2]);

/**
 * Find homogeneous steady state for Gray-Scott.
 * Three roots possible:
 *   (u=1, v=0)          — trivial (no reaction, always stable)
 *   (u₁,v₁), (u₂,v₂)    — non-trivial steady states
 *
 * f(u,v)=0, g(u,v)=0 => F*(1-u)=u*v^2, u*v^2=(F+k)*v
 * If v≠0: v = F*(1-u)/(k) ... solve quadratic
 */
int grayscott_find_hss(const ModelParams *p, HomogeneousSteadyState *hss);

/** Validate Gray-Scott parameter ranges (F > 0, k > 0). */
int grayscott_validate_params(const ModelParams *p);


/* ===========================================================================
 * Model 1: GIERER-MEINHARDT (Activator-Inhibitor)
 *
 * Equations:
 *   ∂a/∂t = ρ_a · a²/h - μ_a·a + c_a + D_a ∇²a    (activator)
 *   ∂h/∂t = ρ_h · a² - μ_h·h + c_h + D_h ∇²h      (inhibitor)
 *
 * General form:
 *   ∂u/∂t = ρ_u · u²/(v(1 + κ_u·u²)) - μ_u·u + D_u ∇²u
 *   ∂v/∂t = ρ_v · u² - μ_v·v + D_v ∇²v
 *
 * Key requirement: D_v ≫ D_u (long-range inhibition, short-range activation)
 *
 * Reference: Gierer, A. & Meinhardt, H. (1972) Kybernetik 12(1):30-39
 *            Meinhardt, H. (1982) Models of Biological Pattern Formation
 * ===========================================================================*/

/** Gierer-Meinhardt reaction rates (simplified form). */
int giere_meinhardt_reaction(double u, double v, const ModelParams *p,
                              double *f_out, double *g_out);

/** Gierer-Meinhardt Jacobian. */
int giere_meinhardt_jacobian(double u, double v, const ModelParams *p,
                              double J[2][2]);

/** Find HSS for Gierer-Meinhardt. */
int giere_meinhardt_find_hss(const ModelParams *p, HomogeneousSteadyState *hss);

/** Validate Gierer-Meinhardt parameters. */
int giere_meinhardt_validate_params(const ModelParams *p);


/* ===========================================================================
 * Model 2: BRUSSELATOR (Prigogine & Lefever 1968)
 *
 * Chemical scheme:
 *   A → X                (constant supply)
 *   B + X → Y + D        (conversion)
 *   2X + Y → 3X          (autocatalytic trimolecular step)
 *   X → E                (removal)
 *
 * Dimensionless PDEs:
 *   ∂u/∂t = A - (B+1)u + u²v + D_u ∇²u
 *   ∂v/∂t = B·u - u²v + D_v ∇²v
 *
 * Turing condition: B > 1 + A² (supercritical Hopf bifurcation interacts)
 *
 * Reference: Prigogine, I. & Lefever, R. (1968) J. Chem. Phys. 48(4):1695-1700
 *            Nicolis, G. & Prigogine, I. (1977) Self-Organization in Nonequilibrium Systems
 * ===========================================================================*/

/** Brusselator reaction rates. */
int brusselator_reaction(double u, double v, const ModelParams *p,
                         double *f_out, double *g_out);

/** Brusselator Jacobian. */
int brusselator_jacobian(double u, double v, const ModelParams *p,
                         double J[2][2]);

/** Find HSS for Brusselator: u* = A, v* = B/A. */
int brusselator_find_hss(const ModelParams *p, HomogeneousSteadyState *hss);

/** Validate Brusselator parameters (A > 0, B > 0). */
int brusselator_validate_params(const ModelParams *p);


/* ===========================================================================
 * Model 3: SCHNAKENBERG (1979)
 *
 * Simplest realistic tri-molecular reaction scheme:
 *   ∂u/∂t = γ(a - u + u²v) + D_u ∇²u
 *   ∂v/∂t = γ(b - u²v) + D_v ∇²v
 *
 * where γ is a scaling factor, a, b are dimensionless parameters.
 * This model separates diffusion and reaction timescales cleanly.
 *
 * Reference: Schnakenberg, J. (1979) J. Theor. Biol. 81(3):389-400
 *            Murray, J.D. (2003) Mathematical Biology II, Ch. 2
 * ===========================================================================*/

/** Schnakenberg reaction rates. */
int schnakenberg_reaction(double u, double v, const ModelParams *p,
                          double *f_out, double *g_out);

/** Schnakenberg Jacobian. */
int schnakenberg_jacobian(double u, double v, const ModelParams *p,
                          double J[2][2]);

/** Find HSS for Schnakenberg: u* = a+b, v* = b/(a+b)². */
int schnakenberg_find_hss(const ModelParams *p, HomogeneousSteadyState *hss);

/** Validate Schnakenberg parameters. */
int schnakenberg_validate_params(const ModelParams *p);


/* ===========================================================================
 * Model 4: FITZHUGH-NAGUMO (Excitable Media)
 *
 * Originally a simplification of Hodgkin-Huxley neuron model.
 * In reaction-diffusion context: produces traveling waves, spiral waves.
 *
 *   ∂u/∂t = u - u³/3 - v + I_ext + D_u ∇²u
 *   ∂v/∂t = ε(u + β - γv) + D_v ∇²v
 *
 * Reference: FitzHugh, R. (1961) Biophys. J. 1(6):445-466
 *            Nagumo, J. et al. (1962) Proc. IRE 50(10):2061-2070
 * ===========================================================================*/

/** FitzHugh-Nagumo reaction rates. */
int fitzhugh_nagumo_reaction(double u, double v, const ModelParams *p,
                             double *f_out, double *g_out);

/** FitzHugh-Nagumo Jacobian. */
int fitzhugh_nagumo_jacobian(double u, double v, const ModelParams *p,
                             double J[2][2]);

/** Find HSS for FitzHugh-Nagumo. */
int fitzhugh_nagumo_find_hss(const ModelParams *p, HomogeneousSteadyState *hss);

/** Validate FitzHugh-Nagumo parameters. */
int fitzhugh_nagumo_validate_params(const ModelParams *p);


/* ===========================================================================
 * Model 5: LENGYEL-EPSTEIN (CIMA Reaction)
 *
 * Models the Chlorite-Iodide-Malonic Acid (CIMA) reaction,
 * the first experimental realization of Turing patterns (1990).
 *
 *   ∂u/∂t = a - u - 4uv/(1+u²) + D_u ∇²u
 *   ∂v/∂t = σ[b(u - uv/(1+u²))] + D_v ∇²v
 *
 * where σ = c (starch indicator complexation).
 *
 * Reference: Lengyel, I. & Epstein, I.R. (1991) Science 251(4994):650-652
 *            Castets, V. et al. (1990) Phys. Rev. Lett. 64(24):2953-2956
 * ===========================================================================*/

/** Lengyel-Epstein reaction rates (CIMA model). */
int lengyel_epstein_reaction(double u, double v, const ModelParams *p,
                             double *f_out, double *g_out);

/** Lengyel-Epstein Jacobian. */
int lengyel_epstein_jacobian(double u, double v, const ModelParams *p,
                             double J[2][2]);

/** Find HSS for Lengyel-Epstein: u* = a/5, v* = 1 + (a/5)². */
int lengyel_epstein_find_hss(const ModelParams *p, HomogeneousSteadyState *hss);

/** Validate Lengyel-Epstein parameters. */
int lengyel_epstein_validate_params(const ModelParams *p);


/* ===========================================================================
 * Model 6: THOMAS (1975) — Substrate-Inhibition Kinetics
 *
 * Experimentally motivated by the uricase-urate-O₂ reaction.
 *
 *   ∂u/∂t = a - u - ρ·R(u,v) + D_u ∇²u
 *   ∂v/∂t = α(b - v) - ρ·R(u,v) + D_v ∇²v
 *
 * where R(u,v) = uv/(1 + u + Ku²)
 *
 * Reference: Thomas, D. (1975) in "Analysis and Control of Immobilized
 *            Enzyme Systems", Springer
 *            Murray, J.D. (2003) Mathematical Biology II, §2.4
 * ===========================================================================*/

/** Thomas model reaction rates. */
int thomas_reaction(double u, double v, const ModelParams *p,
                    double *f_out, double *g_out);

/** Thomas model Jacobian. */
int thomas_jacobian(double u, double v, const ModelParams *p,
                    double J[2][2]);

/** Find HSS for Thomas model. */
int thomas_find_hss(const ModelParams *p, HomogeneousSteadyState *hss);

/** Validate Thomas model parameters. */
int thomas_validate_params(const ModelParams *p);


/* ===========================================================================
 * Universal Model Dispatcher
 *
 * Routes to the correct model-specific function based on ModelType.
 * Each returns 0 on success, -1 if model type is not implemented.
 * ===========================================================================*/

/** Universal reaction rate dispatcher. */
int model_reaction(ModelType type, double u, double v, const ModelParams *p,
                   double *f_out, double *g_out);

/** Universal Jacobian dispatcher. */
int model_jacobian(ModelType type, double u, double v, const ModelParams *p,
                   double J[2][2]);

/** Universal HSS finder dispatcher. */
int model_find_hss(ModelType type, const ModelParams *p,
                   HomogeneousSteadyState *hss);

/** Universal parameter validator dispatcher. */
int model_validate_params(ModelType type, const ModelParams *p);

/** Get model name string. */
const char *model_name(ModelType type);

/** Get model description string. */
const char *model_description(ModelType type);

/** Check if a model supports Turing instability (theoretical criterion). */
int model_supports_turing(ModelType type, const ModelParams *p,
                          double Du, double Dv);

#endif /* TURING_MODELS_H */
