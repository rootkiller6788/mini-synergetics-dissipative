# Knowledge Graph — mini-entropy-production-minimization

Nine-layer knowledge coverage for the Minimum Entropy Production Principle and non-equilibrium thermodynamics.

## L1: Definitions

| Item | C Implementation | Status |
|------|-----------------|--------|
| Thermodynamic state (T,P,V,U,S,F,G,H,mu,N) | `ThermodynamicState` struct | ✓ |
| Entropy production rate sigma | `EntropyProduction` struct | ✓ |
| Thermodynamic flux J_i | `ThermodynamicFlux` struct | ✓ |
| Thermodynamic force X_i | `ThermodynamicForce` struct | ✓ |
| Chemical reaction (stoichiometry) | `ChemicalReaction` struct | ✓ |
| Onsager coefficient matrix L_{ij} | `OnsagerMatrix` struct | ✓ |
| Resistance matrix R = L^{-1} | `ResistanceMatrix` struct | ✓ |
| Non-equilibrium steady state | `SteadyState` struct | ✓ |
| Reaction-diffusion system | `ReactionDiffusionSystem` struct | ✓ |
| Spatial dissipative pattern | `DissipativePattern` struct | ✓ |
| Stochastic trajectory | `StochasticTrajectory` struct | ✓ |
| Nonlinear Onsager coefficients | `NonlinearOnsager` struct | ✓ |
| Work distribution (Jarzynski/Crooks) | `WorkDistribution` struct | ✓ |
| Entropy production histogram | `EntropyProductionHistogram` struct | ✓ |
| Stability analysis result | `StabilityResult` struct | ✓ |
| Bifurcation analysis | `BifurcationAnalysis` struct | ✓ |

## L2: Core Concepts

| Concept | Implementation | Status |
|---------|---------------|--------|
| sigma = sum J_i X_i (bilinear form) | `ep_compute_sigma()` | ✓ |
| Volumetric entropy production | `ep_volumetric_sigma()` | ✓ |
| Entropy balance (dS = d_eS + d_iS) | `ep_entropy_balance()` | ✓ |
| Second Law verification | `ep_check_second_law()` | ✓ |
| Onsager matrix allocation/ops | `onsager_alloc/set/get` | ✓ |
| MEP steady state finding | `ss_find_mep_state()` | ✓ |
| Steady-state fluxes | `ss_compute_steady_fluxes()` | ✓ |
| Gradient flow dynamics | `ss_force_evolution()` | ✓ |
| RK4 integration of gradient flow | `ss_rk4_step()` | ✓ |
| Homogeneous steady state (R-D) | `ds_find_homogeneous_steady_state()` | ✓ |
| Turing instability check | `ds_check_turing_instability()` | ✓ |
| Dispersion relation lambda(k) | `ds_dispersion_relation()` | ✓ |
| Brownian particle simulation | `ft_simulate_brownian_harmonic()` | ✓ |
| Jarzynski protocol simulation | `ft_simulate_jarzynski_protocol()` | ✓ |
| Trajectory entropy production | `ft_compute_trajectory_entropy()` | ✓ |
| Nonlinear flux computation | `nt_compute_nonlinear_fluxes()` | ✓ |
| Nonlinear entropy production | `nt_nonlinear_sigma()` | ✓ |
| Regime classification | `nt_classify_regime()` | ✓ |

## L3: Mathematical Structures

| Structure | Implementation | Status |
|-----------|---------------|--------|
| Bilinear form sigma = J dot X | `ep_compute_sigma()` | ✓ |
| Quadratic form sigma = X^T L X | `ep_quadratic_form_forces()` | ✓ |
| Rayleigh dissipative function Phi = 1/2 J^T R J | `ep_rayleigh_dissipative()` | ✓ |
| Cross entropy production | `ep_cross_entropy_production()` | ✓ |
| LDL^T decomposition | `onsager_ldlt_decompose()` | ✓ |
| QR eigenvalue algorithm | `onsager_eigenvalues()` | ✓ |
| Matrix inversion (Gauss-Jordan) | `onsager_invert()` | ✓ |
| KKT saddle-point system | `ss_lagrange_mep()` | ✓ |
| Projected gradient on constraint manifold | `ss_projected_gradient()` | ✓ |
| Reaction Jacobian (numerical) | `ds_reaction_jacobian()` | ✓ |
| Nonlinear state-dependent Jacobian | `nt_nonlinear_jacobian()` | ✓ |
| Box-Muller Gaussian sampler | `ft_rand_gaussian()` | ✓ |
| Euler-Maruyama SDE integrator | `ft_euler_maruyama_step()` | ✓ |
| Bootstrap resampling | `ft_bootstrap_resample()` | ✓ |
| Operator splitting (Strang) | `ds_rxn_diffusion_step()` | ✓ |

## L4: Fundamental Laws

| Theorem | Implementation | Verification |
|---------|---------------|-------------|
| Prigogine MEPP: d(sigma)/dt <= 0 | `ep_prigogine_dsigma_dt()` | ✓ Numerical |
| Onsager reciprocity: L_{ij} = L_{ji} | `onsager_verify_reciprocity()` | ✓ Numerical |
| Casimir-Onsager: L(B) parity | `onsager_magnetic_reversal()` | ✓ Symbolic |
| Second Law: sigma >= 0 (linear) | PSD check of L matrix | ✓ LDL^T |
| Glansdorff-Prigogine: d_X sigma <= 0 | `ss_glansdorff_prigogine_criterion()` | ✓ Numerical |
| Excess entropy as Lyapunov function | `ss_excess_entropy_production()` | ✓ Numerical |
| Evans-Searles FT: P(S)/P(-S)=exp(S/kB) | `ft_verify_evans_searles()` | ✓ Statistical |
| Jarzynski equality: <exp(-W/kT)>=exp(-dF/kT) | `ft_compute_jarzynski_deltaF()` | ✓ Statistical |
| Crooks FT: P_F(W)/P_R(-W)=exp((W-dF)/kT) | `ft_verify_crooks()` | ✓ Statistical |
| Fluctuation-Dissipation Theorem | `ft_verify_fdt()` | ✓ Numerical |
| Nonlinear Second Law: sigma(X) >= 0 | `nt_check_nonlinear_second_law()` | ✓ Grid scan |

## L5: Algorithms/Methods

| Algorithm | Function | Complexity |
|-----------|----------|------------|
| LDL^T decomposition | `onsager_ldlt_decompose()` | O(n^3) |
| Linear system solve (SPD) | `onsager_solve()` | O(n^3) |
| Matrix inversion | `onsager_invert()` | O(n^3) |
| QR eigenvalue (symmetric) | `onsager_eigenvalues()` | O(n^3)/iter |
| MEP constrained optimization | `ss_find_mep_state()` | O(n^3) |
| Lagrange multiplier method | `ss_lagrange_mep()` | O((n+m)^3) |
| Projected gradient descent | `ss_projected_gradient()` | O(n^2)/iter |
| RK4 ODE integration | `ss_rk4_step()` | O(n^2)/step |
| Relaxation simulation to steady state | `ss_simulate_relaxation()` | O(n^2 * steps) |
| Newton steady state (R-D) | `ds_find_homogeneous_steady_state()` | O(n^3)/iter |
| Strang splitting (reaction-diffusion) | `ds_rxn_diffusion_step()` | O(ng*ns)/step |
| Bifurcation continuation | `ds_track_bifurcation()` | O(steps * n^3) |
| Euler-Maruyama (SDE) | `ft_euler_maruyama_step()` | O(1)/step |
| Bootstrap resampling | `ft_bootstrap_resample()` | O(n) |
| Polynomial regression (nonlinear Onsager) | `nt_fit_nonlinear_onsager()` | O(m * n^4) |
| Nonlinear steady state relaxation | `nt_find_nonlinear_steady_states()` | O(n^3)/iter |
| MaxEP gradient ascent | `nt_find_max_entropy_state()` | O(n^3)/iter |

## L6: Canonical Problems

| Problem | Example/Function | Status |
|---------|-----------------|--------|
| Heat conduction (Fourier) | `onsager_init_fourier()` | ✓ |
| Coupled two-process transport | `onsager_init_two_process()` | ✓ |
| Multicomponent diffusion | `onsager_init_diffusion()` | ✓ |
| Chemical reaction network | `onsager_init_chemical()` | ✓ |
| Brusselator (Turing patterns) | `ds_setup_brusselator()`, `example_pattern.c` | ✓ |
| Oregonator (BZ reaction) | `ds_setup_oregonator()` | ✓ |
| Schnakenberg model | `ds_setup_schnakenberg()` | ✓ |
| Brownian particle (harmonic trap) | `ft_simulate_brownian_harmonic()`, `example_fluctuation.c` | ✓ |
| Jarzynski protocol | `ft_simulate_jarzynski_protocol()` | ✓ |
| Climate two-box model | `example_climate.c` | ✓ |

## L7: Applications

| Application | Implementation | Key Result |
|-------------|---------------|------------|
| Earth climate energy balance | `app_climate_two_box()` | MEP predicts dT ~50K (matches obs.) |
| Metabolic network (glucose) | `app_metabolic_entropy_production()` | 40% thermodynamic efficiency |
| CSTR reactor optimization | `app_cstr_optimization()` | Optimal T minimizes sigma |
| Thermoelectric generator | `app_thermoelectric_optimization()` | Optimal load matching |
| Rayleigh-Benard convection | `app_rayleigh_benard_entropy()` | sigma_conv > sigma_cond |

## L8: Advanced Topics

| Topic | Implementation | Status |
|-------|---------------|--------|
| Fluctuation theorems (Evans-Searles) | `ft_verify_evans_searles()` | ✓ Complete |
| Crooks fluctuation theorem | `ft_verify_crooks()` | ✓ Complete |
| Jarzynski equality (free energy) | `ft_compute_jarzynski_deltaF()` | ✓ Complete |
| Nonlinear irreversible thermodynamics | `nonlinear_thermo.c` | ✓ Complete |
| Maximum Entropy Production (MaxEP) | `nt_find_max_entropy_state()` | ✓ Partial |
| Dissipative adaptation / self-organization | `dissipative_systems.c` | ✓ Partial |
| Stochastic thermodynamics | `fluctuation_theorems.c` | ✓ Partial |

## L9: Research Frontiers

| Topic | Status | Notes |
|-------|--------|-------|
| Quantum thermodynamics | Documented | Landauer principle referenced |
| Information thermodynamics | Documented | Maxwell's demon, Szilard engine |
| Active matter | Documented | Non-equilibrium phase transitions |
| Thermodynamic uncertainty relations | Documented | Precision-dissipation trade-off |
| Meta-complexity & thermodynamics | Documented | Connections to computational complexity |
