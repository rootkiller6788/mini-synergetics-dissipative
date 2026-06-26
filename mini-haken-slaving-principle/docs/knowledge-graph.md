# Knowledge Graph — Haken Slaving Principle

## L1: Definitions (Complete)

| # | Term | C Implementation | Lean Formalization |
|---|------|-----------------|-------------------|
| 1 | Order Parameter | `Haken_OrderParam` struct | `OrderParameter` |
| 2 | Slaving Relation | `Haken_SlavingRelation` struct | `SlavingRelation` |
| 3 | Synergetic System | `Haken_System` struct | `SynergeticSystem` |
| 4 | Mode (Eigenmode) | `Haken_Mode` struct | `Mode` |
| 5 | Control Parameter | `control_param` field | `ControlParam` |
| 6 | Adiabatic Elimination | `haken_adiabatic_eliminate()` | `adiabatic_elimination_existence` |
| 7 | Ginzburg-Landau Parameters | `Haken_GLParameters` struct | (in theorem block) |

## L2: Core Concepts (Complete)

| # | Concept | Implementation |
|---|---------|---------------|
| 1 | Timescale Separation | `haken_timescale_separation()` |
| 2 | Spectral Gap | `haken_spectral_gap()` |
| 3 | Fast/Slow Mode Decomposition | `haken_classify_modes()` |
| 4 | Critical Slowing Down | `haken_relaxation_times()` |
| 5 | Adiabatic Manifold | `haken_evaluate_slaving()` |
| 6 | Self-Organization | (via order parameter emergence) |
| 7 | Bifurcation (Pitchfork/Hopf) | `haken_detect_bifurcation()` |
| 8 | Center Manifold | `center_manifold_dimension` (Lean) |

## L3: Mathematical Structures (Complete)

| # | Structure | Implementation |
|---|-----------|---------------|
| 1 | Linear Stability Matrix | `StabilityMatrix` (Lean), `L(α)` in system |
| 2 | Eigenvalue Decomposition | `haken_compute_eigenmodes()` |
| 3 | Quadratic Coupling Tensor | `Haken_QuadraticCoupling` |
| 4 | Cubic Coupling Tensor | `Haken_CubicCoupling` |
| 5 | Hessenberg Form | `haken_hessenberg()` |
| 6 | QR Decomposition | `haken_qr_decompose()` |
| 7 | Schur Form | `Haken_SchurForm` |
| 8 | Jacobian Matrix | `haken_compute_jacobian()` |

## L4: Fundamental Laws (Complete)

| # | Theorem | C Verification | Lean Theorem |
|---|---------|--------------|--------------|
| 1 | Slaving Principle (Haken 1977 §7.2) | `haken_adiabatic_eliminate()` | `adiabatic_elimination_existence` |
| 2 | Spectral Gap Theorem | `haken_spectral_gap()` | `spectral_gap_strict_separation` |
| 3 | Stable Mode Relaxation | `haken_relaxation_times()` | `stable_mode_finite_relaxation` |
| 4 | Bifurcation Necessity Condition | `haken_detect_bifurcation()` | `critical_mode_necessary_for_bifurcation` |
| 5 | GL Free Energy Minimum | `haken_gl_free_energy()` | `gl_free_energy_minimum` |
| 6 | Reduced Dynamics Preservation | `haken_reduced_dynamics()` | `reduced_dynamics_preserves_bifurcation` |
| 7 | Newton Convergence | `haken_newton_slaving()` | `newton_slaving_convergence` |

## L5: Algorithms (Complete)

| # | Algorithm | Implementation |
|---|-----------|---------------|
| 1 | QR Eigenvalue Algorithm | `haken_compute_eigenmodes()` |
| 2 | Adiabatic Elimination (Perturbation) | `haken_adiabatic_eliminate()` |
| 3 | Newton-Kantorovich Slaving | `haken_newton_slaving()` |
| 4 | Runge-Kutta 4 Integration | `haken_step_rk4()` |
| 5 | Adaptive RK4(5) (Dormand-Prince) | `haken_step_rk45_adaptive()` |
| 6 | Split-Step Integrator (Strang) | `haken_step_split()` |
| 7 | Backward Euler (Stiff) | `haken_step_backward_euler()` |
| 8 | PCA Order Parameter | `haken_identify_pca()` |
| 9 | Predictor-Corrector Continuation | `haken_continuation()` |
| 10 | Inverse Iteration (Eigenvector) | `haken_inverse_iteration()` |

## L6: Canonical Problems (Complete)

| # | Problem | Implementation |
|---|---------|---------------|
| 1 | Haken Prototype Model | `haken_create_prototype_model()` + ex1 |
| 2 | Lorenz Near First Instability | `haken_create_lorenz_instability()` + ex3 |
| 3 | Ginzburg-Landau Phase Transition | `haken_gl_fixed_points()` + ex2 |
| 4 | Rayleigh-Bénard Convection | (theorem in Lean) |
| 5 | Laser Threshold (Haken 1964) | (theorem in Lean) |
| 6 | Pattern Formation | `haken_fourier_mode_amplitude()` |

## L7: Applications (Partial+)

| # | Application | Implementation |
|---|------------|---------------|
| 1 | Laser Physics (Haken Laser Theory) | Documented in Lean L7 block + GL analysis |
| 2 | Rayleigh-Bénard Convection | `haken_create_lorenz_instability()` |
| 3 | Chemical Pattern Formation | `haken_spatial_correlation()` |
| 4 | Neural Synchronization | (documented in docs/gap-report.md) |

## L8: Advanced Topics (Partial+)

| # | Topic | Implementation |
|---|-------|---------------|
| 1 | Stochastic Slaving | `stochastic_slaving_noise_renormalization` (Lean) |
| 2 | Center Manifold Theory | `center_manifold_dimension` (Lean) |
| 3 | Time-Dependent Order Parameters | `haken_continuation()` |
| 4 | Spatially Extended Systems | `haken_spatial_correlation()` |

## L9: Research Frontiers (Partial)

| # | Topic | Status |
|---|-------|--------|
| 1 | Non-Equilibrium Quantum Phase Transitions | Documented |
| 2 | ML for Order Parameter Discovery | Documented |
| 3 | Slaving in Biological Morphogenesis | Documented |
