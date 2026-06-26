# Knowledge Graph — Turing Pattern Formation

## L1: Definitions ✓ Complete
| # | Concept | C Implementation | Lean Formalization |
|---|---------|-----------------|-------------------|
| 1 | Reaction-Diffusion System | `ReactionDiffusion2D` struct | `ReactionKinetics` structure |
| 2 | Scalar Field (2D) | `ScalarField2D` struct | — |
| 3 | Homogeneous Steady State | `HomogeneousSteadyState` struct | `HomogeneousSteadyState` structure |
| 4 | Jacobian Matrix | `JacobianFunc` typedef | `Jacobian` structure |
| 5 | Activator-Inhibitor Model | `MODEL_GIERER_MEINHARDT` | — |
| 6 | Gray-Scott Model | `MODEL_GRAY_SCOTT` | — |
| 7 | Brusselator Model | `MODEL_BRUSSELATOR` | — |
| 8 | Schnakenberg Model | `MODEL_SCHNAKENBERG` | — |
| 9 | FitzHugh-Nagumo Model | `MODEL_FITZHUGH_NAGUMO` | — |
| 10 | Lengyel-Epstein (CIMA) | `MODEL_LENGYEL_EPSTEIN` | — |
| 11 | Thomas Model | `MODEL_THOMAS` | — |
| 12 | Boundary Conditions | `BoundaryCondition` enum | — |
| 13 | Laplacian Stencil | `LaplacianStencil` enum | — |
| 14 | Pattern Types | `PatternType` enum | — |
| 15 | Turing Conditions | `TuringConditions` struct | `TuringConditions` structure |
| 16 | Dispersion Relation | `DispersionRelation` struct | — |
| 17 | Turing Space | `TuringSpace` struct | — |
| 18 | Diffusion Coefficients | `DiffusionCoefficients` struct | `DiffusionCoefficients` structure |

## L2: Core Concepts ✓ Complete
| # | Concept | Implementation |
|---|---------|---------------|
| 1 | Diffusion-Driven Instability | `turing_conditions_compute()` |
| 2 | Short-Range Activation | `giere_meinhardt_reaction()` — low Du |
| 3 | Long-Range Inhibition | `giere_meinhardt_reaction()` — high Dv |
| 4 | Pattern Selection | `turing_dispersion_compute()` — fastest growing mode |
| 5 | Turing Bifurcation | `find_turing_bifurcations()` |
| 6 | Hopf Bifurcation | `is_hopf_unstable()` |
| 7 | Turing-Hopf Interaction | `detect_turing_hopf_interaction()` |
| 8 | Pattern Wavelength | `turing_predicted_wavelength()` |
| 9 | Critical Wavenumber | `turing_critical_wavenumber()` |
| 10 | Fourier Mode Analysis | `laplacian_spectral()` |

## L3: Mathematical Structures ✓ Complete
| # | Structure | Implementation |
|---|-----------|---------------|
| 1 | 2×2 Jacobian Eigenvalues | `matrix_2x2_eigenvalues()` |
| 2 | 5-Point Laplacian Stencil | `LAPLACIAN_5PT` in `laplacian_compute()` |
| 3 | 7-Point Isotropic Stencil | `LAPLACIAN_7PT` |
| 4 | 9-Point Mehrstellen (O(h⁴)) | `LAPLACIAN_9PT` |
| 5 | Spectral Laplacian (FFT) | `laplacian_spectral()` |
| 6 | Dispersion Polynomial | `dispersion_trace`/`dispersion_det` (Lean) |
| 7 | Reaction Rate Functions | 7 models × `_reaction()` |
| 8 | Jacobian Functions | 7 models × `_jacobian()` |

## L4: Fundamental Laws ✓ Complete
| # | Theorem | Implementation |
|---|---------|---------------|
| 1 | Turing Condition 1 (tr < 0) | `tc.cond1`, test verified |
| 2 | Turing Condition 2 (det > 0) | `tc.cond2`, test verified |
| 3 | Turing Condition 3 (cross > 0) | `tc.cond3`, test verified |
| 4 | Turing Condition 4 (disc > 0) | `tc.cond4`, test verified |
| 5 | Turing Instability Theorem | `turing_unstable` in Lean (15 theorems) |
| 6 | Activator-Inhibitor Theorem | `activator_inhibitor_trace_negative` (Lean) |
| 7 | Diffusion Ratio Theorem | `cross_diffusion_inequality` (Lean) |
| 8 | Decidability of Turing Conditions | `turing_unstable_decidable` (Lean) |
| 9 | Schnakenberg Example Verified | `schnakenberg_turing_verified` (Lean) |
| 10 | Generic Example Verified | `generic_turing_example_verified` (Lean) |

## L5: Algorithms ✓ Complete
| # | Algorithm | Implementation |
|---|-----------|---------------|
| 1 | Explicit Euler Integration | `step_euler_explicit()` |
| 2 | RK2 Midpoint Integration | `step_rk2_midpoint()` |
| 3 | Classical RK4 Integration | `step_rk4()` |
| 4 | Semi-Implicit (Strang Split) | `step_semi_implicit()` |
| 5 | IMEX Euler | `step_imex_euler()` |
| 6 | ETD-RK4 | `step_etd_rk4()` |
| 7 | Newton's Method (HSS) | `newton_find_hss()` |
| 8 | Multi-Start Newton | `newton_find_all_hss()` |
| 9 | Jacobi Iteration (Helmholtz) | `jacobi_solve_helmholtz()` |
| 10 | Parameter Space Scanning | `turing_space_compute()` |
| 11 | Bisection Bifurcation Detection | `find_turing_bifurcations()` |
| 12 | Sensitivity Analysis | `turing_sensitivity_analysis()` |
| 13 | Adaptive Time-Stepping | `solver_adjust_dt()` |
| 14 | FFT (Radix-2 Cooley-Tukey) | `fft_1d()` (self-contained) |

## L6: Canonical Problems ✓ Complete
| # | Problem | Implementation |
|---|---------|---------------|
| 1 | Gray-Scott Pattern Formation | `example_grayscott.c` — spots/stripes |
| 2 | Gierer-Meinhardt Morphogenesis | `example_giere_meinhardt.c` — activator-inhibitor |
| 3 | Turing Instability Analysis | `example_stability_analysis.c` — parameter sweeps |
| 4 | Gray-Scott Pearson Classification | `field_classify_pattern()` — α,β,γ,δ,ε types |
| 5 | CIMA Reaction (Lengyel-Epstein) | `lengyel_epstein_reaction()` |
| 6 | FitzHugh-Nagumo Excitable Media | `fitzhugh_nagumo_reaction()` |

## L7: Applications ✓ Partial+ (3 applications)
| # | Application | Implementation |
|---|------------|---------------|
| 1 | Biological Morphogenesis | `example_giere_meinhardt.c` — limb/fingerprint/coat patterns |
| 2 | Chemical Pattern Formation | CIMA model — first experimental Turing patterns (1990) |
| 3 | Zebrafish/Animal Coat Patterns | Pattern classification in `field_compute_metrics()` |

## L8: Advanced Topics ✓ Partial+ (2 topics)
| # | Topic | Implementation |
|---|-------|---------------|
| 1 | Stochastic Turing Patterns | `turing_noise.c` — additive, multiplicative, colored noise |
| 2 | Euler-Maruyama SDE Integration | `euler_maruyama_step()` |
| 3 | Noise-Induced Pattern Selection | `noise_corrected_dispersion()` |
| 4 | Growing Domain Pattern Formation | Documented in course-tree.md |

## L9: Research Frontiers ✓ Partial (documented)
| # | Topic | Status |
|---|-------|--------|
| 1 | Synthetic Morphogenesis | Documented in gap-report.md |
| 2 | Machine Learning for Pattern Prediction | Documented in gap-report.md |
| 3 | Turing Patterns in Ecology | Documented in course-alignment.md |
