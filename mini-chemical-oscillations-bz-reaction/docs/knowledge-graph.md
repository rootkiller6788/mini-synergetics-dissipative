# Knowledge Graph — mini-chemical-oscillations-bz-reaction

## Nine-Layer Knowledge Coverage

### L1: Definitions — COMPLETE
| # | Definition | Implementation |
|---|-----------|---------------|
| 1 | BZ reaction (Belousov-Zhabotinsky) | `bz_reaction.h` — `bz_state_t`, `bz_params_t` |
| 2 | FKN mechanism (Field-Koros-Noyes) | `bz_reaction.h` — `bz_fkn_species_t`, `bz_fkn_process_t` |
| 3 | Oregonator model (Field-Noyes 1974) | `bz_reaction.h` — `bz_oregonator_state_t`, `bz_oregonator_params_t` |
| 4 | Autocatalysis | `bz_reaction.h` — `BZ_FEEDBACK_POSITIVE` |
| 5 | Chemical oscillation | `bz_reaction.h` — `BZ_MODE_OSCILLATORY` |
| 6 | Bistability | `bz_bistability.h` — `BZ_MODE_BISTABLE` |
| 7 | Excitability | `bz_reaction.h` — `BZ_MODE_EXCITABLE` |
| 8 | Reaction-diffusion system | `bz_reaction_diffusion.h` — `bz_rd_1d_t`, `bz_rd_2d_t` |
| 9 | Turing instability | `bz_reaction_diffusion.h` — `bz_turing_unstable()` |
| 10 | Spiral wave | `bz_waves.h` — `bz_spiral_seed_2d()`, `bz_find_spiral_tip()` |

### L2: Core Concepts — COMPLETE
| # | Concept | Implementation |
|---|---------|---------------|
| 1 | Positive feedback (autocatalysis) | `bz_reaction.h` — Process B kinetics |
| 2 | Negative feedback (inhibition) | `bz_reaction.h` — Process C Br- production |
| 3 | Time-scale separation | `bz_oregonator.h` — `bz_oregonator_relaxation_ratio()` |
| 4 | Limit cycle oscillation | `bz_oregonator.h` — `bz_oregonator_period_approx()` |
| 5 | Hopf bifurcation | `bz_stability.h` — `bz_oregonator_hopf_test()` |
| 6 | Excitable media | `bz_reaction_diffusion.h` — `bz_barkley_params_t` |
| 7 | Wave propagation | `bz_waves.h` — `bz_wave_front_2d()`, `bz_wave_speed_field()` |
| 8 | Hysteresis | `bz_bistability.h` — `bz_hysteresis_t`, `bz_compute_hysteresis()` |
| 9 | Relaxation oscillation | `bz_oregonator.h` — `bz_oregonator_relaxation_ratio()` |
| 10 | Synchronization of oscillators | `bz_waves.h` — `bz_cross_correlation()` |

### L3: Mathematical Structures — COMPLETE
| # | Structure | Implementation |
|---|-----------|---------------|
| 1 | Stoichiometric matrix (Z^(n x m)) | `bz_reaction.h` — `bz_stoichiometric_t` |
| 2 | Rate law ODE system | `bz_reaction.c` — `bz_fkn_rhs()` |
| 3 | Reaction network graph | `bz_reaction.h` — `bz_reaction_network_t` |
| 4 | Dimensionless ODE (Oregonator) | `bz_oregonator.c` — `bz_oregonator_rhs()` |
| 5 | Jacobian matrix (3x3) | `bz_oregonator.c` — `bz_oregonator_jacobian()` |
| 6 | Characteristic polynomial | `bz_oregonator.h` — `bz_char_poly_3x3()` |
| 7 | Cubic eigenvalue problem | `bz_oregonator.c` — `solve_cubic()` |
| 8 | Finite difference Laplacian (1D/2D) | `bz_reaction_diffusion.c` — `bz_laplacian_1d/2d()` |
| 9 | Eikonal equation | `bz_waves.h` — `bz_eikonal_fast_marching()` |
| 10 | Poincare map | `bz_oscillation.h` — `bz_poincare_section()` |

### L4: Fundamental Laws & Theorems — COMPLETE
| # | Theorem | Implementation |
|---|---------|---------------|
| 1 | Hopf bifurcation theorem (Oregonator) | `bz_stability.c` — `bz_oregonator_hopf_test()`, Lean: `is_oscillatory_int` |
| 2 | Turing instability condition | `bz_reaction_diffusion.c` — `bz_turing_unstable()`, Lean: `turing_cond1-4` |
| 3 | Catalyst conservation law | `lean`: `catalyst_conservation_invariant` |
| 4 | Steady state existence (quadratic discriminant) | `lean`: `discriminant_nonneg_example` |
| 5 | Field-Koros-Noyes mechanism | `bz_reaction.c` — `bz_fkn_rhs()` |
| 6 | Poincare-Bendixson (limit cycles) | `bz_oscillation.c` — `bz_poincare_section()` |
| 7 | Linear stability (Lyapunov indirect) | `bz_stability.c` — `bz_is_linearly_stable()` |
| 8 | Floquet theory (limit cycle stability) | `bz_stability.c` — `bz_floquet_stability()` |

### L5: Algorithms & Methods — COMPLETE
| # | Algorithm | Implementation |
|---|-----------|---------------|
| 1 | RK4 fixed-step integration | `bz_integration.c` — `bz_rk4_step()` |
| 2 | RK45 adaptive step (Dormand-Prince) | `bz_integration.c` — `bz_rk45_adaptive_step()` |
| 3 | Implicit Euler (semi-implicit) | `bz_integration.c` — `bz_implicit_euler_oregonator()` |
| 4 | Cubic equation solver (trigonometric) | `bz_oregonator.c` — `solve_cubic()` |
| 5 | Eigenvalue computation (3x3) | `bz_oregonator.c` — `bz_eigenvalues_3x3()` |
| 6 | Peak detection algorithm | `bz_oscillation.c` — `bz_count_peaks()`, Lean: `peak_count` |
| 7 | DFT / Fourier analysis | `bz_oscillation.c` — `bz_dft_magnitude()` |
| 8 | Fast marching method (eikonal) | `bz_waves.c` — `bz_eikonal_fast_marching()` |
| 9 | Saddle-node bisection | `bz_bistability.c` — `bz_find_saddle_node()` |
| 10 | Basin of attraction sampling | `bz_bistability.c` — `bz_basin_grid()` |
| 11 | Forward Euler RD time-stepping | `bz_reaction_diffusion.c` — `bz_rd_1d/2d_step()` |

### L6: Canonical Problems — COMPLETE
| # | Problem | Implementation |
|---|---------|---------------|
| 1 | Oregonator oscillation | `examples/example_oregonator_oscillation.c` |
| 2 | BZ spiral wave formation | `examples/example_spiral_wave.c` |
| 3 | Bistability and hysteresis | `examples/example_bistability.c` |
| 4 | Period estimation from peaks | `tests/test_oscillation.c` |
| 5 | Stability classification | `bz_stability.c` |

### L7: Applications — PARTIAL (3/5)
| # | Application | Implementation |
|---|-------------|---------------|
| 1 | Bistable chemical switch (memory) | `bz_bistability.c` — `bz_chemical_switch()` |
| 2 | Cardiac arrhythmia modeling | `bz_waves.c` — spiral wave detection |
| 3 | Biochemical oscillation (glycolysis) | Oregonator as prototype oscillator |
| 4 | Chemical computing (logic gates) | Documented, not implemented |
| 5 | Circadian rhythm modeling | Documented, not implemented |

### L8: Advanced Topics — PARTIAL (2/5)
| # | Topic | Implementation |
|---|-------|---------------|
| 1 | Bogdanov-Takens bifurcation | `bz_stability.c` — `bz_detect_bogdanov_takens()` |
| 2 | Hopf bifurcation curve tracking | `bz_stability.c` — `bz_track_hopf_curve()` |
| 3 | Stochastic BZ reaction | Documented, not implemented |
| 4 | Coupled oscillator synchronization | Referenced in `bz_cross_correlation()` |
| 5 | Complex Ginzburg-Landau equation | Documented, not implemented |

### L9: Research Frontiers — PARTIAL (documented)
| # | Front | Status |
|---|-------|--------|
| 1 | Chemical Turing machine | Documented |
| 2 | Reaction-diffusion memory | Documented |
| 3 | BZ-based soft robotics | Documented |
| 4 | Synthetic biochemical oscillators | Documented |
