# Knowledge Graph — Order Parameters & Emergence

## L1: Definitions
| Entry | C Struct | Lean Definition | Status |
|-------|----------|-----------------|--------|
| Order parameter (scalar) | `OEPScalarOP` | `OrderParameter` | Complete |
| Order parameter (complex) | `OEPComplexOP` | - | Complete |
| Order parameter (multi-component) | `OEPMultiOP` | - | Complete |
| Emergence metrics | `OEPEmergenceResult` | `EmergenceMetrics` | Complete |
| Phase transition type | `OEPPhaseTransitionType` | `PhaseTransitionType` | Complete |
| Slaving system | `OEPSlavingSystem` | `SlavingSystem` | Complete |
| Bifurcation diagram | `OEPBifDiag` | - | Complete |

## L2: Core Concepts
| Entry | Implementation | Status |
|-------|---------------|--------|
| Slaving (subordination) principle | `oep_slaving_principle.c` | Complete |
| Adiabatic elimination | `oep_slaving_adiabatic_elimination` | Complete |
| Time-scale separation | `oep_slaving_check_timescale_separation` | Complete |
| Center manifold reduction | `oep_slaving_center_manifold_h2/h3` | Complete |
| Critical slowing down | `oep_slaving_critical_slowing` | Complete |
| Circular causality | `oep_op_circular_causality` | Complete |
| Self-organization | `oep_self_organization_index` | Complete |
| Thermodynamic branch instability | `oep_phase_transition.c` | Complete |

## L3: Mathematical Structures
| Entry | Implementation | Status |
|-------|---------------|--------|
| Order parameter ODE (dξ/dt = αξ - βξ³) | `oep_op_scalar_rhs` | Complete |
| Landau potential (V = -αξ²/2 + βξ⁴/4) | `oep_op_scalar_potential` | Complete |
| Complex Ginzburg-Landau | `oep_op_complex_rhs` | Complete |
| Time-dependent GL (TDGL) | `oep_lg_local_rhs` | Complete |
| Fokker-Planck stationary distribution | `oep_op_scalar_stationary_pdf` | Complete |
| Stochastic differential equation (Euler-Maruyama) | `oep_op_scalar_evolve_stochastic` | Complete |
| Free energy functional (spatial) | `oep_field_free_energy` | Complete |
| Bifurcation normal forms | `oep_slaving_normal_form` | Complete |
| Floquet exponents | `oep_slaving_floquet_exponents` | Complete |

## L4: Fundamental Laws/Theorems
| Theorem | C Implementation | Lean Statement | Status |
|---------|-----------------|----------------|--------|
| Bifurcation point (α=0) | `oep_op_scalar_equilibria` | `bifurcation_condition` | Complete |
| Equilibrium multiplicity | `oep_op_scalar_equilibria` | `equilibrium_count` | Complete |
| Stability exchange | `oep_op_scalar_stability` | `stability_exchange` | Complete |
| Critical slowing down | `oep_slaving_critical_slowing` | `critical_slowing_down` | Complete |
| Adiabatic elimination validity | `oep_slaving_validate_reduction` | `adiabatic_elimination_valid` | Complete |
| Kuramoto transition (K_c = 2γ) | `oep_kura_critical_coupling` | `kuramoto_transition` | Complete |
| Correlation length divergence | `oep_lg_correlation_length` | `correlation_length_divergence` | Complete |
| Rushbrooke inequality | `oep_pt_check_rushbrooke` | - | Complete |
| Widom scaling | `oep_pt_check_widom` | - | Complete |
| Fisher scaling | `oep_pt_check_fisher` | - | Complete |
| Ginzburg criterion | `oep_pt_ginzburg_criterion` | - | Complete |
| Maxwell construction | `oep_pt_maxwell_construction` | - | Complete |

## L5: Algorithms/Methods
| Algorithm | Implementation | Status |
|-----------|---------------|--------|
| PCA for OP identification | `oep_op_pca` | Complete |
| Order parameter detection from data | `oep_op_identify_from_data` | Complete |
| Runge-Kutta 4 integration | `oep_ode_rk4_step` | Complete |
| Power iteration (eigenvalues) | `oep_op_multi_normal_modes` | Complete |
| Metropolis Monte Carlo (Ising) | `oep_ising_metropolis_step` | Complete |
| Kuramoto mean-field integration | `oep_kura_step` | Complete |
| Vicsek alignment dynamics | `oep_vicsek_step` | Complete |
| Binder cumulant (T_c estimation) | `oep_ising_binder_cumulant` | Complete |
| Hurst exponent (R/S analysis) | `oep_app_market_regime_analyze` | Complete |
| Lyapunov exponent (Rosenstein) | `oep_ts_lyapunov_exponent` | Complete |

## L6: Canonical Problems
| Problem | Example/Implementation | Status |
|---------|----------------------|--------|
| Pitchfork bifurcation | `example1_pitchfork_bifurcation.c` | Complete |
| Kuramoto synchronization | `example2_kuramoto_sync.c` | Complete |
| Landau-Ginzburg domain walls | `example3_landau_ginzburg_domain.c` | Complete |
| Superconductivity (Type I/II) | `oep_lg_is_type_II` | Complete |
| Ising magnetization curve | `oep_ising_magnetization_curve` | Complete |
| Wilson-Cowan neural oscillations | `oep_wc_detect_oscillation` | Complete |
| Structural phase transition | `oep_lg_strain_equilibria` | Complete |

## L7: Applications
| Application | Implementation | Domain Keywords | Status |
|-------------|---------------|-----------------|--------|
| ENSO climate detection | `oep_app_enso_analyze` | climate, SST, El Niño | Complete |
| Brain coherence (consciousness) | `oep_app_brain_coherence_analyze` | EEG, neural, consciousness | Complete |
| Opinion polarization | `oep_app_opinion_polarization_analyze` | social dynamics, polarization | Complete |
| Ecosystem resilience / tipping | `oep_app_ecosystem_resilience_analyze` | ecology, resilience, tipping | Complete |
| Market regime detection | `oep_app_market_regime_analyze` | finance, Hurst, market | Complete |
| Traffic flow phase transition | `oep_app_traffic_flow_analyze` | traffic, congestion, phase | Complete |
| BZ pattern formation | `oep_app_bz_pattern_analyze` | BZ, Turing, pattern | Complete |

## L8: Advanced Topics
| Topic | Implementation | Status |
|-------|---------------|--------|
| Center manifold reduction | `oep_slaving_center_manifold_h2/h3` | Complete |
| Integrated information Φ | `oep_em_integrated_information_bipartition` | Complete |
| Causal emergence | `oep_em_causal_emergence` | Complete |
| Floquet theory (periodic) | `oep_slaving_floquet_exponents` | Complete |
| Transfer entropy | `oep_em_transfer_entropy` | Complete |
| Stochastic resonance (noise) | `oep_op_scalar_evolve_stochastic` | Complete |
| Renormalization (finite-size scaling) | `oep_pt_finite_size_scaling` | Complete |
| Multiscale complexity | `oep_em_multiscale_complexity` | Complete |
| Ginzburg criterion (mean-field breakdown) | `oep_pt_ginzburg_criterion` | Complete |

## L9: Research Frontiers
| Topic | Documentation | Status |
|-------|--------------|--------|
| Optimal coarse-graining (causal emergence) | `oep_em_optimal_causal_emergence` | Partial |
| Structure factor analysis | `oep_structure_factor` | Partial |
| Machine learning for OP discovery | documented in README | Partial |
| Quantum order parameters (superconductivity) | `oep_lg_superconductor_gap` | Partial |
| Topological order parameters | documented | Partial |
| Active matter order parameters | referenced in Vicsek model | Partial |
| Meta-complexity & emergence | `oep_em_analyze_emergence` | Partial |

## Summary
- **L1 Definitions**: Complete (7 entries)
- **L2 Core Concepts**: Complete (8 entries)
- **L3 Mathematical Structures**: Complete (9 entries)
- **L4 Fundamental Laws**: Complete (12 entries)
- **L5 Algorithms/Methods**: Complete (10 entries)
- **L6 Canonical Problems**: Complete (7 entries)
- **L7 Applications**: Complete (7 entries)
- **L8 Advanced Topics**: Complete (9 entries)
- **L9 Research Frontiers**: Partial (7 entries)
