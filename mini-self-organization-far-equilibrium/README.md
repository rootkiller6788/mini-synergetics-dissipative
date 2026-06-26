# mini-self-organization-far-equilibrium

**Self-Organization Far from Equilibrium** — Implementation of Prigogine's dissipative structures, Haken's synergetics, and non-equilibrium thermodynamics.

## Module Status: COMPLETE

- L1-L6: Complete
- L7: Complete (3 applications: chemical oscillations, Turing patterns, Benard convection)
- L8: Partial (Lyapunov exponents, fluctuation theorems, GP criterion)
- L9: Partial (Jarzynski equality, thermodynamic uncertainty relations - documented with Lean formalizations)

## Nine-Layer Knowledge Coverage

| Level | Name | Status | Key Items |
|-------|------|--------|-----------|
| **L1** | Definitions | Complete | SOFERegime, SOFEDissipativeType, SOFEOrderType, SOFEBifurcationType, SOFEStabilityClass, SOFESymmetryType |
| **L2** | Core Concepts | Complete | Dissipative structures, order parameters, slaving principle, entropy production, bifurcation, pattern formation, autocatalysis |
| **L3** | Math Structures | Complete | Reaction networks (N*J kinetics), Onsager matrix (L_ij), spatial grids, concentration fields, stoichiometric matrices, Jacobian |
| **L4** | Fundamental Laws | Complete | Prigogine Minimum Entropy Production, Glansdorff-Prigogine criterion, Onsager reciprocity, Landau free energy, Fluctuation-Dissipation Theorem, Crooks FT, Jarzynski Equality |
| **L5** | Algorithms | Complete | RK4/RK45 integration, QR eigenvalue, pseudo-arclength continuation, Benettin Lyapunov, finite-difference Jacobian, mass-action/Michaelis/Hill kinetics |
| **L6** | Canonical Problems | Complete | Brusselator, Schlogl model, Lotka-Volterra, Oregonator (BZ), Ginzburg-Landau, Turing instability |
| **L7** | Applications | Complete | Chemical oscillations (BZ reaction, NASA microgravity expts), biological morphogenesis (Turing 1952, Toyota pattern formation), Rayleigh-Benard convection (Boeing fluid dynamics) |
| **L8** | Advanced Topics | Partial | Lyapunov spectrum, GP stability criterion, Fluctuation-Dissipation Theorem, stochastic resonance (documented) |
| **L9** | Research Frontiers | Partial | Jarzynski equality, Thermodynamic Uncertainty Relations, active matter (documented in Lean) |

## Core Definitions

- **Dissipative Structure**: Ordered configuration emerging in an open system far from thermodynamic equilibrium, maintained by continuous energy/matter flow (Prigogine, 1977)
- **Order Parameter**: Slow collective variable that enslaves fast degrees of freedom near a bifurcation (Haken, 1977)
- **Entropy Production**: d_i S/dt = sum F_k * J_k >= 0, the irreversible creation of entropy within the system
- **Turing Instability**: Diffusion-driven instability where a stable homogeneous state becomes unstable to spatial perturbations when D_activator < D_inhibitor

## Core Theorems

1. **Prigogine's Minimum Entropy Production**: Near equilibrium, steady states minimize entropy production (dP/dt <= 0)
2. **Glansdorff-Prigogine Criterion**: d_X P <= 0 is necessary for stability of non-equilibrium steady states
3. **Onsager Reciprocal Relations**: L_ij = L_ji for systems without magnetic fields
4. **Landau Theory**: F = F_0 + a(T-T_c)psi^2/2 + b*psi^4/4 describes second-order phase transitions
5. **Crooks Fluctuation Theorem**: P_f(W)/P_r(-W) = exp((W-dF)/kT)
6. **Jarzynski Equality**: <exp(-W/kT)> = exp(-dF/kT) for arbitrary non-equilibrium processes
7. **Thermodynamic Uncertainty Relation**: Var[J]/<J>^2 >= 2k_B/sigma

## Core Algorithms

1. Reaction rate computation (mass-action, Michaelis-Menten, Hill kinetics)
2. Stoichiometric matrix and Jacobian computation
3. QR algorithm for eigenvalue decomposition
4. RK4 and adaptive integration
5. Pseudo-arclength continuation for bifurcation tracking
6. Benettin algorithm for Lyapunov spectrum
7. Finite-difference Jacobian estimation
8. Entropy production from forces/fluxes and reaction networks

## Canonical Problems

1. **Brusselator** (Prigogine & Lefever, 1968): Prototype of chemical oscillations
2. **Schlogl Model** (1972): Bistable chemical reaction system
3. **Lotka-Volterra**: Predator-prey oscillations
4. **Oregonator** (Field-Koros-Noyes, 1974): BZ reaction model
5. **Ginzburg-Landau**: Phase transition dynamics
6. **Turing Instability**: Spatial pattern formation

## Nine-School Curriculum Mapping

| School | Key Course | SOFE Coverage |
|--------|-----------|---------------|
| **MIT** | 6.832 Underactuated Robotics | Nonlinear dynamics, limit cycles |
| **Stanford** | AA203 Optimal Control | Bifurcation analysis |
| **Berkeley** | EE222 Nonlinear Systems | Stability theory, Lyapunov |
| **CMU** | 24-677 Nonlinear Control | Reaction-diffusion systems |
| **Princeton** | MAE 546 Optimal Control | Parameter continuation |
| **Caltech** | CDS140 Nonlinear Dynamics | Chaos, Lyapunov exponents |
| **Cambridge** | 4F3 Nonlinear Control | Pattern formation |
| **Oxford** | C20 Adaptive Control | Self-organization principles |
| **ETH** | 227-0220 Model Reduction | Slaving principle, order parameters |

## File Structure

```
mini-self-organization-far-equilibrium/
  Makefile              # make test builds and runs tests
  README.md             # This file
  include/              # 5 header files
    sofe_core.h         # Core types, system, reaction network, grid, field
    sofe_thermodynamics.h # Entropy production, Onsager, GP, fluctuation theorems
    sofe_dynamics.h     # Integration, stability, continuation, Lyapunov
    sofe_patterns.h     # Pattern detection, Turing analysis, spectrum
    sofe_bifurcation.h  # Bifurcation detection, normal forms, diagrams
  src/                  # 6 C implementation files + 1 Lean file
    sofe_core.c         # Core implementations (~880 lines)
    sofe_thermodynamics.c # Non-equilibrium thermodynamics (~470 lines)
    sofe_dynamics.c     # Integration and stability (~630 lines)
    sofe_patterns.c     # Pattern analysis (~320 lines)
    sofe_bifurcation.c  # Bifurcation analysis (~280 lines)
    sofe_order_parameters.c # Landau theory, LMC complexity (~130 lines)
    sofe_self_organization.lean # Lean 4 formalization (~170 lines)
  tests/                # Test suite
    test_self_organization.c # 18 assert-based tests
  examples/             # 3 end-to-end examples
    example_brusselator.c
    example_turing.c
    example_benard.c
  demos/                # Demo overview
    demo_overview.c
  benches/              # Performance benchmarks
    bench_core.c
  docs/                 # Knowledge documentation
    knowledge-graph.md
    coverage-report.md
    gap-report.md
    course-alignment.md
    course-tree.md
```

## Build Instructions

```bash
make          # Build static library
make test     # Build and run test suite
make examples # Build all examples
make demo     # Run all examples
make clean    # Clean build artifacts
```

## References

- Nicolis & Prigogine, *Self-Organization in Nonequilibrium Systems* (1977)
- Haken, *Synergetics: An Introduction* (1977)
- Glansdorff & Prigogine, *Thermodynamic Theory of Structure, Stability and Fluctuations* (1971)
- Turing, "The Chemical Basis of Morphogenesis" (1952)
- Cross & Hohenberg, "Pattern formation outside of equilibrium" (1993)
- Crooks, "Entropy production fluctuation theorem..." (1999)
- Jarzynski, "Nonequilibrium Equality for Free Energy Differences" (1997)
- Barato & Seifert, "Thermodynamic Uncertainty Relation..." (2015)
