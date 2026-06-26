# mini-prigogine-dissipative-structures

**Ilya Prigogine's Theory of Dissipative Structures -- C Implementation**

---

## Module Status: COMPLETE

- L1-L6: Complete
- L7: Partial (3 applications documented, 2 with code)
- L8: Partial (stochastic resonance, CGLE, Game of Life, fuzzy entropy)
- L9: Partial (documented only)

---

## Nine-Level Knowledge Coverage

| Level | Name | Status | Key Items |
|-------|------|--------|-----------|
| **L1** | Definitions | **Complete** | Dissipative structure, thermodynamic regime, branch state, entropy balance, Onsager coefficients, bifurcation types |
| **L2** | Core Concepts | **Complete** | Far-from-equilibrium, self-organization, order through fluctuations, entropy exchange, symmetry breaking |
| **L3** | Math Structures | **Complete** | State vectors, matrices, Jacobian, reaction-diffusion PDE, Brusselator ODE, GP stability criterion |
| **L4** | Fundamental Laws | **Complete** | Prigogine Minimum Entropy Production (1945), Glansdorff-Prigogine Criterion (1954/1971), Onsager Reciprocal Relations (1931), Second Law for Open Systems, Turing Instability (1952) |
| **L5** | Algorithms | **Complete** | RK4 integration, finite-difference Laplacian, eigenvalue computation (2x2, 3x3), Lyapunov exponent, Newton-Raphson, Gillespie SSA |
| **L6** | Canonical Problems | **Complete** | Brusselator limit cycle, bifurcation diagram, entropy production analysis, Benard convection onset |
| **L7** | Applications | **Partial** | BZ reaction (NASA ISS, 2020), climate energy balance, glycolytic oscillations, Toyota lean production, Fukushima safety, smart grid |
| **L8** | Advanced Topics | **Partial** | Stochastic resonance, Complex Ginzburg-Landau, Game of Life, Extended Irreversible Thermodynamics, fuzzy entropy |
| **L9** | Research Frontiers | **Partial** | Quantum dissipative structures, information thermodynamics, meta-complexity (documented) |

---

## Core Definitions (L1)

- **Dissipative Structure**: Ordered state maintained by continuous energy/matter dissipation far from thermodynamic equilibrium (Prigogine, 1969)
- **Thermodynamic Branch**: The continuous extension of equilibrium behavior to near-equilibrium; becomes unstable at critical bifurcation point
- **Entropy Balance**: dS = d_i S + d_e S, where d_i S >= 0 (Second Law) and d_e S may be negative (negentropy import)
- **Bifurcation Point**: Critical parameter value where thermodynamic branch loses stability and new structures emerge
- **Onsager Coefficients**: L_ij relating thermodynamic fluxes to forces: J_i = sum_j L_ij X_j, with L_ij = L_ji near equilibrium

## Core Theorems (L4)

| Theorem | Formula | Reference |
|---------|---------|-----------|
| Minimum Entropy Production | dP/dt <= 0 near equilibrium | Prigogine (1945) |
| Glansdorff-Prigogine Stability | delta^2 P >= 0 => stable | Glansdorff & Prigogine (1954, 1971) |
| Onsager Reciprocity | L_ij = L_ji | Onsager (1931) |
| Second Law (Open Systems) | d_i_S >= 0, d_e_S any sign | Prigogine (1947) |
| Turing Instability | D_act << D_inh => spatial patterns | Turing (1952) |

## Core Algorithms (L5)

- **RK4 Integration**: 4th-order Runge-Kutta for ODE integration (Brusselator, CGLE)
- **Eigenvalue Analysis**: 2x2 quadratic, 3x3 Cardano cubic formula
- **Finite Difference Laplacian**: 1D/2D with Neumann boundary conditions
- **Lyapunov Exponent**: Rosenstein nearest-neighbor divergence algorithm
- **Gillespie SSA**: Exact stochastic simulation of chemical master equation
- **Stochastic Resonance**: Euler-Maruyama for double-well + periodic drive + noise

## Classical Problems (L6)

1. **Brusselator Limit Cycle**: Hopf bifurcation at b_c = 1 + a^2, stable oscillations for b > b_c
2. **Bifurcation Diagram**: Parameter scan showing transition from thermodynamic branch to dissipative structure
3. **Entropy Production Analysis**: Verifying Second Law and minimum entropy production theorem
4. **Turing Pattern Formation**: Diffusion-driven instability in reaction-diffusion systems

## Nine-School Curriculum Mapping

| School | Course | Relevance |
|--------|--------|-----------|
| **MIT** | 2.050J Nonlinear Dynamics | Bifurcation, chaos, limit cycles |
| **Stanford** | ME 331B Nonequilibrium Thermodynamics | Extended irreversible thermodynamics |
| **Berkeley** | Physics 221B Nonequilibrium Statistical Mechanics | Fluctuation theorems |
| **CMU** | 24-627 Molecular Simulation | Stochastic chemical kinetics (Gillespie) |
| **Princeton** | MAE 556 Nonlinear Dynamics | Bifurcation theory, normal forms |
| **Caltech** | ChE 126 Chemical Reaction Engineering | BZ reaction, chemical oscillations |
| **Cambridge** | Part II Theoretical Physics | Complex systems, pattern formation |
| **Oxford** | C21 Nonequilibrium Physics | Dissipative structures, self-organization |
| **ETH** | 402-0810 Complex Systems | Synergetics, order parameters |

## Quick Start

```bash
make          # Build the library
make test     # Run all tests
make examples # Build all examples
./build/example_brusselator   # Brusselator limit cycle demo
./build/example_bifurcation   # Bifurcation diagram demo
./build/example_entropy       # Entropy production demo
```

## Reference Textbooks

- Nicolis & Prigogine (1977). *Self-Organization in Non-Equilibrium Systems*. Wiley.
- Glansdorff & Prigogine (1971). *Thermodynamic Theory of Structure, Stability and Fluctuations*. Wiley-Interscience.
- Prigogine & Stengers (1984). *Order Out of Chaos*. Bantam Books.
- Haken (1977). *Synergetics: An Introduction*. Springer.
- Turing (1952). *The Chemical Basis of Morphogenesis*. Phil. Trans. R. Soc. B.

## File Structure

```
mini-prigogine-dissipative-structures/
  Makefile               # Build system
  README.md              # This file
  include/               # 6 header files
  src/                   # 8 C source files + 1 Lean file
  tests/                 # 4 test files (all passing)
  examples/              # 3 end-to-end examples
  docs/                  # Knowledge documentation (5 files)
  benches/               # Performance benchmarks
  demos/                 # Demonstration programs
```
