# mini-entropy-production-minimization

**Minimum Entropy Production Principle (MEPP) in Non-Equilibrium Thermodynamics**

Implementation of the Prigogine Minimum Entropy Production Theorem, Onsager reciprocal relations, dissipative structures, fluctuation theorems, and nonlinear irreversible thermodynamics.

---

## Module Status: COMPLETE ✅

- **L1 Definitions**: Complete (6 struct/typedef: ThermodynamicState, EntropyProduction, ThermodynamicFlux, ThermodynamicForce, ChemicalReaction, OnsagerMatrix, SteadyState, BifurcationAnalysis, DissipativePattern, StochasticTrajectory, NonlinearOnsager)
- **L2 Core Concepts**: Complete (sigma computation, entropy balance, Onsager matrix operations, MEP steady state, Turing instability, Brownian simulation)
- **L3 Mathematical Structures**: Complete (bilinear/quadratic forms, gradient flow dynamics, LDL^T decomposition, nonlinear Jacobian, dispersion relation)
- **L4 Fundamental Laws**: Complete (Prigogine theorem, Onsager reciprocity, Glansdorff-Prigogine criterion, fluctuation-dissipation theorem, Evans-Searles FT, Jarzynski equality)
- **L5 Algorithms/Methods**: Complete (LDL^T, QR eigenvalue, Newton-Raphson steady state, RK4, Euler-Maruyama, projected gradient, Langrage multiplier, Strang splitting)
- **L6 Canonical Problems**: Complete (heat conduction, coupled transport, Brusselator, Oregonator, Rayleigh-Benard, Brownian particle)
- **L7 Applications**: Complete (5 applications: climate model, metabolic network, CSTR reactor, thermoelectric generator, Rayleigh-Benard convection)
- **L8 Advanced Topics**: Partial+ (fluctuation theorems, Crooks relation, Jarzynski equality, nonlinear Onsager, Turing patterns, MEP vs MaxEP debate)
- **L9 Research Frontiers**: Partial (quantum thermodynamics connection documented, active matter references)

---

## Code Metrics

| Category | Files | Lines |
|----------|-------|-------|
| Headers (include/) | 6 | 1,988 |
| Sources (src/) | 7 | 4,171 |
| **Total (include/ + src/)** | **13** | **6,159** |
| Tests (tests/) | 3 | 400+ |
| Examples (examples/) | 3 | 300+ |

---

## Core Definitions (L1)

| Definition | Type | Description |
|------------|------|-------------|
| `ThermodynamicState` | struct | T, P, V, U, S, F, G, H, mu_k, N_k |
| `EntropyProduction` | struct | sigma_total, component breakdown, Second Law check |
| `ThermodynamicFlux` | struct | J_q, J_k, r_rho, Pi, i_electric |
| `ThermodynamicForce` | struct | X_q, X_k, A_rho, X_v, X_e |
| `ChemicalReaction` | struct | Stoichiometry, rate constants, activation energy |
| `OnsagerMatrix` | struct | L_{ij} matrix, symmetry/PSD verification |
| `SteadyState` | struct | Forces, fluxes, sigma at steady state |
| `NonlinearOnsager` | struct | L + M_{ijk} + N_{ijkl} coefficients |
| `StochasticTrajectory` | struct | Brownian particle state over time |

---

## Core Theorems (L4)

### 1. Prigogine Minimum Entropy Production Theorem (1945)
For linear constitutive relations J_i = sum_j L_{ij} X_j with time-independent boundary conditions:
```
d(sigma)/dt = 2 sum_{i,j} L_{ij} X_i (dX_j/dt) <= 0
sigma -> min at steady state (subject to constraints)
```
Implementation: `ep_prigogine_dsigma_dt()`, `ss_find_mep_state()`

### 2. Onsager Reciprocal Relations (1931)
```
L_{ij} = L_{ji}  (time-reversal symmetry)
L_{ij}(B) = epsilon_i epsilon_j L_{ji}(-B)  (Casimir-Onsager)
```
Implementation: `onsager_verify_reciprocity()`, `onsager_magnetic_reversal()`

### 3. Glansdorff-Prigogine Universal Evolution Criterion (1971)
```
d_X sigma = sum_i J_i (dX_i/dt) <= 0
```
Implementation: `ss_glansdorff_prigogine_criterion()`

### 4. Evans-Searles Fluctuation Theorem (1994)
```
P(Sigma) / P(-Sigma) = exp(Sigma / k_B)
```
Implementation: `ft_verify_evans_searles()`

### 5. Jarzynski Equality (1997)
```
<exp(-W / k_B T)> = exp(-Delta F / k_B T)
```
Implementation: `ft_compute_jarzynski_deltaF()`

### 6. Crooks Fluctuation Theorem (1999)
```
P_F(W) / P_R(-W) = exp((W - Delta F) / k_B T)
```
Implementation: `ft_verify_crooks()`

---

## Core Algorithms (L5)

| Algorithm | Function | Complexity |
|-----------|----------|------------|
| LDL^T Decomposition | `onsager_ldlt_decompose()` | O(n^3) |
| Linear System Solve | `onsager_solve()` | O(n^3) |
| Matrix Inversion | `onsager_invert()` | O(n^3) |
| QR Eigenvalue | `onsager_eigenvalues()` | O(n^3) per iteration |
| MEP Constrained Optimization | `ss_find_mep_state()` | O(n^3) |
| Lagrange Multiplier MEP | `ss_lagrange_mep()` | O((n+m)^3) |
| Projected Gradient Descent | `ss_projected_gradient()` | O(n^2) per iteration |
| RK4 Integration | `ss_rk4_step()` | O(n^2) per step |
| Euler-Maruyama (SDE) | `ft_euler_maruyama_step()` | O(1) per step |
| Strang Splitting (R-D) | `ds_rxn_diffusion_step()` | O(ng * ns) per step |
| Newton Steady State | `ds_find_homogeneous_steady_state()` | O(n^3) per iteration |
| Polynomial Regression | `nt_fit_nonlinear_onsager()` | O(m * n^4) |

---

## Canonical Problems (L6)

1. **Heat conduction between reservoirs** — `ep_compute_sigma()`, `onsager_init_fourier()`
2. **Coupled two-process transport** — `onsager_init_two_process()`, thermoelectric effect
3. **Brusselator pattern formation** — `ds_setup_brusselator()`, Turing instability analysis
4. **Oregonator (BZ reaction)** — `ds_oregonator_reactions()`, chemical oscillations
5. **Brownian particle in harmonic trap** — `ft_simulate_brownian_harmonic()`
6. **Rayleigh-Benard convection onset** — `app_rayleigh_benard_entropy()`

---

## Known Gaps (L8-L9)

| Gap | Priority | Notes |
|-----|----------|-------|
| PCP Theorem connection | Low | MEPP as variational principle for complexity |
| Quantum thermodynamics | Medium | Landauer principle, quantum fluctuation relations |
| Active matter | Medium | Non-equilibrium phase transitions in active systems |
| Thermodynamic uncertainty relations | Low | Precision-dissipation trade-off |
| Full 2D pattern simulation | Medium | Current: 1D; extension to 2D hexagonal patterns |

---

## Build and Test

```bash
make          # Build static library libepm.a
make test     # Build and run all tests (36 tests, 3 suites)
make examples # Build example programs
make clean    # Remove build artifacts
```

---

## File Structure

```
mini-entropy-production-minimization/
├── Makefile                     # Build system
├── README.md                    # This file
├── include/                     # Header files (6 files, 1988 lines)
│   ├── entropy_production.h     # Core definitions, sigma computation
│   ├── onsager_relations.h      # Onsager matrix and reciprocity
│   ├── steady_state.h           # MEP steady state analysis
│   ├── dissipative_systems.h    # Reaction-diffusion, patterns
│   ├── fluctuation_theorems.h   # FT, Jarzynski, Crooks
│   └── nonlinear_thermo.h       # Nonlinear extensions, MaxEP
├── src/                         # Source files (7 files, 4171 lines)
│   ├── entropy_production.c     # Core implementations
│   ├── onsager_relations.c      # Matrix operations, LDL^T, solve
│   ├── steady_state.c           # MEP optimization, stability
│   ├── dissipative_systems.c    # Turing patterns, BZ reaction
│   ├── fluctuation_theorems.c   # Stochastic simulation, FT verification
│   ├── nonlinear_thermo.c       # Nonlinear fluxes, regime classification
│   └── applications.c           # 5 real-world applications
├── tests/                       # Test suites (3 files)
│   ├── test_core.c              # 23 core tests
│   ├── test_onsager.c           # 6 Onsager tests
│   └── test_steady_state.c      # 7 steady-state tests
├── examples/                    # Example programs (3 files)
│   ├── example_climate.c        # Earth climate two-box model
│   ├── example_pattern.c        # Brusselator Turing patterns
│   └── example_fluctuation.c    # Fluctuation theorem verification
├── docs/
│   ├── knowledge-graph.md
│   ├── coverage-report.md
│   ├── gap-report.md
│   ├── course-alignment.md
│   └── course-tree.md
├── demos/
└── benches/
```

---

## References

### Primary
- Prigogine, I. (1945) "Moderation et transformations irreversibles des systemes ouverts", Acad. Roy. Belg., Bull. Classe Sci., 31, 600-606
- Onsager, L. (1931) "Reciprocal Relations in Irreversible Processes I, II", Phys. Rev. 37, 405; 38, 2265
- de Groot, S.R. & Mazur, P. (1984) "Non-Equilibrium Thermodynamics", Dover
- Glansdorff, P. & Prigogine, I. (1971) "Thermodynamic Theory of Structure, Stability and Fluctuations", Wiley
- Nicolis, G. & Prigogine, I. (1977) "Self-Organization in Nonequilibrium Systems", Wiley

### Fluctuation Theorems
- Evans, D.J. & Searles, D.J. (1994) Phys. Rev. E 50, 1645
- Jarzynski, C. (1997) Phys. Rev. Lett. 78, 2690
- Crooks, G.E. (1999) Phys. Rev. E 60, 2721
- Gallavotti, G. & Cohen, E.G.D. (1995) Phys. Rev. Lett. 74, 2694

### Pattern Formation
- Turing, A.M. (1952) "The Chemical Basis of Morphogenesis", Phil. Trans. R. Soc. B 237, 37
- Prigogine, I. & Lefever, R. (1968) J. Chem. Phys. 48, 1695
- Cross, M.C. & Hohenberg, P.C. (1993) Rev. Mod. Phys. 65, 851

### Synergetics
- Haken, H. (1977) "Synergetics: An Introduction", Springer
- Kondepudi, D. & Prigogine, I. (2014) "Modern Thermodynamics", 2nd ed., Wiley

### Climate & MEP
- Paltridge, G.W. (1975) "Global dynamics and climate", Q. J. R. Meteorol. Soc. 101, 475
- Kleidon, A. & Lorenz, R.D. (2005) "Non-equilibrium Thermodynamics and the Production of Entropy", Springer

---

## Knowledge Coverage Summary

| Level | Name | Rating | Score |
|-------|------|--------|-------|
| L1 | Definitions | **Complete** | 2 |
| L2 | Core Concepts | **Complete** | 2 |
| L3 | Mathematical Structures | **Complete** | 2 |
| L4 | Fundamental Laws | **Complete** | 2 |
| L5 | Algorithms/Methods | **Complete** | 2 |
| L6 | Canonical Problems | **Complete** | 2 |
| L7 | Applications | **Complete** (5) | 2 |
| L8 | Advanced Topics | **Partial+** (3/6) | 1 |
| L9 | Research Frontiers | **Partial** (documented) | 1 |
| **Total** | | | **16/18** |

**Rating: COMPLETE ✅** (≥16/18, L1≠Missing, L4≠Missing, ≥6 layers Complete)
