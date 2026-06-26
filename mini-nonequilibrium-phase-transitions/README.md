# mini-nonequilibrium-phase-transitions

Non-Equilibrium Phase Transitions — Theory, Computation, and Applications

## Module Status: COMPLETE

- L1-L6: Complete
- L7: Complete (3 applications: BZ reaction, Turing patterns, Schlogl bistability)
- L8: Partial (Ginzburg criterion, centre manifold reduction, stochastic resonance)
- L9: Partial (documented: directed percolation universality, non-equilibrium criticality)

---
## Knowledge Coverage Summary

| Level | Name | Status | Count |
|-------|------|--------|-------|
| L1 | Definitions | Complete | 18 struct/enum definitions |
| L2 | Core Concepts | Complete | 8 conceptual modules |
| L3 | Mathematical Structures | Complete | Landau-Ginzburg, normal forms, FPE matrices |
| L4 | Fundamental Theorems | Complete | Landau, scaling relations, Goldstone, Kramers |
| L5 | Algorithms/Methods | Complete | RK4, Crank-Nicolson, Euler-Maruyama, data collapse, FSS |
| L6 | Canonical Problems | Complete | BZ oscillator, Turing instability, Schlogl bistability |
| L7 | Applications | Complete | BZ reaction, Turing patterns, Schlogl bistability |
| L8 | Advanced Topics | Partial | Ginzburg criterion, centre manifold, stochastic resonance |
| L9 | Research Frontiers | Partial | Directed percolation, non-equilibrium universality |

---
## Core Definitions (L1)

- Phase Transition Types: First-order, Continuous, BKT, Quantum, Dynamical, Noise-induced
- System Regimes: Equilibrium, Near-equilibrium, Far-equilibrium, Driven steady, Quenched
- Ehrenfest Classification: Order of singularity in free energy derivatives
- Symmetry Types: Z2 (Ising), ZN (Potts), O(N), Translational, Gauge
- Order Parameter: Scalar, Complex, O(2), O(3), Tensor, Fourier mode, Chemical concentration
- Phase Diagram: Regions, transitions, tricritical points, critical endpoints
- NESS: Entropy production sigma, entropy flux, Lyapunov time
- Correlation Function: G(r), structure factor S(k), correlation length xi
- Landau Coefficients: a (eta^2), b (eta^4), c (eta^6)
- Bifurcation Types: Saddle-node, Pitchfork, Hopf, Homoclinic, Period-doubling
- Critical Exponents: alpha, beta, gamma, delta, nu, eta (Fisher), z (dynamic)
- Fokker-Planck, Langevin, Master equation structures
- Kramers Escape, Stochastic Resonance parameters

---
## Core Theorems (L4)

1. Landau Theory of Phase Transitions (Landau 1937)
2. Scaling Relations: Rushbrooke, Widom, Fisher, Josephson
3. Goldstone Theorem (Goldstone 1961, Nambu 1960)
4. Prigogine Minimum Entropy Production (Prigogine 1945)
5. Glansdorff-Prigogine Evolution Criterion (1971)
6. Ginzburg Criterion for mean-field validity (Ginzburg 1960)
7. Kramers Escape Rate (Kramers 1940)
8. Turing Instability condition (Turing 1952)

---
## Core Algorithms (L5)

1. RK4 Integration for ODE systems
2. Crank-Nicolson implicit FPE solver
3. Euler-Maruyama SDE integration
4. Milstein SDE scheme
5. Cholesky decomposition for correlated noise
6. Thomas tridiagonal algorithm
7. Power iteration for eigenvalues
8. Newton continuation for bifurcation branches
9. Cubic root solver (trigonometric method)
10. Log-log regression for exponent estimation
11. Data collapse with residual minimization
12. Finite-size scaling analysis
13. Benettin algorithm for Lyapunov exponents
14. Box-Muller normal random number generation

---
## Canonical Problems (L6)

1. Mean-Field Ising Transition
2. BZ Reaction (Oregonator oscillator)
3. Turing Pattern Formation
4. Schlogl Model (first-order NEPT)
5. Brusselator oscillatory instability
6. Swift-Hohenberg pattern selection
7. Critical exponent estimation from data

---
## Applications (L7)

1. BZ Chemical Oscillator — dissipative structure, temporal self-organization
2. Turing Morphogenesis — developmental biology, animal coat patterns
3. Schlogl Bistability — bistable chemical reactions, hysteresis

---
## Build and Test

make          # Build static library libnept.a
make test     # Run 39 assert-based tests (all pass)
make examples # Build 4 end-to-end examples
make demo     # Run all examples
make clean    # Remove build artifacts

---
## File Structure

include/ (6 headers): phase_transition.h, order_parameter.h, landau_theory.h, bifurcation.h, stochastic_dynamics.h, critical_phenomena.h
src/ (9 files): phase_transition.c, order_parameter.c, landau_theory.c, bifurcation.c, stochastic_dynamics.c, critical_phenomena.c, pattern_formation.c, reaction_diffusion.c, nonequilibrium_pt.lean
tests/ (1 file): test_nept.c (39 tests)
examples/ (4 files): landau_transition, bz_oregonator, turing_pattern, critical_exponents
docs/ (5 files): knowledge-graph.md, coverage-report.md, gap-report.md, course-alignment.md, course-tree.md

Total include/ + src/ lines: 4912 (exceeds 3000 minimum)

---
## Reference Textbooks

| Topic | Textbook | Author |
|-------|----------|--------|
| Synergetics | Synergetics: An Introduction | Haken (1977) |
| Dissipative Structures | Self-Organization in Nonequilibrium Systems | Nicolis & Prigogine (1977) |
| Statistical Physics | Statistical Physics, Part 1 | Landau & Lifshitz (1980) |
| Critical Phenomena | Introduction to Phase Transitions... | Stanley (1971) |
| Pattern Formation | Pattern Formation Outside Equilibrium | Cross & Hohenberg (1993) |
| Stochastic Processes | Stochastic Processes in Physics and Chemistry | van Kampen (2007) |
| Bifurcation Theory | Elements of Applied Bifurcation Theory | Kuznetsov (2004) |
| Nonlinear Dynamics | Nonlinear Dynamics and Chaos | Strogatz (2015) |
