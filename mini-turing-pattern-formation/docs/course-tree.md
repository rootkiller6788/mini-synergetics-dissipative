# Course Tree — Turing Pattern Formation

## Prerequisites
```
mini-turing-pattern-formation/
  depends on:
    ├── [Calc III]           — Partial derivatives, gradient, divergence
    ├── [Linear Algebra]     — Eigenvalues, eigenvectors, matrix analysis
    ├── [ODE Theory]         — Stability, phase plane, bifurcations
    ├── [PDE Theory]         — Diffusion equation, separation of variables
    ├── [Numerical Methods]  — Finite differences, Runge-Kutta, Newton's method
    ├── [Signal Processing]  — Fourier transform, power spectrum
    └── [Stochastic Processes] — Langevin equation, SDEs (L8)
```

## Dependency Graph (within module)

```
L1: Definitions ─────────────────────────────────────────────┐
  │                                                           │
  ├── L2: Core Concepts (diffusion-driven instability) ──────┤
  │     │                                                     │
  │     ├── L3: Math Structures (Jacobian, dispersion) ──────┤
  │     │     │                                               │
  │     │     ├── L4: Theorems (4 Turing conditions) ────────┤
  │     │     │     │                                         │
  │     │     │     ├── L5: Algorithms (solvers, analysis) ──┤
  │     │     │     │     │                                   │
  │     │     │     │     ├── L6: Problems (Gray-Scott, GM) ─┤
  │     │     │     │     │     │                             │
  │     │     │     │     │     ├── L7: Applications ────────┘
  │     │     │     │     │     │
  │     │     │     │     │     └── L8: Advanced (stochastic)
  │     │     │     │     │
  │     │     │     │     └── L9: Frontiers (documented)
```

## Knowledge Level Flow

1. **L1→L2**: Define reaction-diffusion systems, then explain Turing instability conceptually
2. **L2→L3**: Formalize the mathematics (Jacobian, eigenvalues, dispersion)
3. **L3→L4**: Prove the 4 Turing condition theorem
4. **L4→L5**: Implement numerical algorithms based on the mathematical structure
5. **L5→L6**: Apply algorithms to canonical pattern formation problems
6. **L6→L7**: Connect to real-world biological and chemical applications
7. **L7→L8**: Extend to stochastic and advanced regimes
8. **L8→L9**: Identify research frontiers

## Research Frontiers (L9)

1. **Synthetic Morphogenesis** — Engineering artificial gene circuits for Turing patterns
2. **ML-Predicted Patterns** — Using CNNs to predict pattern type from parameters
3. **Ecological Turing Patterns** — Vegetation patterns as Turing instability
4. **3D Organoid Morphogenesis** — Turing patterns in 3D tissue engineering
5. **Quantum Reaction-Diffusion** — Pattern formation in open quantum systems
