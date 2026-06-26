# Course Dependency Tree — Haken Slaving Principle

## Prerequisites (what this module depends on)

```
Linear Algebra (eigenvalues, eigenvectors, matrices)
  └── Differential Equations (ODEs, dynamical systems)
       └── Stability Theory (Lyapunov, linear stability)
            └── Bifurcation Theory (pitchfork, Hopf, saddle-node)
                 └── Synergetics (order parameters, slaving)  ← THIS MODULE
```

### Internal Dependency Graph

```
haken_core.{h,c}              ← System lifecycle, fundamental types
  ├── haken_mode_decomp.{h,c} ← Eigenvalue computation, spectral analysis
  ├── haken_dynamics.{h,c}    ← Numerical integration, coupling tensors
  ├── haken_adiabatic.{h,c}   ← Adiabatic elimination, slaving relations
  ├── haken_order_param.{h,c} ← Order parameter identification, GL theory
  ├── haken_stability.{h,c}   ← Fixed point, Jacobian, bifurcation analysis
  └── haken_slaving.lean      ← Formal verification (Lean 4)
```

### Downstream Dependencies (modules that need this)

```
mini-haken-slaving-principle (this module)
  ├── mini-nonequilibrium-phase-transitions
  │     (slaving near critical points)
  ├── mini-order-parameters-emergence
  │     (classification and identification)
  ├── mini-self-organization-far-equilibrium
  │     (slaving as mechanism of self-organization)
  └── mini-prigogine-dissipative-structures
        (thermodynamic foundation of slaving)
```

## Knowledge Layer Dependencies

| This Layer | Prerequisite Layers |
|------------|-------------------|
| L1 Definitions | — (self-contained) |
| L2 Core Concepts | L1 |
| L3 Math Structures | L1, L2 |
| L4 Fundamental Laws | L1, L2, L3 |
| L5 Algorithms | L1, L2, L3, L4 |
| L6 Canonical Problems | L1-L5 |
| L7 Applications | L1-L6 |
| L8 Advanced Topics | L1-L7 |
| L9 Research Frontiers | L1-L8 |

## Research Frontier Dependencies (L9)

- Non-equilibrium quantum phase transitions depend on: quantum mechanics, open quantum systems, Lindblad formalism
- ML for order parameter discovery depends on: manifold learning, autoencoders, diffusion maps
- Slaving in biological morphogenesis depends on: Turing pattern formation, reaction-diffusion, developmental biology
