# Course Tree — Order Parameters & Emergence

## Prerequisites

```
Calculus → ODE → Nonlinear Dynamics → Bifurcation Theory → Order Parameters
Statistics → Probability → Stochastic Processes → Fokker-Planck → Critical Phenomena
Linear Algebra → PCA → Order Parameter Identification
Statistical Mechanics → Phase Transitions → Landau Theory → Ginzburg-Landau
```

## Internal Dependency Tree

```
oep_core.h/c                  ← Core (vector, matrix, ODE solver, RNG)
  ├── oep_order_parameter.h/c ← Scalar/complex/multi OP dynamics
  ├── oep_slaving_principle.h/c ← Slaving, adiabatic elimination, center manifold
  ├── oep_phase_transition.h/c  ← Bifurcation diagrams, critical exponents
  ├── oep_emergence_metrics.h/c ← Shannon/MI/Φ/CE emergence metrics
  ├── oep_collective_dynamics.h/c ← Kuramoto, Vicsek, Ising, Wilson-Cowan
  ├── oep_landau_ginzburg.h/c  ← Spatial TDGL, domain walls, superconductivity
  └── oep_applications.h/c     ← Real-world applications (ENSO, brain, etc.)
```

## Learning Path

1. **L1 → L2**: Start with core definitions (oep_core.h), then understand the slaving principle
2. **L2 → L3**: Move to mathematical structures (OP equation, GL free energy)
3. **L3 → L4**: Study fundamental theorems (bifurcation, critical slowing down)
4. **L4 → L5**: Apply algorithms (PCA for OP identification, RK4 integration)
5. **L5 → L6**: Solve canonical problems (pitchfork, Kuramoto, domain walls)
6. **L6 → L7**: Extend to real applications (climate, brain, ecology)
7. **L7 → L8**: Explore advanced topics (Φ, causal emergence, Floquet)
8. **L8 → L9**: Research frontiers (optimal coarse-graining, quantum OPs)

## Concept Map

```
                    ┌──────────────────────────┐
                    │   ORDER PARAMETER (ξ)    │
                    │  dξ/dt = αξ - βξ³ + ... │
                    └────────────┬─────────────┘
                                 │
            ┌────────────────────┼────────────────────┐
            │                    │                    │
    ┌───────▼──────┐   ┌────────▼────────┐  ┌────────▼──────────┐
    │  SLAVING     │   │   BIFURCATION   │  │   PHASE TRANSITION│
    │  PRINCIPLE   │   │   DIAGRAM       │  │   CRITICAL EXP.   │
    │  τ_f ≪ τ_s   │   │   Normal Forms  │  │   Universality    │
    └───────┬──────┘   └────────┬────────┘  └────────┬──────────┘
            │                    │                    │
    ┌───────▼────────────────────▼────────────────────▼──────────┐
    │                     COLLECTIVE DYNAMICS                    │
    │  Kuramoto (sync) | Vicsek (flocking) | Ising (magnetism)   │
    └────────────────────────┬───────────────────────────────────┘
                             │
    ┌────────────────────────▼───────────────────────────────────┐
    │                 SPATIAL PATTERNS                           │
    │  Landau-Ginzburg (TDGL, domain walls, correlation length)  │
    └────────────────────────┬───────────────────────────────────┘
                             │
    ┌────────────────────────▼───────────────────────────────────┐
    │                  REAL-WORLD APPLICATIONS                   │
    │  Climate (ENSO) | Brain (Φ) | Ecology (tipping)           │
    │  Markets (regime) | Traffic (congestion) | BZ (Turing)    │
    └────────────────────────────────────────────────────────────┘
```
