# Course Dependency Tree — mini-chemical-oscillations-bz-reaction

## Prerequisites
```
Basic Chemistry (stoichiometry, rate laws)
  └── Chemical Kinetics (mass-action, ODEs)
       └── Nonlinear Dynamics (bifurcations, limit cycles)
            └── BZ Reaction / Oregonator [THIS MODULE]
                 ├── → mini-prigogine-dissipative-structures (Prigogine)
                 ├── → mini-turing-pattern-formation (Turing)
                 ├── → mini-self-organization-far-equilibrium
                 ├── → mini-haken-slaving-principle (Synergetics)
                 └── → mini-order-parameters-emergence
```

## Internal Dependencies
```
bz_reaction.h (core types)
  ├── bz_oregonator.h (reduced model)
  │     ├── bz_oscillation.h (analysis)
  │     └── bz_stability.h (stability)
  ├── bz_reaction_diffusion.h (RD coupling)
  │     └── bz_waves.h (wave analysis)
  └── bz_bistability.h (multistability)
```

## Downstream Modules
- mini-prigogine-dissipative-structures: BZ is a canonical dissipative structure
- mini-turing-pattern-formation: BZ spiral waves as Turing-type patterns
- mini-self-organization-far-equilibrium: BZ demonstrates self-organization
- mini-nonequilibrium-phase-transitions: Hopf bifurcation as dynamic phase transition
- mini-entropy-production-minimization: Thermodynamics of oscillating reactions

## Knowledge Level Progression
```
L1 (Definitions) → L2 (Concepts) → L3 (Math Structures)
    ↓
L4 (Theorems: Hopf, Turing, Conservation)
    ↓
L5 (Algorithms: RK4, DFT, Eigenvalues, Fast Marching)
    ↓
L6 (Canonical Problems: Oscillation, Spirals, Bistability)
    ↓
L7 (Applications: Chemical Switch, Cardiac Modeling)
    ↓
L8 (Advanced: Bogdanov-Takens, Bifurcation Tracking)
    ↓
L9 (Research Frontiers: Chemical Computing, Bio-oscillators)
```
