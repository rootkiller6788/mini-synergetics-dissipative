# Course Tree — mini-entropy-production-minimization

Prerequisite dependency tree for the knowledge in this module.

## Prerequisites (What You Need Before)

```
Classical Thermodynamics (L1-L3)
├── First Law (energy conservation)
├── Second Law (entropy, irreversibility)
├── Thermodynamic potentials (F, G, H)
├── Chemical potential, phase equilibrium
└── Ideal gas thermodynamics

Statistical Mechanics (L2-L4)
├── Equilibrium ensembles (microcanonical, canonical)
├── Partition function, free energy
├── Boltzmann distribution
├── Fluctuations and response
└── Stochastic processes (Langevin, Fokker-Planck)

Linear Algebra (L2-L3)
├── Matrix operations, eigenvalues
├── Quadratic forms, positive definiteness
├── LDL^T decomposition
└── Linear systems

Calculus / Differential Equations (L2-L3)
├── Multivariable calculus (gradients, Hessians)
├── ODEs (stability, phase plane)
├── PDEs (diffusion equation, reaction-diffusion)
└── Numerical methods (RK4, finite differences)

Optimization (L3-L5)
├── Constrained optimization (KKT conditions)
├── Gradient descent
└── Lagrange multipliers
```

## This Module

```
mini-entropy-production-minimization
│
├── Core Entropy Production (L1-L3)
│   ├── sigma = J dot X (bilinear form)
│   ├── sigma = X^T L X (quadratic form)
│   ├── Rayleigh dissipative function
│   └── Entropy balance equation
│
├── Onsager Relations (L1-L5)
│   ├── L_{ij} matrix (phenomenological coefficients)
│   ├── Onsager reciprocity L_{ij} = L_{ji}
│   ├── Casimir-Onsager extension (magnetic fields)
│   ├── LDL^T decomposition, linear solve
│   └── Data-driven estimation
│
├── Minimum Entropy Production (L2-L5)
│   ├── Prigogine theorem: d(sigma)/dt <= 0
│   ├── Constrained MEP optimization
│   ├── Gradient flow dynamics
│   ├── Stability analysis (Lyapunov)
│   └── Glansdorff-Prigogine criterion
│
├── Dissipative Structures (L2-L6)
│   ├── Reaction-diffusion systems
│   ├── Turing instability
│   ├── Brusselator / Oregonator models
│   ├── Pattern formation
│   └── Bifurcation analysis
│
├── Fluctuation Theorems (L2-L8)
│   ├── Stochastic thermodynamics
│   ├── Langevin equation (Brownian motion)
│   ├── Evans-Searles FT
│   ├── Jarzynski equality
│   ├── Crooks FT
│   └── Fluctuation-Dissipation Theorem
│
├── Nonlinear Thermodynamics (L2-L8)
│   ├── Nonlinear constitutive relations
│   ├── Generalized entropy production
│   ├── State-dependent Jacobian
│   ├── MaxEP vs MEPP debate
│   └── Regime classification
│
└── Applications (L6-L7)
    ├── Climate energy balance
    ├── Metabolic networks
    ├── Chemical reactor optimization
    ├── Thermoelectric generator
    └── Rayleigh-Benard convection
```

## What Comes After (Future Modules)

```
Advanced Non-Equilibrium Thermodynamics
├── Quantum thermodynamics
├── Information thermodynamics (Landauer)
├── Active matter / self-propelled particles
├── Thermodynamic uncertainty relations
└── Stochastic thermodynamics of small systems

Complex Systems / Synergetics (Haken)
├── Slaving principle, order parameters
├── Self-organization, emergence
├── Critical phenomena far from equilibrium
└── Applications to biology, sociology, economics

Computational Thermodynamics
├── Molecular dynamics entropy production
├── Monte Carlo for rare events (Jarzynski sampling)
├── Machine learning for Onsager coefficients
└── Optimal control of dissipative systems
```

## Research Frontiers (L9)

This module connects to active research in:
- Non-equilibrium steady states of quantum systems
- Thermodynamics of information processing
- Maximum entropy production in Earth systems science
- Dissipative adaptation and the origin of life
- Precision-dissipation trade-offs in biological sensors
