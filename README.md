# Mini Synergetics and Dissipative Structures

A collection of **from-scratch, zero-dependency C implementations** of synergetics and dissipative structure theory — the study of self-organization, order parameter emergence, and pattern formation in systems driven far from thermodynamic equilibrium. Each module brings the foundational works of Haken (Synergetics), Prigogine (Dissipative Structures), Turing (Morphogenesis), and Landau (Phase Transitions) to life as runnable, self-contained C code.

## Sub-Modules

| Sub-Module | Topics | Key Courses |
|------------|--------|-------------|
| [mini-chemical-oscillations-bz-reaction](mini-chemical-oscillations-bz-reaction/) | Belousov-Zhabotinsky reaction, Field-Koros-Noyes (FKN) mechanism, Oregonator model, mass-action kinetics, Hopf bifurcation, limit cycle stability, bistability analysis, spiral waves, target patterns, reaction-diffusion coupling, Turing instability in BZ systems | MIT 5.60, MIT 18.385J |
| [mini-entropy-production-minimization](mini-entropy-production-minimization/) | Entropy production rate σ, Prigogine Minimum Entropy Production Principle (MEPP), Onsager reciprocal relations, fluctuation theorems (Crooks, Jarzynski), nonlinear irreversible thermodynamics, Glansdorff-Prigogine stability criterion, steady-state analysis, dissipative systems thermodynamics | MIT 10.40, Stanford ChemE 340 |
| [mini-haken-slaving-principle](mini-haken-slaving-principle/) | Haken's slaving principle (Versklavungsprinzip), adiabatic elimination of fast variables, mode decomposition by timescale hierarchy, spectral gap condition, order parameter identification (spectral/statistical/physical methods), nonlinear synergetic dynamics, stability and bifurcation near critical points | Santa Fe Institute CSSS, MIT 18.385J |
| [mini-nonequilibrium-phase-transitions](mini-nonequilibrium-phase-transitions/) | Bifurcation theory as nonequilibrium phase transitions, Landau-Ginzburg mean-field theory, critical phenomena (power-law singularities, universality classes), order parameters and correlation functions, stochastic dynamics (Master equation, Fokker-Planck, Langevin), symmetry breaking far from equilibrium | MIT 8.334, Cornell PHYS 7654 |
| [mini-order-parameters-emergence](mini-order-parameters-emergence/) | Order parameter dynamics (dξ/dt = αξ − βξ³ + noise), emergence quantification metrics (excess entropy, ε-machine complexity), collective dynamics (Kuramoto synchronization, Vicsek flocking), Haken's slaving principle, Landau-Ginzburg free energy, nonequilibrium phase transitions, ENSO application | Santa Fe Institute CSSS, MIT 8.334 |
| [mini-prigogine-dissipative-structures](mini-prigogine-dissipative-structures/) | Prigogine-Lefever Brusselator model, reaction-diffusion PDE systems, linear/nonlinear stability analysis, bifurcation detection (saddle-node, pitchfork, Hopf, transcritical), entropy production in dissipative structures, Onsager reciprocal relations, Glansdorff-Prigogine stability criterion, extended irreversible thermodynamics | MIT 10.40, MIT 5.60 |
| [mini-self-organization-far-equilibrium](mini-self-organization-far-equilibrium/) | Self-organization mechanisms far from equilibrium, bifurcation analysis (normal forms, codimension), nonlinear dynamics integration, spatiotemporal pattern classification (stripes, spots, hexagons, spirals, traveling waves), nonequilibrium thermodynamics (entropy, free energy, chemical potential), Prigogine-Haken-Nicolis synthesis | Santa Fe Institute CSSS, MIT 18.385J |
| [mini-turing-pattern-formation](mini-turing-pattern-formation/) | Turing instability (activator-inhibitor mechanism), 7 classical reaction-diffusion models (Gierer-Meinhardt, Gray-Scott, Brusselator, Schnakenberg, FitzHugh-Nagumo, Lengyel-Epstein, Thomas), dispersion relation λ(k), Turing space construction, numerical PDE solvers (explicit/implicit/RK4), pattern classification (spots, stripes, labyrinths) | MIT 18.354J, MIT 18.385J |

## Design Philosophy

- **Zero external dependencies** — pure C (C99/C11), only `libc` and `libm`
- **Self-contained modules** — each directory has its own `Makefile`, `include/`, `src/`, `examples/`, `demos/`, `tests/`
- **Theory-to-code mapping** — every struct and function maps to Haken, Prigogine, Turing, or Landau foundational works
- **From-scratch implementations** — no existing libraries; each concept is built from thermodynamic, synergetic, and reaction-diffusion primitives

## Building

Each module is standalone. Navigate to a module directory and run:

```bash
cd mini-chemical-oscillations-bz-reaction
make all    # build everything
make test   # run tests
```

Requires **GCC** and **GNU Make**.

## Project Structure

```
mini-synergetics-dissipative/
├── mini-chemical-oscillations-bz-reaction/  # BZ reaction, Oregonator, FKN mechanism, spiral waves
├── mini-entropy-production-minimization/    # MEPP, Onsager relations, fluctuation theorems
├── mini-haken-slaving-principle/            # Slaving principle, adiabatic elimination, order parameters
├── mini-nonequilibrium-phase-transitions/   # Landau theory, critical phenomena, stochastic dynamics
├── mini-order-parameters-emergence/         # Order parameter dynamics, emergence metrics, collective dynamics
├── mini-prigogine-dissipative-structures/   # Brusselator, reaction-diffusion, bifurcation, thermodynamics
├── mini-self-organization-far-equilibrium/  # Self-organization mechanisms, pattern classification
└── mini-turing-pattern-formation/           # Turing instability, 7 RD models, numerical PDE solvers
```

## License

MIT
