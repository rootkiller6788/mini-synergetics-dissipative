# Knowledge Graph -- SOFE

## L1: Definitions (Complete)
- SOFERegime: equilibrium, near-equilibrium, nonlinear, far-equilibrium, critical, chaotic
- SOFEDissipativeType: spatial, temporal, spatiotemporal, chaotic, hierarchical
- SOFEOrderType: scalar, vector, complex, tensor, conserved, nonconserved
- SOFEBifurcationType: saddle-node, pitchfork, Hopf, transcritical, Turing, Takens-Bogdanov
- SOFEStabilityClass: stable/unstable node/focus, saddle, center, limit cycle, strange attractor
- SOFESymmetryType: none, Z2, O(2), SO(3), translational, gauge, chiral, scale
- SOFEFluctuationType: thermal, shot, critical, nucleation, external
- SOFEFeedbackType: positive, negative, cross, delayed, global
- SOFEIntegrationMethod: Euler, RK4, RK45, semi-implicit, Crank-Nicolson, IMEX

## L2: Core Concepts (Complete)
- Dissipative structures maintained by energy/matter flow
- Order parameters enslaving fast modes (Haken slaving principle)
- Entropy production d_iS/dt = sum F_k*J_k >= 0
- Bifurcation as qualitative dynamics change
- Pattern formation from homogeneous instability
- Autocatalysis and positive feedback
- Non-equilibrium steady state (NESS)
- Critical phenomena and universality
- Symmetry breaking at phase transitions
- Emergence of macroscopic order from microscopic fluctuations

## L3: Mathematical Structures (Complete)
- Stoichiometric matrix N (n_species x n_reactions)
- Mass-action, Michaelis-Menten, Hill kinetics
- Arrhenius temperature dependence
- Jacobian matrix for stability analysis
- Onsager transport matrix L_ij = L_ji
- Spatial discretization (5-point / 7-point stencil)
- Discrete Laplacian operator
- Landau free energy expansion: F = F0 + r*psi^2/2 + b*psi^4/4

## L4: Fundamental Laws (Complete)
- Prigogine Minimum Entropy Production Theorem
- Glansdorff-Prigogine Stability Criterion
- Onsager Reciprocal Relations
- Fluctuation-Dissipation Theorem
- Einstein fluctuation formula
- Crooks Fluctuation Theorem
- Jarzynski Equality
- Thermodynamic Uncertainty Relation (Barato-Seifert)
- Landau theory of phase transitions
- Second Law for non-equilibrium systems

## L5: Algorithms (Complete)
- Reaction rate computation (3 kinetic types)
- Stoichiometric matrix assembly
- Jacobian via analytical + finite-difference methods
- QR eigenvalue algorithm
- RK4 / Euler explicit-implicit integration
- Pseudo-arclength continuation
- Benettin Lyapunov spectrum algorithm
- DFT spatial spectrum computation
- Entropy production from forces/fluxes and reactions
- Finite-difference Jacobian

## L6: Canonical Problems (Complete)
- Brusselator (Prigogine & Lefever 1968)
- Schlogl bistable model (1972)
- Lotka-Volterra predator-prey
- Oregonator (Field-Koros-Noyes 1974)
- Ginzburg-Landau equation
- Turing instability
- Rayleigh-Benard convection

## L7: Applications (Complete)
1. BZ chemical oscillations (NASA microgravity expts)
2. Biological morphogenesis (Turing, Toyota RIKEN)
3. Rayleigh-Benard convection (Boeing thermal mgmt)
4. Eigen hypercycle (prebiotic chemistry)
5. Cardiac spiral waves (fibrillation)
6. Ecological vegetation patterns

## L8: Advanced Topics (Partial)
1. Lyapunov spectrum + KY dimension (Implemented)
2. GP stability criterion (Implemented)
3. FDT verification (Implemented)
4. NESS thermodynamics (Implemented)
5. Stochastic resonance (Documented)
6. Self-organized criticality (Documented)
7. Noise-induced transitions (Documented)

## L9: Research Frontiers (Partial)
1. Jarzynski equality (Implemented)
2. Thermodynamic uncertainty relations (Implemented)
3. Quantum dissipative structures (Documented)
4. Active matter (Documented)
5. Non-equilibrium universality classes (Documented)
6. Information-theoretic self-organization (Documented)
