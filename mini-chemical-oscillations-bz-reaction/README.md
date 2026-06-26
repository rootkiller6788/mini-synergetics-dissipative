# mini-chemical-oscillations-bz-reaction

**Belousov-Zhabotinsky (BZ) Reaction -- Chemical Oscillations Far from Equilibrium**

A complete C implementation of the BZ reaction kinetics, Oregonator model,
reaction-diffusion wave propagation, stability analysis, and bistability.

## Module Status: COMPLETE ✅

| Level | Status | Score |
|-------|--------|-------|
| L1: Definitions | COMPLETE | 2 |
| L2: Core Concepts | COMPLETE | 2 |
| L3: Math Structures | COMPLETE | 2 |
| L4: Fundamental Laws | COMPLETE | 2 |
| L5: Algorithms | COMPLETE | 2 |
| L6: Canonical Problems | COMPLETE | 2 |
| L7: Applications | PARTIAL (3/5) | 1 |
| L8: Advanced Topics | PARTIAL (2/5) | 1 |
| L9: Research Frontiers | PARTIAL (documented) | 1 |

**Total Score: 15/18**
**L1-L6 Complete. L7-L9 Partial.**

---

## Quick Start

```
make test         # Build and run all tests
make examples     # Build all examples
make run-examples # Run all examples
make clean        # Clean build artifacts
```

## Code Size

| Directory | Files | Lines |
|-----------|-------|-------|
| include/  | 7 headers | ~956 |
| src/      | 8 C + 1 Lean | ~2621 |
| tests/    | 3 tests | ~518 |
| examples/ | 3 examples | ~366 |
| **Total include/ + src/** | | **3577** |

---

## Core Definitions (L1)

### Oregonator Model (Field-Noyes 1974)

The Oregonator is a 3-variable reduced model of the BZ reaction:

```
eps   * dx/dt = q*y - x*y + x*(1 - x)
eps'  * dy/dt = -q*y - x*y + f*z
        dz/dt = x - z
```

Where:
- **x** = [HBrO2] -- bromous acid concentration (autocatalytic, fast)
- **y** = [Br-] -- bromide ion concentration (control intermediate, medium)
- **z** = [Ce4+] -- cerium(IV) concentration (catalyst redox, slow)

### FKN Mechanism

The Field-Koros-Noyes mechanism involves 3 major processes:

| Process | Reaction | Role |
|---------|----------|------|
| **A** | BrO3- + 2Br- + 3H+ -> 3HOBr | Br- consumption (slow, induction) |
| **B** | BrO3- + HBrO2 + 2Ce3+ + 3H+ -> 2HBrO2 + 2Ce4+ + H2O | Autocatalytic HBrO2 (fast, oscillation engine) |
| **C** | 10Ce4+ + MA + 6H2O -> 10Ce3+ + 2CO2 + 2HCOOH + 10H+ | Ce4+ reduction (slow, Br- regeneration) |

### Key Parameters

| Parameter | Symbol | Typical Value | Role |
|-----------|--------|----------------|------|
| Stoichiometric factor | f | 0.5-3.0 | Br- production; oscillations for f > 0.5 |
| Rate ratio | q | ~8e-5 | BrO3-/Br- rate ratio |
| Time-scale ratio | eps | ~1e-2 | Fast autocatalytic / slow consumption |
| Second time-scale | eps' | ~1e-5 | Extremely slow process C vs fast process B |

---

## Core Theorems (L4)

### 1. Hopf Bifurcation Theorem (Oregonator)

The Oregonator undergoes a supercritical Hopf bifurcation at f = f_H(q)
where a stable limit cycle (oscillation) emerges from the steady state.

Critical value: f_H(q) = 0.5 + 0.5*sqrt(1+4q) for small q.

Implementation: `bz_oregonator_hopf_critical_f()`, `bz_oregonator_hopf_test()`

### 2. Turing Instability Condition

A homogeneous steady state of a 2-species RD system is Turing-unstable iff:

1. fu + gv < 0  (stable without diffusion)
2. fu*gv - fv*gu > 0  (positive determinant)
3. Dv*fu + Du*gv > 0  (cross-diffusion condition)
4. (Dv*fu+Du*gv)^2 > 4*Du*Dv*(fu*gv-fv*gu)  (sufficient)

Implementation: `bz_turing_unstable()`, `bz_turing_critical_wavenumber()`

### 3. Catalyst Conservation Law

Ce3+ + Ce4+ = constant (the Oregonator invariant).

Implementation: `catalyst_conservation_invariant` (Lean formalization)

### 4. Steady-State Existence

For all f > 0, q > 0, the Oregonator has at least one non-negative steady state.
Discriminant: Delta = (q+f-1)^2 + 4q(1+f) >= 0 always.

---

## Core Algorithms (L5)

| Algorithm | Function | Complexity |
|-----------|----------|------------|
| RK4 fixed-step integration | bz_rk4_step() | O(n) per step |
| RK45 adaptive (Dormand-Prince) | bz_rk45_adaptive_step() | O(n), adaptive |
| Cubic eigenvalue solver | solve_cubic() | O(1) analytical |
| DFT Fourier analysis | bz_dft_magnitude() | O(n^2) |
| Fast marching (eikonal) | bz_eikonal_fast_marching() | O(N log N) |
| Saddle-node bisection | bz_find_saddle_node() | O(log(1/eps)) |

---

## Canonical Problems (L6)

1. **Oregonator Oscillation** (example_oregonator_oscillation.c)
   - Integrate Oregonator ODEs with RK4
   - Detect oscillation peaks, measure period and amplitude
   - Compare with analytical period prediction

2. **BZ Spiral Wave** (example_spiral_wave.c)
   - Set up 2D reaction-diffusion domain
   - Seed spiral wave via wave break
   - Simulate with Barkley kinetics
   - Detect spiral tip and count arms

3. **Bistability and Hysteresis** (example_bistability.c)
   - Compute steady states in bistable regime
   - Detect hysteresis loop boundaries
   - Demonstrate chemical switch (flip-flop)

---

## Applications (L7)

| Application | Status |
|-------------|--------|
| Bistable chemical switch | Implemented (bz_chemical_switch) |
| Cardiac arrhythmia modeling | Spiral wave detection |
| Biochemical oscillation | Oregonator as prototype |
| Chemical logic gates | Documented |
| Circadian rhythm modeling | Documented |

---

## Nine-School Curriculum Mapping

| School | Course | Topic |
|--------|--------|-------|
| MIT | 5.60 Thermal Kinetics | BZ mechanism, far-from-equilibrium |
| MIT | 18.385 Nonlinear Dynamics | Hopf bifurcation, limit cycles |
| Stanford | CME 214 Numerical Methods | RK4/45 integration |
| Berkeley | CHEM 221A Adv Kinetics | FKN mechanism |
| CMU | 24-677 Nonlinear Control | Limit cycle analysis |
| Princeton | CHM 515 Chemical Dynamics | Oregonator prototype |
| Caltech | CDS 140 Nonlinear Dynamics | Bifurcation theory |
| Cambridge | Part II Math Biology | Turing patterns, RD |
| ETH | 529-0474 Reaction-Diffusion | BZ waves, pattern formation |

---

## References

- Field, R.J. and Noyes, R.M. (1974) J. Chem. Phys. 60, 1877-1884
- Field, R.J., Koros, E. and Noyes, R.M. (1972) JACS 94, 8649-8664
- Tyson, J.J. (1985) Oscillations and Traveling Waves in Chemical Systems
- Murray, J.D. (2002) Mathematical Biology, Vol I and II
- Strogatz, S.H. (2015) Nonlinear Dynamics and Chaos
- Winfree, A.T. (1987) When Time Breaks Down
- Barkley, D. (1991) Physica D 49, 61-70
- Kuznetsov, Y.A. (2004) Elements of Applied Bifurcation Theory
- Epstein, I.R. and Pojman, J.A. (1998) Nonlinear Chemical Dynamics
