# Course Alignment — Haken Slaving Principle

## Nine-School Curriculum Mapping

### MIT
- **6.241J Dynamic Systems & Control**: Linear stability analysis, eigenvalue methods → `haken_compute_eigenmodes()`, `haken_classify_modes()`
- **16.323 Optimal Control**: Reduced-order modeling → `haken_adiabatic_eliminate()`
- **6.832 Underactuated Robotics**: Center manifold theory → `center_manifold_dimension` (Lean)

### Stanford
- **AA203 Optimal Control**: Model reduction, singular perturbation → `haken_reduced_dynamics()`
- **EE363 Convex Optimization**: Lyapunov theory for stability → `haken_gl_free_energy()`
- **AA274 Multi-agent Systems**: Synchronization, order parameters → `haken_mean_field_op()`

### Berkeley
- **EE221A Linear Systems**: Eigenvalue decomposition, stability → `haken_compute_eigenmodes()`
- **EE222 Nonlinear Systems**: Bifurcation theory, normal forms → `haken_detect_bifurcation()`

### CMU
- **18-771 Linear Systems**: Matrix analysis, spectral theory → `haken_qr_decompose()`
- **24-677 Nonlinear Control**: Sliding mode, bifurcation → `haken_continuation()`
- **24-654 Systems Thinking**: Emergence, self-organization → (order parameter concept)

### Princeton
- **MAE 546 Optimal Control & Estimation**: Reduced models → `haken_adiabatic_eliminate()`
- **ELE 530 Estimation & Identification**: PCA methods → `haken_identify_pca()`

### Caltech
- **CDS110 Introduction to Control**: Stability analysis → `haken_analyze_fixed_point()`
- **CDS140 Nonlinear Dynamics**: Bifurcation, normal forms → `haken_pitchfork_normal_form()`

### Cambridge
- **4F3 Nonlinear & Predictive Control**: Model reduction → `haken_adiabatic_eliminate()`
- **4F2 Robust Control**: Stability margins → `haken_spectral_gap()`

### Oxford
- **B4 Predictive Control**: Reduced-order models → `haken_reduced_dynamics()`
- **C20 Adaptive Control**: Parameter estimation → `haken_fit_gl_parameters()`

### ETH
- **227-0216 System Identification**: PCA, order parameter ID → `haken_identify_pca()`
- **227-0220 Model Reduction**: Balanced truncation, slaving → `haken_adiabatic_eliminate()`

## Textbook Mapping

| Textbook | Chapter | Implementation |
|----------|---------|---------------|
| Haken (1977) Synergetics | §5.1-5.3 Prototype model | `haken_create_prototype_model()` |
| Haken (1977) Synergetics | §7.2-7.4 Slaving principle | `haken_adiabatic_eliminate()` |
| Haken (1983) Advanced Synergetics | §1.3-1.5 Dynamics | `haken_dynamics.c` |
| Cross & Greenside (2009) Pattern Formation | §5 GL equation | `haken_gl_fixed_points()` |
| Kuznetsov (2004) Bifurcation Theory | §5 Normal forms | `haken_pitchfork_normal_form()` |
| Golub & Van Loan (2013) Matrix Computations | §7 QR algorithm | `haken_compute_eigenmodes()` |
| Strogatz (2015) Nonlinear Dynamics | §3 Bifurcations | `haken_detect_bifurcation()` |
| Landau & Lifshitz (1980) Statistical Physics | §146 Phase transitions | `haken_gl_free_energy()` |
