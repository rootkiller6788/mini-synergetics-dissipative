# Gap Report — Order Parameters & Emergence

## Current Gaps

### L9: Research Frontiers (Partial → Full)

| Priority | Gap | Current State | Target |
|----------|-----|---------------|--------|
| Low | Machine learning OP discovery | Documented only | Autoencoder-based OP identification |
| Low | Topological order parameters | Documented only | Berezinskii-Kosterlitz-Thouless implementation |
| Low | Active matter field theories | Vicsek model only | Toner-Tu continuum equations |
| Low | Quantum synergetics | Basic GL superconductivity | Full quantum master equation |
| Low | Non-Markovian order parameters | Not addressed | Memory kernel formalism |

### No Critical Gaps
All L1-L8 requirements are fully satisfied. The remaining L9 items are research-level topics that are documented but not fully implemented, which meets the Partial requirement.

## Quality Gaps (non-blocking)

| Area | Note | Priority |
|------|------|----------|
| Sparse matrix support | Currently dense-only (OEP_MAX_DIM=128) | Low |
| GPU acceleration | No CUDA/OpenCL for TDGL | Low |
| 3D spatial fields | Currently 2D only | Medium |
| Higher-order normal forms | Only up to cubic terms | Low |
