# Finite Element Methods for Geotechnical Consolidation (Pore Pressure & Settlement)

## Overview
This repository explores numerical methods for geotechnical consolidation, focusing on:
- **Terzaghi 1D consolidation** (diffusion PDE for excess pore pressure)
- **Biot consolidation** (coupled displacement–pore pressure system, planned)

The goal is to build **verified, reproducible implementations** that can be extended to:
layered soils, spatially varying parameters, and eventually coupled consolidation.

---

## What’s implemented
### Terzaghi 1D consolidation (working)
- Finite difference / FEM-style discretisation using **NumPy** (lightweight + fast to iterate)
- **Analytical solution comparison** for verification
- Basic error evaluation (e.g., RMSE / L2-style measures)

> Note: 1D Terzaghi is implemented without FEniCSx on purpose.  
> The core PDE is simple enough to verify cleanly in NumPy first, then scale up to more complex cases.

---

## What’s planned (Roadmap)
### Terzaghi extensions
- [ ] Layered soils (piecewise parameters by depth, e.g., \(c_v(z)\), \(m_v(z)\))
- [ ] More boundary/loading cases (drainage conditions, staged loading)

### Biot consolidation
- [ ] Biot 1D (coupled \(u\)-\(p\) formulation)
- [ ] Biot 2D using **FEniCSx/DOLFINx** (mixed formulation)

---

## Verification philosophy
Numerical results are checked against known solutions and basic sanity checks:
- Analytical vs numerical comparison (Terzaghi 1D)
- Error trends with refinement (mesh/time-step sensitivity where applicable)

(As the project expands, verification scripts/results will be added and standardised.)

---

## FEniCSx / DOLFINx environment
For FEniCSx-based problems (planned), two approaches are included:

- **Conda (recommended):** more reliable and easier to maintain for most users
- **Docker (optional):** useful to avoid local installs, but currently may require fixes due to version mismatch issues

---

## Project structure
```text
.
├── src/
│   ├── terzaghi_1d/      # 1D consolidation (NumPy-based, verified vs analytical)
│   └── biot/             # planned: coupled consolidation
├── .devcontainer/        # optional development container setup
└── README.md