# Finite Element Methods for Geotechnical Consolidation (One Dimensional)

This repository implements FEniCSx to model excess pore pressure dissipation in soils over time, due to consolidation. The project focuses on finite element solvers, verification, and streamlit integration of models.

The repository currently includes:
- **Terzaghi 1D Consolidation (Single Layer)**: Analytical reference solution + verified FEM implementation + Streamlit integration

- **Terzaghi 1D Consolidation (Multi-Layer)**: FEM model with layered, piecewise material properties *(working; verification in progress)* + Streamlit integration

- **Terzaghi 2D Consolidation (Single or Multi-Layer)**: Extension to 2D mesh-based FEM modelling *(under active development)*

- **Biot Consolidation (Planned)**: Future implementation of fully coupled displacement pore pressure consolidation theory

## Repository Structure
- `app.py` + `pages/` – Streamlit user interface

- `docs/` – Notebooks and supporting derivations
 
- `.devcontainer/` + `Dockerfile` – Reproducible development environment  

### Scripts
- `scripts/terzaghi_1d/` – Analytical + FEM solver (single-layer)
  
- `scripts/terzaghi_1d_multilayer/` – FEM solver (multi-layer)
   
- `scripts/terzaghi_2d/` – 2D FEM consolidation *(under development)*  


## References

- Terzaghi, K. (1943). *Theoretical Soil Mechanics*. Wiley.  

- Biot, M. A. (1941). General theory of three-dimensional consolidation.  
  *Journal of Applied Physics*, 12(2), 155–164.  

- FEniCSx Project Documentation:  
  https://docs.fenicsproject.org/

- Larson, M. G., & Bengzon, F. (2013). *The Finite Element Method: Theory, Implementation, and Applications*. Springer.  
  *(Used as supporting background for FEM formulation.)*




