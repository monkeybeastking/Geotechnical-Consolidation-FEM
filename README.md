# Finite Element Methods for Geotechnical Consolidation (One Dimensional)

Finite element solvers and verification work for geotechnical consolidation. it includes:
- **Terzaghi 1D (single-layer)**: analytical reference + FEM (verified) + streamlit integration
- **Terzaghi 1D (multi-layer)**: FEM with layered piesewise properties (working; verification + Streamlit integration in progress)
- **Biot consolidation**: planned (coupled displacement–pore pressure)

## Structure
- `scripts/terazaghi_1d/` – analytical + FEM (single-layer)
- `scripts/terazaghi_1d_multilayer/` – FEM (multi-layer)
- `app.py` + `pages/` – Streamlit UI
- `docs/` – notebooks / supporting work
- `.devcontainer/` + `Dockerfile` – reproducible environment (recommended)



