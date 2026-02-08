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

## Run
```bash
# Streamlit
streamlit run app.py

# Scripts (from repo root)
python scripts/terazaghi_1d/analytical.py
python scripts/terazaghi_1d/fea_fenicsx.py
python scripts/terazaghi_1d_multilayer/mfea_fenics.py

Implemented:
- FEM in **FEniCSx**
- Layering handled via **cell markers** and a DG0 field for material parameters (e.g. \(C_v\))

Status:
- **Runs and produces outputs**
- **Verification and Streamlit integration are in progress** (partially complete)

---

## Repository structure


