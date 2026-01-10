# Finite Element Analysis methods in solving geotechnical displacement / pore pressure

## Project Purpose and Description
Fundamentally, the two governing equations related to consolidation settlement are the **Terzaghi** partial differential equation (heat/diffusion equation) and the **Biot consolidation theory** PDE system.

- **Terzaghi theory** only considers vertical drainage. As a result, extending from 1D Terzaghi to 2D Terzaghi does not necessarily add meaningful complexity or purpose, other than accounting for more tedious loading and stress distributions in soil.
- **Biot consolidation theory** is a coupled PDE system. This has not yet been solved within this project.

**Note:** The current 1D Terzaghi consolidation implementation does **not** use FEniCSx. Since the 1D Terazaghi partial differentrial equation is a simple diffusion differential equation by nature, an implemention using **only** using numpy has been undertaken. For more complex, or multi level / layer geology FEniCS have been incorporated. FEniCSx will be progressively introduced in more advanced cases, such as coupled Biot consolidation or multi-layer Terzaghi consolidation.

## FEniCS / DOLFINx
When using FEniCS, you can work either with Docker or Conda.

- **Docker**: Docker files are provided to allow use of this project with VS Code or other IDEs while avoiding local installation issues with the FEniCS library. However, several errors and version mismatches currently exist in these files, so this method is not recommended at present.
- **Conda (preferred method)**: It is recommended to use Conda. The official FEniCS workshops provide a Conda environment file, which is more reliable and easier to maintain.

## src/
- **1D Terzaghi consolidation analysis**: Implemented using standard Python tools (NumPy/Pandas), including comparison between the FEM solution and the analytical solution with error evaluation.
- **Biot consolidation settlement**: Not yet implemented.

### Project structure
```text
.
├── src/
│   ├── terzaghi_1d/
│   └── biot/        (planned)
├── .devcontainer/  (optional)
└── README.md
