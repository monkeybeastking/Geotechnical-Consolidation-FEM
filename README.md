# Finite Element Analysis methods in solving geotechnical settlement / pore pressure

## Project Purpose and Description
Fundamentally, the two governing equations related to consolidation settlement are the **Terzaghi** partial differential equation (heat/diffusion equation) and the **Biot consolidation theory** PDE system.

- **Terzaghi theory** only considers vertical drainage. As a result, extending from 1D Terzaghi to 2D Terzaghi does not necessarily add meaningful complexity or purpose, other than accounting for more tedious loading and stress distributions in soil.
- **Biot consolidation theory** is a coupled PDE system. This has not yet been solved within this project.

## FEniCS / DOLFINx
When using FEniCS, you can work either with Docker or Conda.

- **Docker**: Docker files are provided to allow use of this project with VS Code or other IDEs while avoiding local installation issues with the FEniCS library. However, several errors and version mismatches currently exist in these files, so this method is not recommended at present.
- **Conda (preferred method)**: It is recommended to use Conda. The official FEniCS workshops provide a Conda environment file, which is more reliable and easier to maintain.

## src/
- **1D Terzaghi consolidation analysis**: Includes comparison between the FEM solution and the analytical solution, along with error evaluation.
- **Biot consolidation settlement**: Not yet implemented.

## References
- FEniCSx tutorial: https://jsdokken.com/dolfinx-tutorial/chapter2/heat_equation.html  
- ParaView: https://www.paraview.org/download/  
- NumPy documentation: https://numpy.org/doc/stable/user/absolute_beginners.html
