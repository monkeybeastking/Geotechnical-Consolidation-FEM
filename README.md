# Finite Element Analysis methods in solving for geotechnical settlement / pore presusre 


Project Purpose and Descriptions
Fundermentally the two governing equation related to consolidation settlement is namely Terazaghi PDE (Heat or Diffusion Equastion) or Biot Consolidation Thoery PDE. 
- Terazaghi theory only considers virtical drainage and thus moving from 1d terazaghi to 2d terazaghi doesnt necessarily and more complexity nor purpose other than loads and stress in soils can become mnore tedious.
- Biot Consolidation theory is a coupled PDE system. This has not as of yet been solved by myself. 


FEniCS / Dolfinx / 
in using FEniCS, you can either use docker or Conda.
- Docker files provided here allow you to use this program and either vscode or whatever other IDE to work with the FEniCS library without causing issue. However several errors and incorrect version are not currently are within this file so it is prefered to use conda (look at FEniCS workshop they provide a conda file in which is usefull).
- Conda (prefeered method) 

.src:
- 1D Terazaghi consolidation analysis. This also compared the FEM model with the analytical solution to get the error.
- Biot Consolidation settlement (yet to be done)



References:
FEnicsX: https://jsdokken.com/dolfinx-tutorial/chapter2/heat_equation.html
Paraview: https://www.paraview.org/download/
numpy: https://numpy.org/doc/stable/user/absolute_beginners.html 
