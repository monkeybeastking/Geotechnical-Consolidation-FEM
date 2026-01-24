from pathlib import Path
from mpi4py import MPI
from petsc4py.PETSc import ScalarType 

import numpy as np
from petsc4py import PETSc
import pyvista as pv
import matplotlib as mpl


import ufl
from dolfinx import fem, io, mesh, plot
from dolfinx.fem import petsc
from dolfinx.fem.petsc import LinearProblem

from dolfinx.fem.petsc import (
    assemble_vector,
    assemble_matrix,
    create_vector,
    apply_lifting,
    set_bc)

t = 0 
T = (60*60*24) * 365# final time (days)
time_steps = 100
H = -5 # depths
nx = 50
dt = T / time_steps  
mv = 0
cv = 0
load = 100 # applied load 

c_1 = 2.0e-3  # m^2/s 


# interval mesh | MUST keep in positional arguemnts i.e. "comm", "nx"
msh = mesh.create_interval(
    comm=MPI.COMM_WORLD,
    nx=nx,
    points=[0.0, abs(H)],
)

initial_condition = lambda x: np.full(x.shape[1], load)

V = fem.functionspace(msh, ("Lagrange", 1))

# Solution functions
u_n = fem.Function(V)
u_n.name = "u_n"
# Initial condition
u_n.interpolate(initial_condition) 
# boundary creation


fdim = msh.topology.dim - 1
boundary_facets = mesh.locate_entities_boundary(
    msh, fdim, 
    marker = lambda x: np.isclose(x[0], H))

dofs = fem.locate_dofs_topological(V=V, entity_dim= fdim, entities= boundary_facets)
bc = fem.dirichletbc( value= PETSc.ScalarType(0), dofs= dofs, V=V)

uh = fem.Function(V)
uh.name = "uh"
uh.interpolate(initial_condition)
# xdmf.write_function(uh,t)

# varational form
u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
 # PETSc make sure the background stuff doesnt explode
c = fem.Constant(msh, c_1)

a = (u * v) * ufl.dx + dt * c * ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
L = (u_n) * v * ufl.dx
bilinear_form = fem.form(a)
linear_form = fem.form(L)

A = assemble_matrix(bilinear_form, bcs = [bc])
A.assemble()
b = petsc.create_vector(fem.extract_function_spaces(linear_form))


# creating linear solver 
solver = PETSc.KSP().create(msh.comm)
solver.setOperators(A)
solver.setType(PETSc.KSP.Type.PREONLY)
solver.getPC().setType(PETSc.PC.Type.LU)



grid = pv.UnstructuredGrid(*plot.vtk_mesh(V.mesh))
plotter = pv.Plotter()
plotter.open_gif("Results/Results.gif", fps=10)

grid.point_data["uh"] = uh.x.array
warped = grid.warp_by_scalar("uh", factor=1)

viridis = mpl.colormaps.get_cmap("viridis").resampled(25)
sargs = dict(
    title_font_size=25,
    label_font_size=20,
    fmt="%.2e",
    color="black",
    position_x=0.1,
    position_y=0.8,
    width=0.8,
    height=0.1,
)

renderer = plotter.add_mesh(
    warped,
    show_edges=True,
    lighting=False,
    cmap=viridis,
    scalar_bar_args=sargs,
    clim=[0, max(uh.x.array)],
)

for i in range(time_steps):
    with b.localForm() as loc_b:
        loc_b.set(0.0)
    assemble_vector(b, linear_form)

    apply_lifting(b ,[bilinear_form], [[bc]])
    b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES,
                  mode=PETSc.ScatterMode.REVERSE)
    set_bc(b, [bc])

    # Solve
    solver.solve(b, uh.x.petsc_vec)

    # Update time-step solution
    u_n.x.array[:] = uh.x.array
    u_n.x.scatter_forward()

        # Update PyVista
    grid.point_data["uh"][:] = uh.x.array
    new_warped = grid.warp_by_scalar("uh", factor=1)
    warped.points[:, :] = new_warped.points
    warped.point_data["uh"][:] = grid.point_data["uh"]

    plotter.write_frame()

plotter.close()
# xdmf.close()
A.destroy()
b.destroy()
solver.destroy


