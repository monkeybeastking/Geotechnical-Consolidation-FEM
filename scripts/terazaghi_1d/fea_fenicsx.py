from pathlib import Path
from mpi4py import MPI
from petsc4py.PETSc import ScalarType 

import pyvista as pv
import numpy as np
from petsc4py import PETSc
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


def initial_condition1(x, load):
    u = np.full(x.shape[1], load, dtype=np.float64)  
    u[np.isclose(x[0], 0.0)] = 0.0                    
    return u

def initial_condition2(x, load, base):
    z = np.maximum(x[0], 1e-12)                       
    u = (2.0 * load / np.pi) * (
        np.arctan(base / (2.0 * z)) +
        (base * z) / (2.0 * z**2 + 0.5 * base**2)
    )
    u[np.isclose(z, 0.0)] = 0.0                       
    return u




def Get_Terazaghi1D_FEA(H:float, num:int, load:float, Tx:float, time_steps:int, Cv:float, base:float, U0=True):
    dt = Tx / (time_steps -1)

        # interval mesh | MUST keep in positional arguemnts i.e. "comm", "nx"
    msh = mesh.create_interval(
        comm=MPI.COMM_WORLD,
        nx=num,
        points=[0.0, abs(H)],
    )

    # Initial condition callback for fem.Function.interpolate
    if U0:
        initial_condition = lambda x : initial_condition1(x, load)
    else:
        initial_condition = lambda x : initial_condition2(x, load, base)


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
        marker=lambda x: np.isclose(x[0], 0.0)
    )

    dofs = fem.locate_dofs_topological(V, fdim, boundary_facets)
    bc = fem.dirichletbc(PETSc.ScalarType(0), dofs, V)

    uh = fem.Function(V)
    uh.name = "uh"
    uh.interpolate(initial_condition)
    # xdmf.write_function(uh,t)

    # varational form
    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
    # PETSc make sure the background stuff doesnt explode
    c = fem.Constant(msh, Cv)

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

    u_hist = np.zeros((time_steps, uh.x.array.size), dtype=float)
    u_hist[0, :] = uh.x.array.copy()   # initial state
    x = V.mesh.geometry.x[:, 0].copy()


    for i in range(time_steps -1):
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

        u_hist[i + 1, :] = uh.x.array.copy()

    # xdmf.close()
    A.destroy()
    b.destroy()
    solver.destroy()

    u0 = u_hist[0, :]                 # initial condition in space
    local_dcons = 1 - u_hist / u0[None,:]
    local_dcons[:,0] = int(1)        # top drained node




    return local_dcons, u_hist

