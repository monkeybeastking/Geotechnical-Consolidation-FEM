from pathlib import Path
from mpi4py import MPI
from petsc4py.PETSc import ScalarType 

import numpy as np
from petsc4py import PETSc
import matplotlib.pyplot as plt
import seaborn as sns


import ufl
from dolfinx import fem, io, mesh, plot
from dolfinx.fem import petsc
from dolfinx.fem.petsc import LinearProblem
from dolfinx.mesh import meshtags
from dolfinx.fem.petsc import (
    assemble_vector,
    assemble_matrix,
    create_vector,
    apply_lifting,
    set_bc)

# Embankment properties
Load = 100 # load (kPa)
Base = 10 # embankment width (m)
 
nx = 50

name = ["Made ground", "Soft clay", "Firm clay"]
depths = [1, 2, 4, 5]
H = np.max(depths)
k = [2e-7,2e-7, 2e-7, 2e-7] # this is Cv 
Mv = [5e-4, 10e-4, 5e-4, 5e-4]

T = (60*60*24) * 365 # final time (days)
time_steps = 100
dt = T / time_steps  

H = max(depths)
z = np.linspace(0, H, nx, dtype= np.float64)# kept as small and not zero to prevent error



# these are made specailly for a FEniCSx Function and not to be used normally
def Boussinesq(z):  # no input paramters required as this is only recalled under the same paramter each time
    """
    Boussinesq initial condition (strip load, centreline) for Terzaghi consolidation.
    Assumptions:
    - uniform embankment pressure q
    - embankment width B
    - 1D column at centreline
    - elastic stress distribution used only for initial condition
    """
    z = np.maximum(z,1e-12)
    u = (2.0 * Load / np.pi) * (np.arctan(Base / (2.0 * z)) + (Base * z) / (2.0 * z**2 + 0.5 * Base**2))
    return u

def uniform(z):
    u = np.full(z.shape[1], Load, dtype=np.float64)   # load everywhere
    u[np.isclose(z[0], 0.0)] = 0.0                    # top boundary x=0 -> 0
    return u

def Get_Terazaghi1d_multilayerfem(H:float, num:int, load:float, Tx:float, time_steps:int, Cv:float, base:float, U0=True):