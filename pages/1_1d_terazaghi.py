import streamlit as st
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns

from scripts.terazaghi_1d.fea_fenicsx import Get_Terazaghi1D_FEA 


# this has been put here as it will only use this for this here else where this isnt helpful
def get_settlement(Mv:float, spacing:float, initial_stress:float, Load:float):
    settlement = 0
    for i in initial_stress:
        if i > Load:
            pass # we are not looking at these types of laods
        elif i<= Load:
            settlement_i = Mv * spacing * i
            settlement += settlement_i 
        else:
            break
    return settlement

def initial_condition(x):
    u = np.full(x.shape[1], load, dtype=np.float64)   # shape (npts,)
    u[np.isclose(x[0], 0.0)] = 0.0                    # enforce u=0 at z=0
    return u  

def Boussinesq_condition(x, load, base):
    z = np.maximum(x, 1e-12)                       # shape (npts,)
    u = (2.0 * load / np.pi) * (
        np.arctan(base / (2.0 * z)) +
        (base * z) / (2.0 * z**2 + 0.5 * base**2)
    )
    u[np.isclose(z, 0.0)] = 0.0                       # optional safety at top
    return u





# setting up Page config for streamlit 
st.set_page_config(layout ="wide")
st.title("1D Terazaghi Consolidation Theory")


# pre text within 
st.subheader("Single-layer consolidation (Finite Element Analysis)")
st.write(
    "This dashboard solves 1D Terzaghi consolidation using the finite element method (FEM). "
    " Choose the parameter and numerical settings, then click **Solve** to compute excess pore pressure"
    " dissipation, local degree of consolidation and settlement over time. The initial excess pore pressure can be set to uniform "
    " (constant with depth) or a non-uniform (Boussinesq).")

st.subheader("Assumptions and scope")
st.write(
    "- 1D vertical consolidation only (variables vary with depth *$z$*; no lateral drainage or 3D effects, this is a 1d consolidation settlement model).\n"
    "- Soil is fully saturated and consolidation is governed by Darcy flow and compressibility of soil skeleton.\n"
    "- Single homogeneous layer with constant parameters over depth.\n"
    "- Loading is represented through the initial excess pore pressure profile at $t=0$ (uniform or non-uniform).\n"
    "- Boundary conditions represent drainage at the surface only and a impermeable base.")

st.subheader("Limitations")
st.write(
    "- This is an 1D model: it ignores radial/3D drainage, stress redistribution, and stress-dependent soil properties.\n"
    "- Secondary compression (creep) and immediate and other non-Terzaghi behaviours are not included.\n"
    "- Results aresensitive to numerical choices:\n"
    "  - early times produce **steep pore pressure gradients at and near the drained boundary,\n"
    "  - coarse meshes or large time steps can smooth gradients or create early-time error.\n"
    "- The Boussinesq option here is a 1D non-uniform initial condition, not a full 3D Boussinesq consolidation analysis."
)




col1, col2 = st.columns([3.5,1])
with col2: 
    H = st.number_input("depth (m)", value=5.0)  # in meters
    num = st.number_input("number of elements", value=10)
    nodes = num + 1
    P = st.number_input("Load applied (kN)", value=100.0) 
    Tx = st.number_input("Final time (days)", value= 365.0)
    Tx = Tx*60*60*24
    time_step = st.number_input("time step", value=100)
    Cv = st.number_input("Cv (1e-7)", value=2)
    Cv = 1e-7 * Cv 
    Mv = st.number_input("Mv (1e-4) (m^2/kN)", value=5)
    Mv = Mv*1e-4
    initial_conditions = st.toggle("Use uniform initial condition (U0)", value=False) 
    base = st.number_input("base of load placed (m)", value =2.5)






# loading in function and formatting
uniform_total_settlement = Mv*P*H

# Solving for fem and getting data pre proccesed 
fem_cdata = Get_Terazaghi1D_FEA(H, num, P, Tx, time_step, Cv, base, initial_conditions)
Z = -np.linspace(0, H, num = nodes)
time = np.linspace(0,(Tx/(60*60*24)), num= time_step)
fem_cdata = pd.DataFrame(fem_cdata, columns = Z, index = time)








if initial_conditions == True:
    fem_settlement = fem_cdata.mean(axis=1)*uniform_total_settlement
else:
    spacing = H/num
    initial_conditions = Boussinesq_condition(-Z, P, base)
    total_settlement = get_settlement(Mv,spacing, initial_conditions, P)
    fem_settlement = fem_cdata.mean(axis=1)*total_settlement

spacing = H/num

boussinesq_initial_conditions = Boussinesq_condition(-Z, P, base)
Boussinesq_settlement = get_settlement(Mv,spacing, boussinesq_initial_conditions, P)


uniform_initial_conditions = np.ones((nodes,1))*P

with col1:
    if st.button("Solve"): 

        # plotting initial conditions 
        st.subheader("Initials Conditions for (Boussinesq & Uniform conditions)")
        st.write("Initial excess pore pressure profile at *t* = 0 for uniform and" \
        " Boussinesq loading, as a function of depth.")

        fig_ini, ax_ini = plt.subplots(figsize = (8,5))
        ax_ini.plot(boussinesq_initial_conditions, Z , label="Boussinesq initial conditions")
        ax_ini.plot(uniform_initial_conditions, Z, label = "Uniform initial conditions")
        ax_ini.set_xlabel("Depths (z)")
        ax_ini.set_ylabel("Initial excess pore pressure, (u0 kPa)")
        ax_ini.legend()
        plt.title("Initial Conditions (U0) over depth (z)" )
        st.pyplot(fig_ini)

        # plotting end settlement 
        st.subheader("Settlement response (FEM)")
        st.write("Settlement is computed from the FEM pore pressure and plotted" \
        " over the selected time for given minitial condition.")
        fig_settl, ax = plt.subplots(figsize = (8,5))
        ax.set_ylim(-np.max(uniform_total_settlement),0)
        ax.plot(time, -fem_settlement, label="FEM Settlement")
        ax.set_xlabel("Time (days)")
        ax.set_ylabel("Settlement in m")
        ax.legend()
        ax.set_title("settlement over time")
        st.pyplot(fig_settl)

        st.write(f"Total settlement **(uniform)** : {uniform_total_settlement:.4f} m")
        st.write(f"Total settlement **(Boussinenesq)**: {Boussinesq_settlement:.4f} m")


        # plotting local degree of consolidation heat map
        st.subheader("Local degree of consolidation (depth–time map)")
        st.write("Heatmap of local consolidation response through time and depth, based" \
        " on the normalised excess pore pressure, within the 1D mesh.")
        fig_cons, ax_cons = plt.subplots(figsize = (8,5))
        kx = max(1, len(time)//10)    # ~8 labels across, auto
        ky = max(1, len(Z)//10)  # ~10 labels down, auto 
        ax_cons = sns.heatmap(fem_cdata.T,
                            annot=False,
                            cmap="coolwarm", 
                            xticklabels=time,
                            yticklabels=Z)
        ax_cons.set_xticks(np.arange(0, len(time), kx) + 0.5)
        ax_cons.set_xticklabels([f"{time[i]:.1f}" for i in range(0, len(time), kx)],
                   rotation=0)
        ax_cons.set_yticks(np.arange(0, len(Z), ky) + 0.5)
        ax_cons.set_yticklabels([f"{Z[i]:.1f}" for i in range(0, len(Z), ky)],
                   rotation=0)      
        ax_cons.set_xlabel("Time (Days)")
        ax_cons.set_ylabel("Depth (m)")
        ax_cons.set_title("Local Degree of consolidation within 1D Mesh")
        st.pyplot(fig_cons)