import streamlit as st
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from scripts.terazaghi_1d.fea_fenicsx import Get_Terazaghi1D_FEA 
from scripts.terazaghi_1d.analytical import Get_Terazaghi1d_Analytical 
from scripts.settlements import get_settlement
from scripts.initial_conditions import Boussinesq_condition

# setting up Page config for streamlit 
st.set_page_config(layout ="wide")
st.title("1D Terazaghi Consolidation Settlement - Singular Layer")
col1, col2 = st.columns([3.5,1])

with col2: 
    H = st.number_input("depth (m)", value=5.0)  # in meters
    num = st.number_input("number of elements", value=10)
    nodes = num + 1
    P = st.number_input("Load applied (kN)", value=100.0) 
    Tx = st.number_input("Final time (days)", value= 365.0)
    Tx = Tx*60*60*24
    time_step = st.number_input("time step", value=1000)
    Cv = st.number_input("Cv (1e-7)", value=2)
    Cv = 1e-7 * Cv 
    Mv = st.number_input("Mv (1e-4) (m^2/kN)", value=5)
    Mv = Mv*1e-4
    N_terms = st.number_input("N terms for analytical solution", value=100)
    initial_conditions = st.toggle("Use uniform initial condition (U0)", value=False) 
    base = st.number_input("base of load placed (m)", value =10.0)




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

# analytical solution plotitng 
analytical_cdata, a_Z, = Get_Terazaghi1d_Analytical(H, Tx, time_step, num, Cv, N_terms)
a_T = np.linspace(0, Tx/(60*60*24), time_step, dtype=float)
analytical_cdata = pd.DataFrame(analytical_cdata, columns= -a_Z, index= a_T)
analytical_settlement = fem_cdata.mean(axis=1)*uniform_total_settlement

error = fem_cdata - analytical_cdata
#RMSE error per time step 
RMSE = np.sqrt((error**2).mean(axis=1).mean(axis=0))


with col1:
    if st.button("1D single layer Solver"): 
        st.subheader("FEM and analytical settlement")
        st.write(f"the total settlement calculare was settlement: {uniform_total_settlement} m")
        spacing = H/num
        initial_conditions = Boussinesq_condition(-Z, P, base)
        A = get_settlement(Mv,spacing, initial_conditions, P)
        st.write(A)
        st.write("RMSE between analytical solution and the fem solution is:",RMSE)
        colsub1, colsub2 = st.columns([2,2])

        with colsub1:
            st.write("settlement from fem (m)", fem_settlement)

        with colsub2:
            st.write("settlement from analytical (m)", analytical_settlement)

        fig, ax = plt.subplots()
        ax.plot(time, -fem_settlement, label="FEM Settlement")
        ax.plot(time, -analytical_settlement, label = "Analytical settlement")
        ax.set_xlabel("Time (days)")
        ax.set_ylabel("Settlement in m")
        ax.legend()
        ax.set_title("settlement over time")
        st.pyplot(fig)