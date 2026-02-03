import streamlit as st
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from scripts.terazaghi_1d.fea_fenicsx import Get_Terazaghi1D_FEA 
from scripts.settlements import get_settlement
from scripts.initial_conditions import Boussinesq_condition

# setting up Page config for streamlit 
st.set_page_config(layout ="wide")
st.title("1D Terazaghi Consolidation Settlement - Singular Layer")
col1, col2 = st.columns([3.5,1])

with col2: 
    H = st.number_input("depth (m)", value=5.0)  # in meters
    num = st.number_input("number of elements", value=100)
    nodes = num + 1
    P = st.number_input("Load applied (kN)", value=100.0) 
    Tx = st.number_input("Final time (days)", value= 365.0)
    Tx = Tx*60*60*24
    time_step = st.number_input("time step", value=10000)
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


with col1:
    if st.button("Solve"): 
        st.subheader("FEM and analytical settlement")
        st.write(f"Total settlement (uniform): {uniform_total_settlement:.4f} m")
        spacing = H/num
        initial_conditions = Boussinesq_condition(-Z, P, base)
        A = get_settlement(Mv,spacing, initial_conditions, P)
        st.write(f"Total settlement (Boussinenesq):{A:.4f} m")
        colsub1, colsub2 = st.columns([2,2])

        st.write("settlement from fem (m)", fem_settlement)

        fig, ax = plt.subplots()
        ax.set_ylim(-np.max(uniform_total_settlement),0)
        ax.plot(time, -fem_settlement, label="FEM Settlement")
        ax.set_xlabel("Time (days)")
        ax.set_ylabel("Settlement in m")
        ax.legend()
        ax.set_title("settlement over time")
        st.pyplot(fig)