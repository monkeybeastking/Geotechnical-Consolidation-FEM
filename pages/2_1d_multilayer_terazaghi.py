import streamlit as st
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

from scripts.terazaghi_1d_multilayer.mfea_fenics import Get_Terazaghi1dMultilayer_FEA

# setting up Page config for streamlit 
st.set_page_config(layout ="wide")
st.title("1D Terazaghi Consolidation Settlement - Multi Layer")

st.write(
    "This dashboard solves 1D Terzaghi consolidation for a **multi-layer soil profile** using the finite element method (FEM). "
    "Each layer can have different consolidation parameters. After selecting the soil profile and numerical settings, "
    "click **Solve** to compute excess pore pressure dissipation and settlement over time. The initial excess pore pressure can be set to uniform "
    " (constant with depth) or a non-uniform (Boussinesq).")

st.write("Please ensure that when entering Depths, Mv, Cv information," \
        "that these are comma seperated")



col1, col2 = st.columns([3.5,1])
with col2: 
    num = st.number_input("number of elements", value=100)

    Load = st.number_input("Load applied (kPa)", value=100.0) 

    Tx = st.number_input("Final time (days)", value= 365.0)
    Tx = Tx*60*60*24

    time_step = st.number_input("time step", value=1000)

    initial_conditions = st.toggle("Use uniform initial condition (U0)", value=False) 
    base = st.number_input("base of load placed (m)", value =2.5)

    depths_text = st.text_input("Depths (comma separated)", "1, 2, 4, 5")
    depths = [float(x.strip()) for x in depths_text.split(",") if x.strip()]
    H = np.max(depths)

    Mv_text = st.text_input("Mv (comma separated)", "5e-4, 10e-4, 5e-4, 5e-4")
    Mv = [float(x.strip()) for x in Mv_text.split(",") if x.strip()]

    Cv_text = st.text_input("Cv (comma separated)", "2e-7,2e-7, 2e-7, 2e-7")
    Cv = [float(x.strip()) for x in Cv_text.split(",") if x.strip()]



with col1:
    if st.button("Solve"): 
        # kept this inside the button so it doesnt try to keep rerunning
        fem_cdata, fem_udata, settlement = Get_Terazaghi1dMultilayer_FEA(depths, num, Load, Tx, time_step, Cv ,Mv, base, U0=initial_conditions)
        Z = -np.linspace(0, H, num + 1 )
        time = np.linspace(0,(Tx/(60*60*24)), num = time_step)

        fem_setdata = np.sum((settlement[:,None] * fem_cdata.T), axis = 0)

        
        # plotting initial conditions 
        st.subheader("Initials Conditions")
        st.write("Initial excess pore pressure profile at *t* = 0.")

        fig_ini, ax_ini = plt.subplots(figsize = (8,5))
        ax_ini.plot(fem_udata.T[:,0], Z , label="initial conditions")
        ax_ini.set_ylabel("Depths (z)")
        ax_ini.set_xlabel("Initial excess pore pressure, (u0 kPa)")
        ax_ini.legend()
        plt.title("Initial Conditions (U0) over depth (z)" )
        st.pyplot(fig_ini)


        # plotting end settlement 
        st.subheader("Settlement response (FEM)")
        st.write("Settlement is computed from the FEM pore pressure and plotted" \
        " over the selected time for given initial condition.")

        st.write(f"Total Consolidation Settlement:  {np.sum(settlement, axis =0):.4f}m")
        st.write(f"Total ConsolidationSettlement after {Tx/(60*60*24)} days:    {np.max(fem_setdata):.4f}m")
        fig_settl, ax = plt.subplots(figsize = (8,5))
        ax.set_ylim(-np.max(np.sum(settlement, axis = 0)),0)
        ax.plot(time, -fem_setdata, label="FEM Settlement")
        ax.set_xlabel("Time (days)")
        ax.set_ylabel("Settlement in m")
        ax.legend()
        ax.set_title("settlement over time")
        st.pyplot(fig_settl)



        # plotting excess pore pressure heat map
        st.subheader("Excess Pore pressure dissipation (depth time map)")
        st.write("Heatmap of Pore pressure dissipation through time and depth.")
        fig_cons, ax_cons = plt.subplots(figsize = (8,5))
        kx = max(1, len(time)//10)    # ~8 labels across, auto
        ky = max(1, len(Z)//10)  # ~10 labels down, auto 
        ax_cons = sns.heatmap(fem_udata.T,
                            annot=False,
                            cmap="Blues", 
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
        ax_cons.set_title("Excess Pore Pressure dissipation over time in a 1D Mesh")
        st.pyplot(fig_cons)