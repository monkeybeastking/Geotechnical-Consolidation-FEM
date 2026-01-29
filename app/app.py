import streamlit as st
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scripts.terazaghi_1d.fea_numpy import Get_Terazaghi1D_Numpy
import seaborn as sns

# preprocessing information 


# setting up Page config for streamlit 
st.set_page_config(layout ="wide")
st.title("Goetechnical Consolidation Settlement Project")
col1, col2 = st.columns([5,1])

with col2: 
    H = st.number_input("depth (m)", value=5.0)  # in meters
    num = st.number_input("number of elements", value=10)
    P = st.number_input("Load applied (kN)", value=100.0) 
    Tx = st.number_input("Final time (days)", value= 365.0)
    Tx = Tx*60*60*24
    time_step = st.number_input("time step", value=10)
    Cv = st.number_input("Cv (1e-7)", value=2)
    Cv = 1e-7 * Cv 
    Mv = st.number_input("Mv (1e-4) (1/kPa or m^2/kN)", value=5)
    Mv = Mv*1e-4



# loading in function and formatting
cdata, total_settlement, Z, time = Get_Terazaghi1D_Numpy(H, num, P, Tx, time_step, Cv, Mv)
# used for plotting
fem_cdata = pd.DataFrame(cdata, columns = Z, index = time)
fem_settlement = fem_cdata.mean(axis=1)*total_settlement

fig, ax = plt.subplots()
ax.plot(time, -fem_settlement, label="FEM Settlement")
ax.set_xlabel("Time (days)")
ax.set_ylabel("Settlement in m")
ax.legend()
ax.set_title("settlement over time")


with col1:
    if st.button("1D single layer Terazaghi Settlement"): 
        st.subheader("Degree of local consolidation from FEA analysis")
        st.write(f"total settlement: {total_settlement} m")
        st.write(fem_cdata.T)
        st.pyplot(fig)