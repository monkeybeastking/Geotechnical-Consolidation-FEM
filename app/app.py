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

p = Get_Terazaghi1D_Numpy(H, num, P, Tx, time_step, Cv)
cdata = 1 - (p / P) 

# used for plotting
Z = -np.linspace(0, H,num+1)
time = np.linspace(0, (Tx/(60*60*24)), time_step)


with col1:
    if st.button("1D local degree of consolidation"): 
        fem_cdata = pd.DataFrame(cdata, columns = Z, index = time).T
        st.subheader("Degree of local consolidation from FEA analysis")
        st.write("X axis = days & Y Axis depths(m) ", fem_cdata)

    if st.button("plot settlement over time"):
        pass 