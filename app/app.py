import streamlit as st
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scripts.terazaghi_1d.fea_numpy import Get_Terazaghi1D_Numpy


st.title("Goetechnical Consolidation Settlement Project")




H = 5
num = 100
P = 100
Tx = 60*60*24*150
time_step = 100
Cv = 2e-7

p = Get_Terazaghi1D_Numpy(H, num, P, Tx, time_step, Cv)

cdata = 1 - (p / P) 
Z = -np.linspace(0, H,num+1)
fem_cdata = pd.DataFrame(cdata, columns = Z, index = np.arange(0, time_step))


st.write("fea data", fem_cdata.head())

