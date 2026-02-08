import streamlit as st
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

from scripts.terazaghi_1d_multilayer.mfea_fenics import Get_Terazaghi1dMultilayer_FEA

# setting up Page config for streamlit 
st.set_page_config(layout ="wide")
st.title("1D Terazaghi Consolidation Settlement - Multi Layer")



