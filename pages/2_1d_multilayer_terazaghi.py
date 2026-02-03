import streamlit as st
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from scripts.terazaghi_1d.fea_fenicsx import Get_Terazaghi1D_FEA 
from scripts.settlements import get_settlement
from scripts.initial_conditions import Boussinesq_condition

# setting up Page config for streamlit 
st.set_page_config(layout ="wide")
st.title("1D Terazaghi Consolidation Settlement - Multi Layer")
