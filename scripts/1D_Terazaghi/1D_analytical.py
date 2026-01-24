import numpy as np
import matplotlib.pyplot as plt

H = 5
n = 100
nodes = n + 1
P = 100
Tx = 60*60*24*150
time_step = 100
dt = Tx / time_step
Cv = 2e-7

time_factor = (Cv * dt) / (H**2) 

def terzaghi(H, Tx, time_step, nodes, Cv):
    Z = np.linspace(0, H, nodes)                 
    Cdata = np.zeros((time_step, nodes), dtype=float)
    T = np.linspace(0, Tx, time_step, dtype=float)
    N_terms = 1000

    for i in np.arange(nodes):
        z = Z[i]
        for j in np.arange(len(T)):
            t = T[j]
            S = 0.0
            for n in range(N_terms):
                m = (2 * n) + 1
                S += (1.0 / m) * np.sin((m * np.pi * z) / (2.0 * H)) * np.exp(-((m**2) * (np.pi**2) * Cv * t) / (4.0 * (H**2)))
            
            Cn = (4.0 / np.pi) * S
            Cdata[j, i] = Cn
    return Cdata

cdata = terzaghi(H, Tx, time_step, nodes, Cv)
cdata = 1 - cdata
print(cdata)
Z = -np.linspace(0, H, nodes)
