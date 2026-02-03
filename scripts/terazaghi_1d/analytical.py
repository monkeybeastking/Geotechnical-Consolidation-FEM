import numpy as np
import matplotlib.pyplot as plt

def terzaghi(H, Tx, time_step, nodes, Cv,n_terms):
    Z = np.linspace(0, H, nodes)                 
    Cdata = np.zeros((time_step, nodes), dtype=float)
    T = np.linspace(0, Tx, time_step, dtype=float)

    for i in np.arange(nodes):
        z = Z[i]
        for j in np.arange(len(T)):
            t = T[j]
            S = 0.0
            for n in range(n_terms):
                m = (2 * n) + 1
                S += (1.0 / m) * np.sin((m * np.pi * z) / (2.0 * H)) * np.exp(-((m**2) * (np.pi**2) * Cv * t) / (4.0 * (H**2)))
            
            Cn = (4.0 / np.pi) * S
            Cdata[j, i] = Cn
    return Cdata, Z, T


def Get_Terazaghi1d_Analytical(H, Tx, time_step, num, Cv, n_terms):
    nodes = num + 1
    data, Z, T  = terzaghi(H, Tx, time_step, nodes, Cv,n_terms)
    cdata = 1 - data
    return cdata, Z  