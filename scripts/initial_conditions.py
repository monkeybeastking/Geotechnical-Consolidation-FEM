import numpy as np 

def initial_condition(x):
    u = np.full(x.shape[1], load, dtype=np.float64)   # shape (npts,)
    u[np.isclose(x[0], 0.0)] = 0.0                    # enforce u=0 at z=0
    return u  


def Boussinesq_condition(x, load, base):
    z = np.maximum(x, 1e-12)                       # shape (npts,)
    u = (2.0 * load / np.pi) * (
        np.arctan(base / (2.0 * z)) +
        (base * z) / (2.0 * z**2 + 0.5 * base**2)
    )
    u[np.isclose(z, 0.0)] = 0.0                       # optional safety at top
    return u


z = np.linspace(0,5,10)

print(Boussinesq_condition(z, 10, 10))