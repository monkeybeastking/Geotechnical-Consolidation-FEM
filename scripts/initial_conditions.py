import numpy as np 


def Boussinesq_condition(x, load, base):
    z = np.maximum(x, 1e-12)                       # shape (npts,)
    u = (2.0 * load / np.pi) * (
        np.arctan(base / (2.0 * z)) +
        (base * z) / (2.0 * z**2 + 0.5 * base**2)
    )
    u[np.isclose(z, 0.0)] = 0.0                       # optional safety at top
    return u

