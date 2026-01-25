import numpy as np
import matplotlib.pyplot as plt


def create_local(Cv, Z):
    he = np.diff(Z)
    if np.round(he[0:1]) == np.round(he[1:2]):
        he = he[0:1]
        me = np.array([[2,1], [1,2]])
        Me = (he / 6) * me
        ke = np.array([[1,-1], [-1,1]])
        Ke = (Cv / he) * ke
        return Me, Ke 
    else:
        print("ERROR - 'he' not uniform")

def create_global(Matrix, num_nodes):
    num_nodes = int(num_nodes)
    Global_M = np.zeros((num_nodes,num_nodes))
    n = np.arange(num_nodes)
    for i in range(np.max(n)):
        i = int(i)
        Global_M[i : (i+2), i : (i+2)] += Matrix
    return Global_M

def boundary_conditions(LF, nodes):
    p = 0
    LF[:,p] = 0
    LF[p,:] = 0
    LF[p,p] = 1
    return LF

def solve_timestep(GMe, B, U0, time_step, nodes):
    p = 0
    GMe[:,p] = 0
    GMe[p,:] = 0
    GMe[p,p] = 1

    U0[p] = 0
    linear_form = GMe @ U0
    loop = np.arange(int(time_step))

    for i in loop:
        if i == 0: 
            pass
        elif i == 1:
            data = U0
            linear_form = GMe @ data
            linear_form[p] = 0
            Un = np.linalg.solve(B, linear_form)
            Un = np.array([Un])
            data = data.reshape(-1, 1).T
            data = np.concatenate((data,Un), axis=0)
        elif i > 1:
            Un = data[(i-1):i,:nodes]
            Un = Un.T
            linear_form = GMe @ Un 
            Un1 = np.linalg.solve(B, linear_form)
            Un1 = Un1.reshape(-1, 1).T
            data = np.concatenate((data,Un1), axis=0)
    return data
          

def Get_Terazaghi1D_Numpy(H, num, P, Tx, time_step, Cv):
    nodes = num + 1
    dt = Tx / time_step

    Z = -np.linspace(0, H, num = nodes)
    U0 = P*np.ones(Z.shape, dtype=float)

    Me, Ke = create_local(Cv, Z)
    GMe, GKe = create_global(Me, nodes), create_global(Ke, nodes)
    bilinear_form = (GMe) + (dt * GKe)

    B = boundary_conditions(bilinear_form, nodes)
    data = solve_timestep(GMe, B, U0, time_step, nodes)

    return data


if __name__ == "__main__":
    H = 5
    num = 100
    P = 100
    Tx = 60*60*24*150
    time_step = 100
    Cv = 2e-7

    p = Get_Terazaghi1D_Numpy(H, num, P, Tx, time_step, Cv)
    print(p)
