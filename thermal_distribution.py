import numpy as np 


kb = 1.380649e-23

def gen_velocities(N, T, m):
    return np.random.normal(0, np.sqrt(kb*T/m), (round(N),3))