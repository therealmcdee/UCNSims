import os 
import numpy as np 
import matplotlib.pyplot as plt 

fname = 'STC_3DMAPS_2025/response_curves/SUD_response.txt'

data = np.loadtxt(fname, delimiter = ' ')

axis = data[(data[:,0] == 0) & (data[:,2] == 0)]

fig = plt.figure()
for i in range(3):
    plt.scatter(axis[:,1], axis[:,3+i])

plt.show()