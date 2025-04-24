import os 
import numpy as np 
import matplotlib.pyplot as plt 

dname = 'STC_3DMAPS_2025/response_curves'
coil = 'C'

fig = plt.figure()

for file in os.listdir(dname):
    coilname = file.split('_')[0]
    if coilname[-1] == coil:
        data = np.loadtxt(dname+'/'+file, delimiter = ' ')
        axis = data[(data[:,0] == 0) & (data[:,2] == 0)]
        plt.scatter(axis[:,1], axis[:,4])


#for i in range(3):
#    plt.scatter(axis[:,1], axis[:,3+i])

plt.show()