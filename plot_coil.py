import numpy as np 
import matplotlib.pyplot as plt 

odata = np.loadtxt('NUB_150mA_subtracted.txt', delimiter = ' ')
data = odata[:(len(odata)-1)]
normalizer = np.zeros(2)
normalizer[0] = 1#max(data[:,1])
magnitude = np.sqrt(data[:,4]**2 + data[:,5]**2 + data[:,6]**2)
normalizer[1] = 1#max(magnitude)
data[:,:3] = data[:,:3]/normalizer[0]
data[:,4:7] = data[:,4:7]/normalizer[1]


xray = 0.01
zray = 0.01

corner = data[(data[:,0]==xray ) & (data[:,2]==zray)]

comps = ['x', 'y', 'z']

fig = plt.figure()
for i in range(3):
    plt.plot(corner[:,1], corner[:,4+i], label = comps[i])
plt.legend()
plt.show()