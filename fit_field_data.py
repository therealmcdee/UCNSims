import numpy as np 
import matplotlib.pyplot as plt 
from scipy.optimize import minimize

from simple_field import generic_xderiv_3, generic_zderiv_3

filename = 'NUB_150mA_subtracted.txt'

def residual(x, positions, true):
    predicted = np.zeros((len(positions), 2))
    predicted[:,0] = generic_xderiv_3(x, positions)
    predicted[:,1] = generic_zderiv_3(x, positions)
    return np.sqrt(sum((predicted[:,0]-true[:,0])**2) + sum((predicted[:,1]-true[:,2])**2))

odata = np.loadtxt(filename, delimiter = ' ')
data = odata[:(len(odata)-1)]
normalizer = np.zeros(2)
normalizer[0] = 1#max(data[:,1])
magnitude = np.sqrt(data[:,4]**2 + data[:,5]**2 + data[:,6]**2)
normalizer[1] = 1#max(magnitude)
data[:,:3] = data[:,:3]/normalizer[0]
data[:,4:7] = data[:,4:7]/normalizer[1]



onaxis = data[(data[:,0]==0) & (data[:,2]==0)]
coefs = np.zeros((len(onaxis),9))
plane_diffs = np.zeros(len(onaxis))
for i in range(len(onaxis)):
    field_data = data[(data[:,1]==onaxis[i,1])]
    res = minimize(residual, x0 = np.zeros(9), args=(field_data[:,:3], field_data[:, 4:7]))
    coefs[i] = res.x
    plane_diffs[i] = residual(coefs[i], field_data[:,:3], field_data[:,4:7])

plane = 10
plane_dat = data[(data[:,1]==onaxis[i,1])]
n = 100
x = np.linspace(-max(data[:,0]), max(data[:,0]), n)
mesh = np.zeros((len(x)**2, 3))
cnt = 0
for i in range(len(x)):
    for j in range(len(x)):
        mesh[cnt][0] = x[i]
        mesh[cnt][1] = plane
        mesh[cnt][2] = x[j]
        cnt += 1
print(coefs[plane])

order = [0, 0, 0, 1, 1, 2, 2, 3, 3]
dcoefs = np.zeros(np.shape(coefs))
for i in range(len(coefs[0])):
    dcoefs[:, i] = coefs[:,i]/(normalizer[0]**order[i])
print(dcoefs[plane])

fig0 = plt.figure()
ax = fig0.add_subplot(projection = '3d')
ct = ax.scatter(1e2*data[:,0]*normalizer[0], 1e2*data[:,1]*normalizer[0], 1e2*data[:,2]*normalizer[0], c = data[:,6]*1e6*normalizer[1])
plt.colorbar(ct, label = r'$B_{z}$ $[\mu T]$')
ax.set_xlabel('x [cm]')
ax.set_ylabel('y [cm]')
ax.set_zlabel('z [cm]')
plt.title(filename)
plt.tight_layout()

fig = plt.figure()

for i in range(len(coefs[0])):
    plt.plot(onaxis[:,1], coefs[:, i], label = r'$U_{}$'.format(order[i]))
plt.legend()
plt.tight_layout()

bzfield = generic_zderiv_3(dcoefs[plane], mesh)
bxfield = generic_xderiv_3(dcoefs[plane], mesh)

fig2, ax2 = plt.subplots(1,2)
c0 = ax2[0].scatter(plane_dat[:,0]*normalizer[0], plane_dat[:,2]*normalizer[0], c = plane_dat[:,6]*normalizer[1]*1e6, vmin = min(bzfield*normalizer[1]*1e6), vmax = max(bzfield*normalizer[1]*1e6))
fig2.colorbar(c0)
#ax2[0].set_title(r'$B_{x}$ $[\mu T]$')
c1 = ax2[1].scatter(mesh[:,0]*normalizer[0], mesh[:,2]*normalizer[0], c = bzfield*normalizer[1]*1e6)
fig2.colorbar(c1)
ax2[1].set_title(r'$B_{z}$ $[\mu T]$')
plt.tight_layout()

fig3, ax3 = plt.subplots(1,2)
c0 = ax3[0].scatter(plane_dat[:,0]*normalizer[0], plane_dat[:,2]*normalizer[0], c = plane_dat[:,4]*normalizer[1]*1e6, vmin = min(bxfield*normalizer[1]*1e6), vmax = max(bxfield*normalizer[1]*1e6))
fig3.colorbar(c0)
#ax2[0].set_title(r'$B_{x}$ $[\mu T]$')
c1 = ax3[1].scatter(mesh[:,0]*normalizer[0], mesh[:,2]*normalizer[0], c = bxfield*normalizer[1]*1e6)
fig3.colorbar(c1)
ax3[1].set_title(r'$B_{x}$ $[\mu T]$')
plt.tight_layout()


positions = data[:,:3]
bzerror = np.zeros(len(onaxis))
bxerror = np.zeros(len(onaxis))
for i in range(len(onaxis)):
    plane = positions[(positions[:,1] == onaxis[i,1])]
    bzrecreate = generic_zderiv_3(dcoefs[i], plane)
    bzerror[i] = np.sqrt(sum((bzrecreate-data[(data[:,1]==onaxis[i,1])][:,6])**2))
    bxrecreate = generic_xderiv_3(dcoefs[i], plane)
    bxerror[i] = np.sqrt(sum((bxrecreate-data[(data[:,1]==onaxis[i,1])][:,4])**2))
    

fig4 = plt.figure()
plt.plot(onaxis[:,1]*normalizer[0], plane_diffs*normalizer[1]*1e6, label = 'Minimizing Function')
plt.plot(onaxis[:,1]*normalizer[0], bzerror*normalizer[1]*1e6, label = r'$B_{z}$')
plt.plot(onaxis[:,1]*normalizer[0], bxerror*normalizer[1]*1e6, label = r'$B_{x}$')
plt.xlabel('y [m]')
plt.title(r'$\sqrt{\left(B_{x}-B_{x, true}\right)^{2} + \left(B_{z}-B_{z,true}\right)^{2}}$ $[\mu T]$')
plt.legend()
plt.tight_layout()
plt.show()