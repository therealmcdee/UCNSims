import numpy as np 
import matplotlib.pyplot as plt 
from scipy.optimize import minimize

from simple_field import generic_xderiv_3, generic_zderiv_3

filename = 'STC_3DMAPS_2025/Acoils/NUA_10mA_subtracted.txt'


def chisqbx(x, positions, true):
    predicted = generic_xderiv_3(x, positions)
    return sum(((predicted - true[:,0])**2/abs(true[:,0])))

def chisqbz(x, positions, true):
    predicted = generic_zderiv_3(x, positions)
    return sum(((predicted - true[:,2])**2/abs(true[:,2])))

def minfunc1(x, positions, true):               
    return np.sqrt(chisqbx(x, positions, true)**2 + chisqbz(x, positions, true)**2)


def minfunc2(x, positions, true):
    return chisqbx(x, positions, true) + chisqbz(x, positions, true) 

odata = np.loadtxt(filename, delimiter = ' ')
data = odata[:(len(odata)-1)]
normalizer = np.zeros(2)
normalizer[0] = 1#max(data[:,1])
magnitude = np.sqrt(data[:,4]**2 + data[:,5]**2 + data[:,6]**2)
normalizer[1] = 1#max(magnitude)
data[:,:3] = data[:,:3]/normalizer[0]
data[:,4:7] = data[:,4:7]/normalizer[1]



onaxis = data[(data[:,0]==0) & (data[:,2]==0)]
coefs1 = np.zeros((len(onaxis),9))
coefs2 = np.zeros((len(onaxis),9))
xcoefs = np.zeros((len(onaxis),9))
zcoefs = np.zeros((len(onaxis),9))
method1err = np.zeros(len(onaxis))
method2err = np.zeros(len(onaxis))
xmethoderr = np.zeros(len(onaxis))
zmethoderr = np.zeros(len(onaxis))
for i in range(len(onaxis)):
    field_data = data[(data[:,1]==onaxis[i,1])]
    res1 = minimize(minfunc2, x0 = np.zeros(9), args=(field_data[:,:3], field_data[:, 4:7]))
    coefs1[i] = res1.x
    method1err[i] = minfunc2(coefs1[i], field_data[:,:3], field_data[:,4:7])
    res2 = minimize(minfunc1, x0 = np.zeros(9), args=(field_data[:,:3], field_data[:, 4:7]))
    coefs2[i] = res2.x
    method2err[i] = minfunc1(coefs2[i], field_data[:,:3], field_data[:, 4:7])
    xres = minimize(chisqbx, x0 = np.zeros(9), args=(field_data[:,:3], field_data[:, 4:7]))
    xcoefs[i] = xres.x
    xmethoderr[i] = chisqbx(xcoefs[i], field_data[:,:3], field_data[:, 4:7])
    zres = minimize(chisqbz, x0 = np.zeros(9), args=(field_data[:,:3], field_data[:, 4:7]))
    zcoefs[i] = zres.x
    zmethoderr[i] = chisqbz(zcoefs[i], field_data[:,:3], field_data[:, 4:7])

plane = 4
plane_dat = data[(data[:,1]==onaxis[plane,1])]
n = 100
x = np.linspace(-max(data[:,0]), max(data[:,0]), n)
mesh = np.zeros((len(x)**2, 3))
cnt = 0
for i in range(len(x)):
    for j in range(len(x)):
        mesh[cnt][0] = x[i]
        mesh[cnt][1] = onaxis[plane,1]
        mesh[cnt][2] = x[j]
        cnt += 1


order = [0, 0, 0, 1, 1, 2, 2, 3, 3]
dcoefs = np.zeros(np.shape(coefs1))
for i in range(len(coefs1[0])):
    dcoefs[:, i] = coefs2[:,i]/(normalizer[0]**order[i])
print(coefs2[plane])


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
guide_r = 0.035
for i in range(len(coefs1[0])-1):
    plt.plot(onaxis[:,1], 1e6*dcoefs[:, i+1]*(guide_r**order[i+1]), label = r'$U_{}$'.format(order[i+1]))
plt.title(f'Moment Contribution at r = {guide_r} m')
plt.ylabel(r'$U_{m}r^{m}$ $[\mu T]$')
plt.xlabel('y[m]')
plt.legend()
plt.tight_layout()

bzfield = generic_zderiv_3(dcoefs[plane], mesh)
bzcompp = generic_zderiv_3(dcoefs[plane], plane_dat[:,:3])
bxfield = generic_xderiv_3(dcoefs[plane], mesh)
bxcompp = generic_xderiv_3(dcoefs[plane], plane_dat[:,:3])

fig2, ax2 = plt.subplots(1,2)
cax2 = ax2[0].scatter(plane_dat[:,0]*normalizer[0], plane_dat[:,2]*normalizer[0], c = plane_dat[:,6]*normalizer[1]*1e6)#, vmin = min(bzfield*normalizer[1]*1e6), vmax = max(bzfield*normalizer[1]*1e6))
plt.colorbar(cax2)
ax2[0].set_title(r'$B_{z, measured}$ $[\mu T]$')
cax22 = ax2[1].scatter(mesh[:,0]*normalizer[0], mesh[:,2]*normalizer[0], c = bzfield*normalizer[1]*1e6)
plt.colorbar(cax22)
ax2[1].set_title(r'$B_{z, fit}$ $[\mu T]$')
plt.tight_layout()

fig3, ax3 = plt.subplots(1,2)
cax3 = ax3[0].scatter(plane_dat[:,0]*normalizer[0], plane_dat[:,2]*normalizer[0], c = plane_dat[:,4]*normalizer[1]*1e6)#, vmin = min(bxfield*normalizer[1]*1e6), vmax = max(bxfield*normalizer[1]*1e6))
plt.colorbar(cax3)
ax3[0].set_title(r'$B_{x, measured}$ $[\mu T]$')
cax33 = ax3[1].scatter(mesh[:,0]*normalizer[0], mesh[:,2]*normalizer[0], c = bxfield*normalizer[1]*1e6)
plt.colorbar(cax33)
ax3[1].set_title(r'$B_{x, fit}$ $[\mu T]$')
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
    

fig5 = plt.figure()
plt.plot(onaxis[:,1]*normalizer[0], method1err*normalizer[1], label = 'Method 1')
plt.plot(onaxis[:,1]*normalizer[0], method2err*normalizer[1], label = 'Method 2')
plt.plot(onaxis[:,1]*normalizer[0], xmethoderr*normalizer[1], label = r'$B_{x}$')
plt.plot(onaxis[:,1]*normalizer[0], zmethoderr*normalizer[1], label = r'$B_{z}$')
plt.xlabel('y [m]')
#plt.title(r'$\sqrt{\frac{\left(B_{x}-\overline{B}_{x, true}\right)^{2}}{N}}$ + $\sqrt{\frac{\left(B_{z}-\overline{B}_{z, true}\right)^{2}}{N}}$ $[\mu T]$', fontsize = 16)
plt.legend()
plt.tight_layout()

fig6 = plt.figure()
plt.plot(np.arange(len(dcoefs[1:])), dcoefs[1:])

fig7, ax7 = plt.subplots(1,2)

c7 = ax7[0].scatter(plane_dat[:,0], plane_dat[:,2], c = bzcompp-plane_dat[:,6])
plt.colorbar(c7)


plt.show()