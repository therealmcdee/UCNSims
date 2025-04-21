import numpy as np 
import matplotlib.pyplot as plt 

from rgi_field import generate_field_interp

fname = 'STC_3DMAPS_2025/Ccoils/SUC_100mA_subtracted.txt'
#fname = 'STC_3DMAPS_2025/response_curves/SUC_response.txt'

rflag = 0
rbool = False
if fname.split('/')[1] == 'response_curves':
    rflag = 1
    rbool = True


odata = np.loadtxt(fname, delimiter = ' ')
data = odata[:(len(odata)-1)]

Bx, By, Bz = generate_field_interp(data, rbool)

#odata = np.loadtxt(fname, delimiter = ' ')
#data = odata[:(len(odata)-1)]
onaxis = data[(data[:,0]==-0.01) & (data[:,2]==-0.01)]

N = 100
test = np.zeros((N**2, 3))
x = np.linspace(-0.01, 0.01, N)
cnt = 0
for i in range(N):
    for j in range(N):
        test[cnt][0] = x[i]
        test[cnt][1] = 0.08
        test[cnt][2] = x[j]
        cnt += 1

bztest = Bz(test)
bxtest = Bx(test)
bytest = By(test)


fig = plt.figure()

ct = plt.scatter(test[:,0], test[:,2], c = bztest)
plt.colorbar(ct, label = r'$B_{z}$')

fig0 = plt.figure()
ct0 = plt.scatter(test[:,0], test[:,2], c = bxtest)
plt.colorbar(ct0, label = r'$B_{x}$')

axis = np.zeros((N, 3))
y = np.linspace(0, 0.5, N)
for i in range(len(y)):
    axis[i][0] = -0.01
    axis[i][1] = y[i]
    axis[i][2] = -0.01

bzonax = Bz(axis)
bxonax = Bx(axis)
byonax = By(axis)

fig1 = plt.figure()
plt.plot(axis[:,1], bzonax, label = r'$B_{z}$')
plt.plot(axis[:,1], bxonax, label = r'$B_{x}$')
plt.plot(axis[:,1], byonax, label = r'$B_{y}$')

plt.scatter(onaxis[:,1], onaxis[:,4 - rflag])
plt.scatter(onaxis[:,1], onaxis[:,5 - rflag])
plt.scatter(onaxis[:,1], onaxis[:,6 - rflag])

plt.legend()

bxcompare = Bx(data[:,:3])
bycompare = By(data[:,:3])
bzcompare = Bz(data[:,:3])

fig2, ax2 = plt.subplots(1, 3)
nbins = 30
ax2[0].hist(bxcompare - data[:,4 - rflag], nbins)
ax2[1].hist(bycompare - data[:,5 - rflag], nbins)
ax2[2].hist(bzcompare - data[:,6 - rflag], nbins)

plt.tight_layout()



plt.show()