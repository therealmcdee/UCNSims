import numpy as np 
import matplotlib.pyplot as plt 
from scipy.interpolate import RegularGridInterpolator as RGI 

odata = np.loadtxt('NUB_150mA_subtracted.txt', delimiter = ' ')
data = odata[:(len(odata)-1)]

new = np.zeros(np.shape(data))

xsort = np.sort(data, axis = 0)

x = np.ndarray.flatten(xsort[:,0])
y = np.ndarray.flatten(xsort[:,1])
z = np.ndarray.flatten(xsort[:,2])

xv = []
yv = []
zv = []
for i in range(len(x)):
    if np.any(xv==x[i]) == True:
        continue
    else:
        xv.append(x[i])
for i in range(len(x)):
    if np.any(yv == y[i]) == True:
        continue
    else:
        yv.append(y[i])
for i in range(len(x)):
    if np.any(zv == z[i]) == True:
        continue
    else:
        zv.append(z[i])

xg, yg, zg = np.meshgrid(xv,yv,zv)

bxgrid = np.zeros((len(xv), len(yv), len(zv)))
bygrid = np.zeros((len(xv), len(yv), len(zv)))
bzgrid = np.zeros((len(xv), len(yv), len(zv)))
grid = np.zeros((len(xv)*len(yv)*len(zv), 3))
cnt = 0
for i in range(len(xv)):
    for j in range(len(yv)):
        for k in range(len(zv)):
            for n in range(len(data)):
                if data[n,0]==xv[i] and data[n, 1]==yv[j] and data[n, 2]==zv[k]:
                    bxgrid[i][j][k] = data[n, 4]
                    bygrid[i][j][k] = data[n, 5]
                    bzgrid[i][j][k] = data[n, 6]
                    grid[cnt][0] = xv[i]
                    grid[cnt][1] = yv[j]
                    grid[cnt][2] = zv[k]
                    cnt += 1
                    

Bx = RGI([xv, yv, zv], bxgrid)
By = RGI([xv, yv, zv], bygrid)
Bz = RGI([xv, yv, zv], bzgrid)

fig = plt.figure()
ax = fig.add_subplot(projection = '3d')
ct = ax.scatter(grid[:,0], grid[:,1], grid[:,2], c = bzgrid)
plt.colorbar(ct)


N = 100
test = np.zeros((N**2, 3))
x = np.linspace(-0.01, 0.01, N)
cnt = 0
for i in range(N):
    for j in range(N):
        test[cnt][0] = x[i]
        test[cnt][1] = 0.1
        test[cnt][2] = x[j]
        cnt += 1

bztest = By(test)

fig = plt.figure()

ct = plt.scatter(test[:,0], test[:,2], c = bztest)
plt.colorbar(ct)


axis = np.zeros((N, 3))
y = np.linspace(0, 0.30, N)
for i in range(len(y)):
    axis[i][0] = 0
    axis[i][1] = y[i]
    axis[i][2] = 0

bzonax = Bz(axis)
fig1 = plt.figure()
plt.scatter(axis[:,1], bzonax)


plt.show()