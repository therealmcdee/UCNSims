import numpy as np 
import matplotlib.pyplot as plt 
from scipy.interpolate import RegularGridInterpolator as RGI 


def generate_field_interp(filename):
    odata = np.loadtxt(filename, delimiter = ' ')
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

    return Bx, By, Bz