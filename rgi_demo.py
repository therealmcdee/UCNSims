import numpy as np 
import matplotlib.pyplot as plt 

from rgi_field import generate_field_interp


Bx, By, Bz = generate_field_interp('NUB_150mA_subtracted.txt')


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

bztest = Bz(test)

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