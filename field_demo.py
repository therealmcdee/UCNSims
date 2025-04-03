import numpy as np 
import matplotlib.pyplot as plt

import simple_field

coef = [0, 0, 0, 0, 1, 0, 0, 0, 0]
N = 100
cnt = 0
plane = np.zeros((N*N, 3))
vals = np.linspace(-1, 1, N)
for i in range(N):
    for j in range(N):
        plane[cnt][0] = vals[i]
        plane[cnt][1] = 0
        plane[cnt][2] = vals[j]
        cnt += 1
plane = np.asarray(plane)
        

bxfield = simple_field.generic_xderiv_3(coef, plane)
bzfield = simple_field.generic_zderiv_3(coef, plane)

fig = plt.figure()

c = plt.scatter(plane[:,0], plane[:,2], c = bxfield)
plt.title('Bx')
plt.colorbar(c)

fig2 = plt.figure()
c2 = plt.scatter(plane[:,0], plane[:,2], c = bzfield)
plt.title('Bz')
plt.colorbar(c2)

plt.show()