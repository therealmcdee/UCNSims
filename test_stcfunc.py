import numpy as np 
import matplotlib.pyplot as plt 

from rgi_field import generate_field_interp
from stc_coil_function import gen_stc_interpolation

coils = 'NU'

curr = (1e-3)*np.array([8, 30, 30, 100])
print(curr)

stcx, stcy, stcz = gen_stc_interpolation(coils, curr)

N = 100
axis = np.zeros((N, 3))
y = np.linspace(0, 0.8, N)
axis[:,0] = np.zeros(N)
axis[:,1] = y
axis[:,2] = np.zeros(N)

print(axis)

print(stcx([0, 0.4, 0]))
print(stcy([0, 0.4, 0]))
print(stcz([0, 0.4, 0]))

bx = stcx(axis)
by = stcy(axis)
bz = stcz(axis)

fig = plt.figure()
plt.plot(y, bx, label = 'bx')
plt.plot(y, by, label = 'by')
plt.plot(y, bz, label = 'bz')
plt.legend()
plt.show()