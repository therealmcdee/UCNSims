import numpy as np 
import matplotlib.pyplot as plt 

from thermal_distribution import gen_velocities

Joule_eV = 1.619e-19

m_neutron = 1.6749286e-27
g_neutron = -1.913
nuclear_mag = 5.051e-27
h = 6.626e-34
neutron_magnetic_dipole = g_neutron*nuclear_mag
neutron_gyromagnetic_ratio = 2*(g_neutron*nuclear_mag)/(h/(2*np.pi))


m_mercury = 3.3e-25
mercury_gyromagnetic_ratio = 2*np.pi*(7.590118e6)


print(mercury_gyromagnetic_ratio)
print(neutron_gyromagnetic_ratio)

neutrons = 10000
mercury = 10000

neutron_vel = gen_velocities(N = neutrons, T = 4e-3, m = m_neutron)
mercury_vel = gen_velocities(N = mercury, T = 293, m = m_mercury)

neutron_speed = np.sqrt(neutron_vel[:,0]**2 + neutron_vel[:,1]**2 + neutron_vel[:,2]**2)
mercury_speed = np.sqrt(mercury_vel[:,0]**2 + mercury_vel[:,1]**2 + mercury_vel[:,2]**2)

neutron_energy = 0.5*m_neutron*(neutron_speed**2)/Joule_eV
mercury_energy = 0.5*m_mercury*(mercury_speed**2)/Joule_eV

fig0 = plt.figure()
bins = 50
plt.hist(neutron_energy, bins)
plt.title(f'UCN <v> = {sum(neutron_speed)/neutrons} m/s')
plt.xlabel('Kinetic Energy [eV]')
#plt.legend()
plt.tight_layout()

fig2 = plt.figure()
plt.hist(mercury_energy, bins)
plt.title(f'199Hg <v> = {sum(mercury_speed)/mercury} m/s')
plt.xlabel('Kinetic Energy [eV]')
#plt.legend()
plt.tight_layout()

plt.show()