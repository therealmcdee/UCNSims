import numpy as np 
import matplotlib.pyplot as plt 

from move_stepper import neutron_stepper

neutrons = 50
temp = 4e-3
m_neutron = 1.6749286e-27

bounds = [0.25]#, 0.15, 0.05]
print(len(bounds))

neutron_log = neutron_stepper(neutrons, temp, m_neutron, t_span = 0.1, bounds = bounds, decay = False, spin_tracking = True)



fig, ax = plt.subplots(2, 1, sharex=True)
spins = [r'$S_{x}$', r'$S_{y}$', r'$S_{z}$']
for i in range(3):
    ax[1].plot(neutron_log[0,:,0], neutron_log[0,:,-(i+1)],label=f'{spins[-(i+1)]}')
ax[1].legend()

ax[0].plot(neutron_log[0,:,0], np.sqrt(neutron_log[0,:,-3]**2 + neutron_log[0,:,-2]**2 + neutron_log[0,:,-1]**2)-1)
ax[0].axhline(0, ls='--',c='k')
ax[0].set_ylabel(r'$|S|-1$')
ax[1].set_ylabel(r'$S_{i}$', fontsize = 14)
ax[1].set_xlabel('Time [s]')
plt.suptitle('Single Particle Tracking')

plt.tight_layout()

fig, ax = plt.subplots(2, 1, sharex=True)
groupspins = [r'$<S_{x}>$', r'$<S_{y}>$', r'$<S_{z}$>']
for i in range(3):
    ax[1].plot(neutron_log[0,:,0], sum(neutron_log[:,:,-(i+1)])/neutrons,label=f'{groupspins[-(i+1)]}')
ax[1].legend()
ax[0].plot(neutron_log[0,:,0], sum(np.sqrt(neutron_log[:,:,-3]**2 + neutron_log[:,:,-2]**2 + neutron_log[:,:,-1]**2)-1)/neutrons)
ax[0].axhline(0, ls='--',c='k')
ax[0].set_ylabel(r'$|S|-1$')
ax[1].set_ylabel(r'$Polarization$', fontsize = 14)
ax[1].set_xlabel('Time [s]')
plt.suptitle(f'n = {neutrons} UCN in {bounds[0]} m guide')
plt.tight_layout()


fig = plt.figure()
for i in range(neutrons):
    plt.plot(neutron_log[i,:,1], neutron_log[i,:,3])
plt.xlabel('x [m]')
plt.ylabel('z [m]')


plt.show()