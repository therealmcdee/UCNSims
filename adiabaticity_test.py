import numpy as np 
import matplotlib.pyplot as plt

from simple_field import generic_field_3, generic_xderiv_3, generic_zderiv_3
from move_stepper import neutron_stepper
from thermal_distribution import gen_velocities

neutrons = 100
temperature = 4e-3
m_neutron = 1.6749286e-27

init_velocities = gen_velocities(N = neutrons, T = temperature, m = m_neutron)
init_velocities[:,1] = abs(init_velocities[:,1])
init_speed = np.sqrt(init_velocities[:,0]**2 + init_velocities[:,1]**2 + init_velocities[:,2]**2)

init_positions = np.zeros((neutrons, 3))
width = 1e-9
init_positions[:,0] = np.random.normal(0, width, np.shape(init_positions[:,0]))#np.random.uniform(-np.sqrt(0.035), np.sqrt(0.035), np.shape(neutrons,)) // np.random.normal(0, width, np.shape(init_positions[:,0]))
init_positions[:,1] = np.zeros(np.shape(init_positions[:,0]))
init_positions[:,2] = np.random.normal(0, width, np.shape(init_positions[:,0]))#np.random.uniform(-np.sqrt(0.035), np.sqrt(0.035), np.shape(neutrons,)) // np.random.normal(0, width, np.shape(init_positions[:,0]))

initial_log = np.zeros((neutrons,6))
initial_log[:,:3] = init_positions
initial_log[:,3:] = init_velocities

neutron_log = []
bounds = [0.05]#, 0, 0]# [0.0375, 0, 0.0375]

d = bounds[0]*2
avg_init_speed = sum(init_speed)/neutrons
avg_collision_time = d/avg_init_speed
print(avg_collision_time)

decay = 0
for i in range(neutrons):
    log, bounces, dnf = neutron_stepper(initial_log[i], t0 = 0, t_span = 0.1, dt = avg_collision_time*1e-3, bounds = bounds, decay = False, spin_tracking = True)
    if dnf == 0:
        neutron_log.append(log)
    else:
        decay += 1
neutron_log = np.asarray(neutron_log)

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