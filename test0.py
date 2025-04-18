import numpy as np 
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


from thermal_distribution import gen_velocities
from move_stepper import neutron_stepper

neutrons = 30
m_neutron = 1.6749286e-27
temperature = 4e-3  # Kelvin
Joule_eV = 1.619e-19


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
bounds = [0.1, 0, 0]# [0.0375, 0, 0.0375]

d = bounds[0]*2
avg_init_speed = sum(init_speed)/neutrons
avg_collision_time = d/avg_init_speed
print(avg_collision_time)

decay = 0
for i in range(neutrons):
    log, bounces, dnf = neutron_stepper(initial_log[i], t0 = 0, t_span = 0.2, dt = avg_collision_time*5e-4, bounds = bounds)#,  decay = False, spin_tracking = False)
    if dnf == 0:
        neutron_log.append(log)
    else:
        decay += 1
neutron_log = np.asarray(neutron_log)

print(decay)


final_points = neutron_log[:, -1, 1:4]
print(len(final_points[(final_points[:,1]>=0.75)])/neutrons)

fig = plt.figure()
plt.hist(init_speed, bins = 50, label = round(sum(init_speed)/neutrons,3))
plt.legend()

fig1 = plt.figure()
for i in range(round(np.ceil(neutrons*1/10))):
    plt.scatter(neutron_log[i,:,1], neutron_log[i,:,3])

fig2 = plt.figure()
avlabs = ['<x>','<y>', '<z>']
neutron_pos = np.zeros((len(neutron_log[0,:,0]),3))
for j in range(len(neutron_log[0,:,0])):
    neutron_pos[j] = sum(neutron_log[:, j, 1:4])/len(neutron_log[:, j, 1:4])
for k in range(3):
    plt.plot(neutron_log[0,:,0], neutron_pos[:,k],label=f'{avlabs[k]}')
plt.legend()
plt.xlabel('Time [s]')
plt.ylabel('Average Position Over Neutrons [m]')


plt.show()