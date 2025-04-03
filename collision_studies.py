import numpy as np 
import matplotlib.pyplot as plt 

from thermal_distribution import gen_velocities
from move_stepper import neutron_stepper

neutrons = 50
m_neutron = 1.6749286e-27
temperature = 4e-3  # Kelvin
Joule_eV = 1.619e-19

width = 1e-9
bounds = [10, 0, 0]# [0.0375, 0, 0.0375]
time_elapse = 1#0.2


init_velocities = gen_velocities(N = neutrons, T = temperature, m = m_neutron)
init_speed = np.sqrt(init_velocities[:,0]**2 + init_velocities[:,1]**2 + init_velocities[:,2]**2)

init_positions = np.zeros((neutrons, 3))
init_positions[:,0] = np.random.normal(0, width, np.shape(init_positions[:,0]))
init_positions[:,1] = np.zeros(np.shape(init_positions[:,0]))
init_positions[:,2] = np.random.normal(0, width, np.shape(init_positions[:,0]))

initial_log = np.zeros((neutrons,6))
initial_log[:,:3] = init_positions
initial_log[:,3:] = init_velocities

neutron_log = []

radii = np.linspace(1, 0.01, 3)
d = radii*2
avg_init_speed = sum(init_speed)/neutrons
avg_collision_time = d/avg_init_speed


sim_collisions = np.zeros((neutrons, len(radii)))
final_points = np.zeros((3,neutrons,len(radii)))

for i in range(len(radii)):
    for j in range(neutrons):
        log, bounces, lost = neutron_stepper(initial_log[j], t0 = 0, t_span = time_elapse, dt = avg_collision_time[i]*1e-3, bounds = [radii[i]])
        sim_collisions[j][i] = bounces
        for k in range(3):
            final_points[k][j][i] = log[-1, k+1]


fig = plt.figure()

plt.scatter(radii, sum(sim_collisions)/(neutrons*time_elapse))
plt.plot(radii, 1/avg_collision_time)
plt.xlabel('Radius R [m]')
plt.ylabel(r'$\overline{\alpha}$ [1/s]')
plt.title('Average Collision Time v. Cylinder Radius')


fig1 = plt.figure()
comps = ['x', 'y', 'z']
for k in range(3):
    rms_distance_covered = np.sqrt(sum(final_points[k]**2)/neutrons)
    plt.scatter(radii, rms_distance_covered, label=f'{comps[k]}')

plt.xlabel('Radius R [m]')
plt.ylabel('RMS Distance')
plt.legend()

plt.show()