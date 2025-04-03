import numpy as np 
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import minimize


from thermal_distribution import gen_velocities
from move_stepper import neutron_stepper

def exponential(x, t):
    return x[0]*np.exp(-t/x[1])

def residual(x, t, data):
    return sum((exponential(x, t) - data)**2)

neutrons = 500
m_neutron = 1.6749286e-27
temperature = 1e-3  # Kelvin
Joule_eV = 1.619e-19


init_velocities = gen_velocities(N = neutrons, T = temperature, m = m_neutron)
init_speed = np.sqrt(init_velocities[:,0]**2 + init_velocities[:,1]**2 + init_velocities[:,2]**2)

init_positions = np.zeros((neutrons, 3))
width = 0.05
init_positions[:,0] = np.random.normal(0, width, np.shape(init_positions[:,0]))
init_positions[:,1] = np.zeros(np.shape(init_positions[:,0]))
init_positions[:,2] = np.random.normal(0, width, np.shape(init_positions[:,0]))

initial_log = np.zeros((neutrons,6))
initial_log[:,:3] = init_positions
initial_log[:,3:] = init_velocities

neutron_log = []
bounds = [10, 0, 0]# [0.0375, 0, 0.0375]

d = bounds[0]*2
avg_init_speed = sum(init_speed)/neutrons
avg_collision_time = d/avg_init_speed
print(avg_collision_time)

nsamps = 30
tsmall = 5
sample_times = np.linspace(tsmall, 2500, nsamps)
survivors = np.zeros(nsamps)
decays = np.zeros(nsamps)

decay = 0
for j in range(nsamps):
    for i in range(neutrons):
        log, bounces, dnf = neutron_stepper(initial_log[i], t0 = 0, t_span = sample_times[j], dt = tsmall, decay = True)#, bounds = bounds)
        if dnf == 0:
            survivors[j] += 1
        else:
            decays[j] += 1
            
res = minimize(residual, x0 = [neutrons, 1000], args = (sample_times, survivors))
print(res)

fig = plt.figure()

plt.scatter(sample_times, survivors)
plt.plot(sample_times, exponential(res.x, sample_times), color = 'r')

plt.show()

