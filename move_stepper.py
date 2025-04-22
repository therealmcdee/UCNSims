import numpy as np 
from scipy.spatial.transform import Rotation
from scipy.interpolate import RegularGridInterpolator as RGI

from simple_field import generic_xderiv_3, generic_zderiv_3
from rgi_field import generate_field_interp
from thermal_distribution import gen_velocities
from stc_coil_function import gen_stc_interpolation

m_neutron = 1.6749286e-27
g_neutron = -1.913
nuclear_mag = 5.051e-27
h = 6.626e-34
neutron_mag = g_neutron*nuclear_mag
gyromagnetic_ratio = 2*(g_neutron*nuclear_mag)/(h/(2*np.pi))

def neutron_stepper(N, T, m, t_span, bounds = None, decay = False, spin_tracking = False, coef = None, fieldfile = None, stcfield = False):
    if stcfield == True:
        Bx, By, Bz  = gen_stc_interpolation('SL', np.array([10, 30, 30, 100]))
    init_velocities = gen_velocities(N, T, m)
    init_speed = np.sqrt(init_velocities[:,0]**2 + init_velocities[:,1]**2 + init_velocities[:,2]**2)
    avg_speed = sum(init_speed)/N
    init_positions = np.zeros((N, 3))
    width = bounds[0]*1e-9
    diameter = bounds[0]*2
    avg_collision_time = diameter/avg_speed
    dt = avg_collision_time*1e-3
    if len(bounds) == 1:
        init_velocities[:,1] = abs(init_velocities[:,1])
        init_positions[:,0] = np.random.uniform(0, width, np.shape(init_positions[:,0]))
        init_positions[:,1] = np.zeros(np.shape(init_positions[:,0]))
        init_positions[:,2] = np.random.uniform(0, width, np.shape(init_positions[:,0]))
    elif len(bounds) != 1:
        init_positions[:,0] = np.random.uniform(0, width/np.sqrt(2), np.shape(init_positions[:,0]))
        init_positions[:,1] = np.random.uniform(0, width/np.sqrt(2), np.shape(init_positions[:,0]))
        init_positions[:,2] = np.random.uniform(bounds[2], bounds[1], np.shape(init_positions[:,0]))
    initial_log = np.zeros((N,6))
    initial_log[:,:3] = init_positions
    initial_log[:,3:] = init_velocities
    neutron_log = []
    decays = 0
    a = np.array([0, 0, -9.81])
    for n in range(N):
        s_initial = np.array([0, 0, -1])
        t = 0
        log = []
        decayed = 0
        initial = initial_log[n]
        if spin_tracking == False:
            next = np.zeros(6)
            log.append([t, initial[0], initial[1], initial[2], initial[3], initial[4], initial[5]])
        elif spin_tracking == True:
            if coef is None:
                coef = np.array([0, 0, 1e-9, 1e-6, 0, 0, 0, 0, 0])
                #coef = np.array([0.00000000e+00,  7.19805747e-08,  8.60186407e-06, -4.91957698e-06, -1.11316291e-05,  1.66550372e-03,  1.39837144e-05,  8.76780834e-06, -5.45416237e-07]
            next = np.zeros(9)
            s_next = np.zeros(3)
            spin_step = dt
            log.append([t, initial[0], initial[1], initial[2], initial[3], initial[4], initial[5], s_initial[0], s_initial[1], s_initial[2]])
        while t < t_span:
            bounces = 0
            if decay == True:
                if np.random.rand()<dt/880:
                    decayed += 1
                    break
            for i in range(3):
                next[i] = initial[i] + initial[i+3]*dt + 0.5*a[i]*(dt**2)
                next[i+3] = initial[i+3] + a[i]*dt
            if bounds != None:
                if len(bounds) == 1:
                    R = bounds[0]
                    if np.sqrt(next[0]**2 + next[2]**2)>=R:
                        r = [initial[0], initial[2]]
                        p = np.sqrt(r[0]**2 + r[1]**2)
                        v = [initial[3], initial[5]]
                        vr = np.sqrt(v[0]**2 + v[1]**2)
                        time_to_hit_wall = (R-p)/vr             #### have to add gravity into time elapsed!!!
                        rwall = r + time_to_hit_wall*vr
                        normal = -rwall/R
                        ndotv = v[0]*normal[0] + v[1]*normal[1]
                        vrefl = v - 2*ndotv*normal
                        vnew = [vrefl[0], initial[4], vrefl[1]]
                        rnew = [rwall[0], initial[1], rwall[1]]
                        for i in range(3):
                            next[i] = rnew[i] + vnew[i]*(dt-time_to_hit_wall) + 0.5*a[i]*(dt-time_to_hit_wall)**2
                            next[i+3] = vnew[i] + a[i]*(dt-time_to_hit_wall)
                        bounces += 1
                elif len(bounds) > 1:
                    R = bounds[0]
                    h1 = bounds[1]
                    h2 = bounds[2]
                    if np.sqrt(next[0]**2 + next[1]**2)>=R:
                        r = [initial[0], initial[1]]
                        p = np.sqrt(r[0]**2 + r[1]**2)
                        v = [initial[3], initial[4]]
                        vr = np.sqrt(v[0]**2 + v[1]**2)
                        time_to_hit_wall = (R-p)/vr            
                        rwall = r + time_to_hit_wall*vr
                        normal = -rwall/R
                        ndotv = v[0]*normal[0] + v[1]*normal[1]
                        vrefl = v - 2*ndotv*normal
                        vnew = [vrefl[0], vrefl[1], initial[5]]
                        rnew = [rwall[0], rwall[1], initial[2]]
                        for i in range(3):
                            next[i] = rnew[i] + vnew[i]*(dt-time_to_hit_wall) + 0.5*a[i]*(dt-time_to_hit_wall)**2
                            next[i+3] = vnew[i] + a[i]*(dt-time_to_hit_wall)
                        bounces += 1
                    elif next[2]>h1:
                        time_to_hit_wall = (h1-initial[2])/initial[5]
                        print(h1-initial[2])
                        r_evo = initial[:3] + initial[3:6]*time_to_hit_wall + 0.5*a*(time_to_hit_wall**2)
                        vnew = [initial[3], initial[4], -initial[5]]
                        for i in range(3):
                            next[i] = r_evo[i] + vnew[i]*(dt-time_to_hit_wall) + 0.5*a[i]*(dt-time_to_hit_wall)**2
                            next[i+3] = vnew[i] + a[i]*(dt-time_to_hit_wall)
                        bounces += 1
                    elif next[2]<h2:
                        time_to_hit_wall = (h2-initial[2])/initial[5]
                        r_evo = initial[:3] + initial[3:6]*time_to_hit_wall + 0.5*a*(time_to_hit_wall**2)
                        vnew = [initial[3], initial[4], -initial[5]]
                        for i in range(3):
                            next[i] = r_evo[i] + vnew[i]*(dt-time_to_hit_wall) + 0.5*a[i]*(dt-time_to_hit_wall)**2
                            next[i+3] = vnew[i] + a[i]*(dt-time_to_hit_wall)
                        bounces += 1
            if spin_tracking == True:
                if stcfield == False:
                    B = np.zeros(3)
                    B[0] = generic_xderiv_3(coef, np.array([next[:3]]))
                    B[2] = generic_zderiv_3(coef, np.array([next[:3]]))
                if next[1] >= 0.8:
                    B = np.zeros(3)
                else:
                    B = np.zeros(3)
                    B[0] = Bx(next[:3])
                    B[1] = By(next[:3])
                    B[2] = Bz(next[:3])
                    print(B)
                dx = np.zeros(3)
                k1 = gyromagnetic_ratio*np.cross(s_initial, B)*spin_step/2
                #r = Rotation.from_rotvec(k1)
                s1 = s_initial + k1*spin_step/2 #np.matmul(r.as_matrix(), s_initial)
                k2 = gyromagnetic_ratio*np.cross(s1, B)*spin_step/2
                #r = Rotation.from_rotvec(k2)
                s2 = s_initial + k2*spin_step/2 #np.matmul(r.as_matrix(), s1)
                k3 = gyromagnetic_ratio*np.cross(s2, B)*spin_step
                #r = Rotation.from_rotvec(k3)
                s3 = s_initial + k3*spin_step #np.matmul(r.as_matrix(), s2)
                k4 = gyromagnetic_ratio*np.cross(s3, B)*spin_step
                dx = (k1 + 2*k2 + 2*k3 + k4)*(spin_step/6)
                #r = Rotation.from_rotvec(dx)
                s_next = s_initial + dx #np.matmul(r.as_matrix(), s_initial)
                length = np.sqrt(s_next[0]**2 + s_next[1]**2 + s_next[2]**2)
                s_next = s_next/length
            t += dt
            if spin_tracking == False:
                log.append([t, next[0], next[1], next[2], next[3], next[4], next[5]])
            elif spin_tracking == True:
                log.append([t, next[0], next[1], next[2], next[3], next[4], next[5], s_next[0], s_next[1], s_next[2]])
            initial = next
            if spin_tracking == True:
                s_initial = s_next
        if decayed == 0:
            neutron_log.append([])
            neutron_log[n] = log
        else:
            continue 
        print(n)
    return np.asarray(neutron_log)
            
    
    