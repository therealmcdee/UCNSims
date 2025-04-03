import numpy as np 
from scipy.spatial.transform import Rotation

from simple_field import generic_xderiv_3, generic_zderiv_3

m_neutron = 1.6749286e-27
g_neutron = -1.913
nuclear_mag = 5.051e-27
h = 6.626e-34
neutron_mag = g_neutron*nuclear_mag
gyromagnetic_ratio = 2*(g_neutron*nuclear_mag)/(h/(2*np.pi))

def neutron_stepper(initial, t0, t_span, dt, bounds = None, decay = False, spin_tracking = False, coef = None):
    log = []
    s_initial = np.array([0, 0, -1])
    t = 0
    a = [0, 0, -9.81]
    if spin_tracking == False:
        next = np.zeros(6)
        log.append([t, initial[0], initial[1], initial[2], initial[3], initial[4], initial[5]])
    elif spin_tracking == True:
        if coef == None:
            coef = np.array([0.00000000e+00, 8.36279257e-03, 9.14379959e-01, -1.60091662e-01, -3.61960986e-01, 1.65867711e+01, 1.39017307e-01, 1.46966811e+01, -9.29917816e-01])
        next = np.zeros(9)
        s_next = np.zeros(3)
        spin_step = dt
        log.append([t, initial[0], initial[1], initial[2], initial[3], initial[4], initial[5], s_initial[0], s_initial[1], s_initial[2]])
    bounces = 0
    lost = 0
    while t < t_span:
        if decay == True:
            if np.random.rand()<dt/880:
                lost = 1
                break
        for i in range(3):
            next[i] = initial[i] + initial[i+3]*dt + 0.5*a[i]*(dt**2)
            next[i+3] = initial[i+3] + a[i]*dt
        if spin_tracking == True:
            B = np.zeros(3)
            B[0] = generic_xderiv_3(coef, np.array([initial[:3]]))
            B[2] = generic_zderiv_3(coef, np.array([initial[:3]]))
            dx = np.zeros(3)
            k1 = gyromagnetic_ratio*np.cross(s_initial, B)*spin_step/2
            r = Rotation.from_rotvec(k1)
            s1 = np.matmul(r.as_matrix(), s_initial)
            k2 = gyromagnetic_ratio*np.cross(s1, B)*spin_step/2
            r = Rotation.from_rotvec(k2)
            s2 = np.matmul(r.as_matrix(), s1)
            k3 = gyromagnetic_ratio*np.cross(s2, B)*spin_step
            r = Rotation.from_rotvec(k3)
            s3 = np.matmul(r.as_matrix(), s2)
            k4 = gyromagnetic_ratio*np.cross(s3, B)*spin_step
            dx = (k1 + 2*k2 + 2*k3 + k4)*(spin_step/6)
            r = Rotation.from_rotvec(dx)
            s_next = np.matmul(r.as_matrix(), s_initial)
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
                        next[i+3] = vnew[i] + a[i]*(dt-time_to_hit_wall)**2
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
                    time_to_hit_wall = (R-p)/vr             #### have to add gravity into time elapsed!!!
                    rwall = r + time_to_hit_wall*vr
                    normal = -rwall/R
                    ndotv = v[0]*normal[0] + v[1]*normal[1]
                    vrefl = v - 2*ndotv*normal
                    vnew = [vrefl[0], vrefl[1], initial[5]]
                    rnew = [rwall[0], rwall[1], initial[2]]
                    for i in range(3):
                        next[i] = rnew[i] + vnew[i]*(dt-time_to_hit_wall) + 0.5*a[i]*(dt-time_to_hit_wall)**2
                        next[i+3] = vnew[i] + a[i]*(dt-time_to_hit_wall)**2
                    bounces += 1
#                if next[2]>h1:
# HERE IS WHERE YOU ARE WORKING IN PROGRESS
        t += dt
        if spin_tracking == False:
            log.append([t, next[0], next[1], next[2], next[3], next[4], next[5]])
        elif spin_tracking == True:
            log.append([t, next[0], next[1], next[2], next[3], next[4], next[5], s_next[0], s_next[1], s_next[2]])
        initial = next
        if spin_tracking == True:
            s_initial = s_next 
    return np.asarray(log), bounces, lost
            
    
    