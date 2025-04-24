import numpy as np 
import matplotlib.pyplot as plt 
import os

from rgi_field import generate_field_interp


def gen_stc_interpolation(coilset, currents):
    parent_dir = 'STC_3DMAPS_2025/response_curves'
        ## A, B, C, D
    Afield = []
    Bfield = []
    Cfield = []
    Dfield = []


    for i in os.listdir(parent_dir):
        if i[:2] == coilset:
            if i[2] == 'A':
                coilstart = 5/100
                start = 0.65
                end = 0.75
                L = 0.1
                response = np.loadtxt(parent_dir+'/'+i, delimiter = ' ')
                Afield = np.zeros_like(response)
                Afield[:,:3] = response[:,:3]
                Afield[:,1] = response[:,1] + start - coilstart
                Afield[:,3:6] = currents[0]*response[:,3:6]
            elif i[2] == 'B':
                coilstart = 7.62/100
                start = 0.5
                end = 0.65
                L = 0.15
                response = np.loadtxt(parent_dir+'/'+i, delimiter = ' ')
                Bfield = np.zeros_like(response)
                Bfield[:,:3] = response[:,:3]
                Bfield[:,1] = response[:,1] + start - coilstart
                Bfield[:,3:6] = currents[1]*response[:,3:6]
            elif i[2] == 'C':
                coilstart = 12.5/100
                start = 0.25
                end = 0.5
                L = 0.25
                response = np.loadtxt(parent_dir+'/'+i, delimiter = ' ')
                Cfield = np.zeros_like(response)
                Cfield[:,:3] = response[:,:3]
                Cfield[:,1] = response[:,1] + start - coilstart
                Cfield[:,3:6] = currents[2]*response[:,3:6]
            elif i[2] == 'D':
                coilstart = 12.5/100
                start = 0.
                end = 0.25
                L = 0.25
                response = np.loadtxt(parent_dir+'/'+i, delimiter = ' ')
                Dfield = np.zeros_like(response)
                Dfield[:,:3] = response[:,:3]
                Dfield[:,1] = response[:,1] + start - coilstart
                Dfield[:,3:6] = currents[3]*response[:,3:6]

    Aonax = Afield[(Afield[:,0]==0) & (Afield[:,2]==0)]
    Bonax = Bfield[(Bfield[:,0]==0) & (Bfield[:,2]==0)]
    Conax = Cfield[(Cfield[:,0]==0) & (Cfield[:,2]==0)]
    Donax = Dfield[(Dfield[:,0]==0) & (Dfield[:,2]==0)]
    print(Donax)
    ai_x, ai_y, ai_z = generate_field_interp(Afield, response_curve = True)
    bi_x, bi_y, bi_z = generate_field_interp(Bfield, response_curve = True)
    ci_x, ci_y, ci_z = generate_field_interp(Cfield, response_curve = True)
    di_x, di_y, di_z = generate_field_interp(Dfield, response_curve = True)

    N2 = 31
    cnt = 0
    same_grid = np.zeros((N2**3, 6))
    yax = np.linspace(min(Donax[:,1]), max(Aonax[:,1]), N2)
    pax = np.linspace(-0.01, 0.01, N2)
    for i in range(len(pax)):
        for j in range(len(yax)):
            for k in range(len(pax)):
                x, y, z = pax[i], yax[j], pax[k]
                same_grid[cnt][0] = x
                same_grid[cnt][1] = y
                same_grid[cnt][2] = z
                if min(Aonax[:,1])<= y <= max(Aonax[:,1]):
                    same_grid[cnt][3] += ai_x([x, y, z])
                    same_grid[cnt][4] += ai_y([x, y, z])
                    same_grid[cnt][5] += ai_z([x, y, z])
                    #print(same_grid[cnt][3:6])
                if min(Bonax[:,1])<= y <= max(Bonax[:,1]):
                    same_grid[cnt][3] += bi_x([x, y, z])
                    same_grid[cnt][4] += bi_y([x, y, z])
                    same_grid[cnt][5] += bi_z([x, y, z])
                    #print(same_grid[cnt][3:6])
                if min(Conax[:,1])<= y <= max(Conax[:,1]):
                    same_grid[cnt][3] += ci_x([x, y, z])
                    same_grid[cnt][4] += ci_y([x, y, z])
                    same_grid[cnt][5] += ci_z([x, y, z])
                    #print(same_grid[cnt][3:6])
                if min(Donax[:,1])<= y <= max(Donax[:,1]):
                    same_grid[cnt][3] += di_x([x, y, z])
                    same_grid[cnt][4] += di_y([x, y, z])
                    same_grid[cnt][5] += di_z([x, y, z])
                    #print(same_grid[cnt][3:6])
                cnt += 1

    #print(same_grid)
    stc_x, stc_y, stc_z = generate_field_interp(same_grid, response_curve = True)

    return stc_x, stc_y, stc_z