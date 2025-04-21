import numpy as np 
import matplotlib.pyplot as plt 
import os

from rgi_field import generate_field_interp

parent_dir = 'STC_3DMAPS_2025/response_curves'

coilset = 'SU'
        ## A, B, C, D
currents = (1e-3)*np.array([10, 30, 30, 100])

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
            Afield[:,1] = response[:,1] + start
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

#Aonax = Afield[(Afield[:,0]==0) & (Afield[:,2]==0)]
Bonax = Bfield[(Bfield[:,0]==0) & (Bfield[:,2]==0)]
Conax = Cfield[(Cfield[:,0]==0) & (Cfield[:,2]==0)]
Donax = Dfield[(Dfield[:,0]==0) & (Dfield[:,2]==0)]

#ai_x, ai_y, ai_z = generate_field_interp(Afield, response_curve = True)
bi_x, bi_y, bi_z = generate_field_interp(Bfield, response_curve = True)
ci_x, ci_y, ci_z = generate_field_interp(Cfield, response_curve = True)
di_x, di_y, di_z = generate_field_interp(Dfield, response_curve = True)

N2 = 31
cnt = 0
same_grid = np.zeros((N2**3, 6))
yax = np.linspace(0, 0.8, N2)
pax = np.linspace(-0.01, 0.01, N2)
for i in range(len(pax)):
    for j in range(len(yax)):
        for k in range(len(pax)):
            x, y, z = pax[i], yax[j], pax[k]
            same_grid[cnt][0] = x
            same_grid[cnt][1] = y
            same_grid[cnt][2] = z
            if min(Bonax[:,1])<= y <= max(Bonax[:,1]):
                same_grid[cnt][3] += bi_x([x, y, z])[0]
                same_grid[cnt][4] += bi_y([x, y, z])[0]
                same_grid[cnt][5] += bi_z([x, y, z])[0]
            if min(Conax[:,1])<= y <= max(Conax[:,1]):
                same_grid[cnt][3] += ci_x([x, y, z])[0]
                same_grid[cnt][4] += ci_y([x, y, z])[0]
                same_grid[cnt][5] += ci_z([x, y, z])[0]
            if min(Donax[:,1])<= y <= max(Donax[:,1]):
                same_grid[cnt][3] += di_x([x, y, z])[0]
                same_grid[cnt][4] += di_y([x, y, z])[0]
                same_grid[cnt][5] += di_z([x, y, z])[0]
            cnt += 1

#Bgrid = same_grid[(same_grid[:,1]<=max(Bonax[:,1])) & (same_grid[:,1]>=min(Bonax[:,1]))][:,:3]
#Bcont = bi_z(Bgrid)
#Binterp = np.zeros((len(Bgrid), 6))
#Binterp[:,:3] = Bgrid
#Binterp[:,5] = Bcont
#bf_z = Binterp[(Binterp[:,0]==0) & (Binterp[:,2]==0)]

#Cgrid = same_grid[(same_grid[:,1]<=max(Conax[:,1])) & (same_grid[:,1]>=min(Conax[:,1]))][:,:3]
#Ccont = ci_z(Cgrid)
#Cinterp = np.zeros((len(Cgrid), 6))
#Cinterp[:,:3] = Cgrid
#Cinterp[:,5] = Ccont
#cf_z = Cinterp[(Cinterp[:,0]==0) & (Cinterp[:,2]==0)]

#Dgrid = same_grid[(same_grid[:,1]<=max(Donax[:,1])) & (same_grid[:,1]>=min(Donax[:,1]))][:,:3]
#Dcont = di_z(Dgrid)
#Dinterp = np.zeros((len(Dgrid), 6))
#Dinterp[:,:3] = Dgrid
#Dinterp[:,5] = Dcont
#df_z = Dinterp[(Dinterp[:,0]==0) & (Dinterp[:,2]==0)]

combo_onax = same_grid[(same_grid[:,0]==0) & (same_grid[:,2]==0)]

fig0, ax0 = plt.subplots(3, 1, sharex=True)
ax0[0].scatter(Bonax[:,1], Bonax[:,3]*1e6, color = 'k')
ax0[0].scatter(Conax[:,1], Conax[:,3]*1e6, color = 'k')
ax0[0].scatter(Donax[:,1], Donax[:,3]*1e6, color = 'k')
ax0[0].plot(combo_onax[:,1], combo_onax[:,3]*1e6, color='r')

ax0[1].scatter(Bonax[:,1], Bonax[:,4]*1e6, color = 'k')
ax0[1].scatter(Conax[:,1], Conax[:,4]*1e6, color = 'k')
ax0[1].scatter(Donax[:,1], Donax[:,4]*1e6, color = 'k')
ax0[1].plot(combo_onax[:,1], combo_onax[:,4]*1e6, color='r')

#plt.scatter(Aonax[:,1], Aonax[:,5])
ax0[2].scatter(Bonax[:,1], Bonax[:,5]*1e6, color = 'k')
ax0[2].scatter(Conax[:,1], Conax[:,5]*1e6, color = 'k')
ax0[2].scatter(Donax[:,1], Donax[:,5]*1e6, color = 'k')
ax0[2].plot(combo_onax[:,1], combo_onax[:,5]*1e6, color = 'r')



plt.show()