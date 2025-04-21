import numpy as np 
import matplotlib.pyplot as plt 
import os 

parent_dir = 'STC_3DMAPS_2025'

for i in os.listdir(parent_dir):
    if i[1:] != 'coils':
        continue
    else:
        print(i)
        for j in os.listdir(parent_dir+'/'+i):
            for k in os.listdir(parent_dir+'/'+i):
                if j==k:
                    continue
                else:
                    sega = j.split('_')
                    segb = k.split('_')
                    if sega[0] == segb[0]:
                        cur1 = float(sega[1].strip('mA'))
                        cur2 = float(segb[1].strip('mA'))
                        if cur1 > cur2:
                            odata_hi = np.loadtxt(parent_dir+'/'+i+'/'+j)
                            odata_lo = np.loadtxt(parent_dir+'/'+i+'/'+k)
                        else:
                            odata_hi = np.loadtxt(parent_dir+'/'+i+'/'+k)
                            odata_lo = np.loadtxt(parent_dir+'/'+i+'/'+j)
                        data_hi = odata_hi[:(len(odata_hi)-1)]
                        data_lo = odata_lo[:(len(odata_lo)-1)]
                    
                        dcur = (data_hi[:,3] - data_lo[:,3])
                        response_curves = np.zeros((len(data_hi),3))
                        response_curves[:,0] = (data_hi[:,4] - data_lo[:,4])/dcur
                        response_curves[:,1] = (data_hi[:,5] - data_lo[:,5])/dcur
                        response_curves[:,2] = (data_hi[:,6] - data_lo[:,6])/dcur
                        print(data_hi[:,0])
                        with open(parent_dir+'/response_curves/'+sega[0]+'_response.txt', 'w+') as newf:
                            for p in range(len(data_hi)):
                                newf.write(f'{data_hi[p,0]} {data_hi[p,1]} {data_hi[p,2]} {response_curves[p,0]} {response_curves[p,1]} {response_curves[p,2]}\n')