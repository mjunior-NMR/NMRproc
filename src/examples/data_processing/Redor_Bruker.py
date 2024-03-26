# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 11:04:35 2024

@author: Marcos
"""

import nmrglue as ng
import matplotlib.pyplot as plt
import ssnmr_proc.bruker as br
import numpy as np

data_dir = r'D:\Marcos\IFSC\Dados\Bruker\2020\2.5mm_NaPO3F_24_03_20\4\pdata\1'


spc = br.ProcData(data_dir)

spc.normalize()

ax = plt.subplot()
ax.plot(spc.ppm_scale,spc.data.real[:,48])
ax.set_xlim([20,-40])



region = (-40,20)

nspec = spc.dic['acqu2s']['TD']

area = np.zeros(nspec)

for i in range(0,nspec):
    data = spc.data.real[:,i]
    p1 = spc.uc(str(region[0])+' ppm') #convert from ppm_scale to points p1 and p2
    p2 = spc.uc(str(region[1])+' ppm')
    reduced_rdata = data.real[min(p1,p2):max(p1,p2)]
    reduced_ppm_scale = spc.ppm_scale[min(p1,p2):max(p1,p2)]
    area[i] = abs(np.trapz(reduced_rdata,reduced_ppm_scale))

DS = np.zeros(int(nspec/2))
NTR = np.zeros(int(nspec/2))

k = 0

for i in range(0,nspec,2):
    DS[k] = (area[i+1]-area[i])/area[i+1]
    NTR[k] = 1/spc.dic['acqus']['CNST'][31]*1e3*i+2
    k = k+1



fig,ax = plt.subplots()
ax.scatter(NTR,DS)
    

a = list(range(0,nspec,2))








# plt.plot(x,data[10,:])


