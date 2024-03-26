# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 08:26:17 2023

@author: Marcos
"""
import numpy as np
import nmrglue as ng
import matplotlib.pyplot as plt
import os



data_dir = r'6\pdata\1'
exp_dic,exp = ng.bruker.read_pdata(data_dir) #load experimental data.

os.system('simpson input.in') # Run simulation
 
dic,spec = ng.simpson.read('input.spe') # read simpson simulation

sim = spec.real.reshape(8192,)

sw = dic['SW']
num_points = dic['NP']

num_points = 8192
sw = 20000

#%%
freq = np.linspace(-sw/2,sw/2,num_points).reshape(num_points,)
ppm = (freq)/exp_dic['acqus']['O1']
# freq
a = np.arange(0,8192)
c = a**2

ax = plt.subplot()

plt.plot(freq,sim)
ax.set_yticks([])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
# plt.xlim(273,276)