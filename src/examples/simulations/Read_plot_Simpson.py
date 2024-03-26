# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 08:26:17 2023

@author: Marcos
"""

import os
import nmrglue as ng
from matplotlib import pyplot as plt
import numpy as np


# os.system('simpson input_1H_13C.in') # Run simulation
 
dic,spec = ng.simpson.read('simpson_output.spe') # read simpson simulation

b = spec.real.reshape(8192,)

sw = dic['SW']
num_points = dic['NP']

num_points = 8192
sw = 20000

freq = np.linspace(-sw/2,sw/2,num_points).reshape(num_points,)
# freq
a = np.arange(0,8192)
c = a**2

plt.plot(freq,b)


# ax.set_yticks([])
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# ax.spines['left'].set_visible(False)