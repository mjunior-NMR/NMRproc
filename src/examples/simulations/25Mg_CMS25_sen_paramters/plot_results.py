# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 11:44:42 2024

@author: mcabe
"""

import csdmpy as cp
import matplotlib.pyplot as plt
from ssnmr_proc import nmrplot

samples = {'CMS25': r'CaMgSi$_2$O$_6$',
            'Na4': r'Na$_4$MgSi$_3$O$_9$',
            'Na2': r'Na$_2$MgSi$_2$O$_6$',
            'K4': r'K$_4$MgSi$_3$O$_9$',
            'K2': r'K2MgSi2O6'}

fig, axs = plt.subplots(1,5, num = 1, subplot_kw={"projection": "csdm"})

spc = list(range(5))

for i,key in enumerate(list(samples.keys())):
    spc[i] = cp.load('simulation_results/' + key + '_ext_fitdata.csdf')