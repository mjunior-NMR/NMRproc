# -*- coding: utf-8 -*-
"""
Created on Fri May 24 11:33:11 2024

Read Json of Agilent data processed on ssNake

@author: Marcos
"""

import os
import json
import csdmpy as cp
import numpy as np
import matplotlib.pyplot as plt
from ssnmr_proc.nmrplot import nmrstyle

data_dir = r'D:\Marcos\IFSC\Dados\Agilent\Renato\A3campinas_31P_BPO4_10kHz.fid'

file_dir = data_dir + r'\\' + os.path.split(data_dir)[1].replace('.fid','.json')

f = f = open(file_dir)

spec = json.load(f)
data = np.array(spec['dataReal'][0])+1j*np.array(spec['dataImag'][0])
np = data.size
Hz_scale = spec['xaxArray'][0]

dv = cp.as_dependent_variable(data, unit="")
dim = cp.LinearDimension(
    count = np, 
    origin_offset = f'{spec["freq"][0]} Hz',  
    #coordinates_offset = f'{spec["metaData"]["Offset [Hz]"]} Hz',
    coordinates_offset = f'{spec["freq"][0]-spec["ref"][0]} Hz',
    increment = f'{Hz_scale[1]-Hz_scale[0]} Hz',
    complex_fft=True,
    label="Frequency",
    reciprocal={'quantity_name': 'time'}
    )
csdm_spec = cp.CSDM(dependent_variables=[dv], dimensions=[dim])
csdm_spec.dimensions[0].to("ppm", "nmr_frequency_ratio") 

x = csdm_spec.x[0].coordinates
y = csdm_spec.y[0].components[0]

plt.figure(figsize=(5, 3.5), num = 1)
ax = plt.subplot(projection="csdm")
ax.plot(csdm_spec)
plt.tight_layout()
plt.show()
nmrstyle(ax)

# plt.plot(x,y)