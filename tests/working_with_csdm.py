# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 17:09:56 2023

@author: Marcos
"""

import numpy as np
import nmrglue as ng
import matplotlib.pyplot as plt 
from ssnmr_proc.bruker import ProcData
import csdmpy as cp

#%% Read data into ProcData object
data_dir = r"D:\Marcos\IFSC\Orientacoes\Silvio\Dados\Complexo_amarelo\106/pdata/1"
pd = ProcData(data_dir)
pd.normalize()


#%% Create csdm with nmrglue
converter = ng.fileio.convert.converter()
converter.from_universal(dic=pd.udic, data = np.flip(pd.data))
# Data is a spectrum, but nmrglue treat it as FID.
temp = converter.to_csdm()

#%% Fourrier transform temp to get reciprocal as main coordinates in spec
spec = temp.fft()
spec.y[0] = temp.y[0] # Add spectrum data into spec.y[0]

# Change units.
spec.dimensions[0].to("ppm", "nmr_frequency_ratio")



#%% Creating csdm by hand


dv = cp.as_dependent_variable(pd.data, unit="")
dim = cp.LinearDimension(
    count = pd.udic[0]['size'], 
    origin_offset = f'{pd.udic[0]["obs"]*1e6} Hz',  
    coordinates_offset = f'{pd.udic[0]["car"]} Hz',
    increment = f'{pd.hz_scale[1]-pd.hz_scale[0]} Hz',
    complex_fft=True,
    label="Frequency",
    reciprocal={'quantity_name': 'time', 'label': '31P'}
    )    
    
spec_1 = cp.CSDM(dependent_variables=[dv], dimensions=[dim])
spec_1.dimensions[0].to("ppm", "nmr_frequency_ratio")
    
#%% Plot
plt.figure(figsize=(5, 3.5))
ax = plt.subplot(projection = 'csdm')
ax.plot(spec.real, label = r'dirty solution', color = 'r')
ax1 = plt.subplot()
ax1.plot(pd.ppm_scale,pd.rdata, label = 'Original pdata', color = 'b')
ax = plt.subplot(projection = 'csdm')
ax.plot(spec_1.real, label = r'by hand', color = 'm', )
ax.set_xlim([-20,10])
plt.legend(loc="upper left")
