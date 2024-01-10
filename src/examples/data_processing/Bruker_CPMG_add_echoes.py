# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 15:03:56 2024

@author: Marcos
"""

import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
import ssnmr_proc.bruker as br
import csdmpy as cp
import nmrglue as ng
import ssnmr_proc.nmrplot as nmrplt

# read in the bruker formatted data
data_dir = r'D:\Marcos\IFSC\Dados\Bruker\2020\3.2mm_Philip_23_12_14_CMS25\3\pdata\1'
dic, data = ng.bruker.read(data_dir)

spec  = br.ProcData(data_dir)
spec.normalize()
spiklet = spec.to_csdm()

# remove the digital filter
data = ng.bruker.remove_digital_filter(dic, data)
data /= max(data)

dw = 1/dic['acqus']['SW_h']
num_points = np.size(data)
D = dic['acqus']['D']
P = [i * 1e-6 for i in dic['acqus']['P']]
CNST = dic['acqus']['CNST']
L = dic['acqus']['L']
aq = dic['acqus']['TD']/2*dw


acqu= D[6]
D[31]
echo1=D[31]*L[3]-P[3]/2-P[4]/2-dw-D[3]
echo2=L[3]*D[31]-P[4]/2
dead=echo2-acqu
echoes=(aq-acqu)/(2*acqu+2*dead+P[4])-0.5
rest=aq-echoes*(2*acqu+2*dead+P[4])-data


time = np.arange(num_points)*dw
# plt.plot(time,data.real)

acqu_time = round(acqu/dw)
dead_time = round(2*dead/dw)
pulse_time = round(P[4]/dw)

# cut = 0
# cut = 3*acqu_time+dead_time+pulse_time

FID = (round(echoes)+1)*[None] #pre-allocation


#%%

# First FID:
cut = acqu_time+dead_time+pulse_time
FID[0] = data[0:cut-dead_time-pulse_time]

for i in np.arange(round(echoes))+1:    
    cut += 2*acqu_time+dead_time+pulse_time
    aux = data[cut:cut+2*acqu_time]
    mid = round(np.size(aux)/2)    
    FID[i] = np.flip(aux[0:mid]) + aux[mid:]

new_TD = FID[1].size

sum_fid = np.zeros(new_TD,dtype = 'complex128')

for j in range(round(echoes)):        
    sum_fid += FID[j]

# plt.plot(time[0:new_TD],sum_fid)

#Zero fill data
z = np.zeros([8192,], dtype = 'complex128')
sum_fid = np.concatenate((sum_fid, z), axis=-1)

dv = cp.as_dependent_variable(sum_fid, unit="s")
dim = cp.LinearDimension(
    count = sum_fid.size, 
    increment = f'{dw} s',    
    quantity_name="Time",
    reciprocal={'origin_offset': f'{spec.udic[0]["obs"]*1e6} Hz',
                'coordinates_offset': f'{spec.udic[0]["car"]} Hz',
                'quantity_name': 'frequency'
                }
    )
csdm_FID = cp.CSDM(dependent_variables=[dv], dimensions=[dim])

envelope = csdm_FID.fft()

spiklet = spiklet/np.max(spiklet.y[0].components).real

envelope.dimensions[0].to("ppm", "nmr_frequency_ratio")
spiklet.dimensions[0].to("ppm", "nmr_frequency_ratio")



teta = 200
envelope *= np.exp(-2j*np.pi*teta/360)
envelope = envelope/np.max(envelope.y[0].components).real

# plot of the frequency domain data after FFT.
ax = plt.subplot(projection="csdm")
ax.plot(envelope.real)
ax.plot(spiklet.real)
nmrplt.nmrstyle(ax)

envelope.save(r'D:\Marcos\IFSC\Python\NMRproc\data\processed\25Mg NMR Phill\CMS25_cpmg_envelope.csdf')
