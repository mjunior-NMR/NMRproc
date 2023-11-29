# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 16:01:13 2023

@author: Marcos
"""

import nmrglue as ng
from nmrglue import proc_base as pd
import matplotlib.pyplot as plt
import numpy as np


data_dir = r'D:\Marcos\IFSC\Orientacoes\Wellington\Dados\Agilent\onepul_FE593-1_6khz.fid'


# read in the Agilent data
dic, fid = ng.varian.read(data_dir)
procpar = ng.varian.read_procpar(data_dir+r'\procpar')

# Set the spectral parameters.
sw = float(procpar['sw']['values'][0])
np = float(procpar['np']['values'][0])
rfl = float(procpar['rfl']['values'][0])
rfp = float(procpar['rfp']['values'][0])
reffrq = float(procpar['reffrq']['values'][0])
sfrq = float(procpar['sfrq']['values'][0])
SR = sw/2-rfl+rfp; # To get the reference position

udic = {'ndim': 1,
        0: {'sw': float(procpar['sw']['values'][0]),
            'complex': True,
            'obs': reffrq,
            'car': SR, # reference position (spectral reference, as from TopSpin)
            'size': float(procpar['np']['values'][0])/2.0,
            'label': 'FID',
            'encoding': 'direct',
            'time': True,
            'freq': False}}

# Processing
lb = float(procpar['lb']['values'][0])/sw
fn = float(procpar['fn']['values'][0])
lsfid = int(procpar['lsfid']['values'][0])

proc_fid = pd.zf(fid,int(procpar['fn']['values'][0]))
if lsfid != 0:
    proc_fid = pd.ls(proc_fid,lsfid)
proc_fid = pd.em(proc_fid,lb)

udic[0]['size'] = proc_fid.size
spec = pd.fft(proc_fid)
spec = pd.ps(spec,p0 = float(procpar['rp']['values'][0])-127,p1 = float(procpar['lp']['values'][0]))

uc = ng.fileiobase.uc_from_udic(udic) 
ppm_scale = uc.ppm_scale() # ppm axis
hz_scale = uc.hz_scale()


plt.plot(hz_scale,spec.real/spec.real.max())
# plt.xlim([5,-3])


