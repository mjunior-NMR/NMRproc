# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 16:07:16 2024

@author: Marcos
"""


import matplotlib.pyplot as plt
from ssnmr_proc import nmrplot
from ssnmr_proc import agilent as ag
from ssnmr_proc import phasespec as ps
import nmrglue as ng
from pybaselines import Baseline
from matplotlib.widgets import Slider, Button




def getpars(procpar_dic):
    pars = dict()
    for key in procpar_dic.keys():
        pars[key]  = procpar_dic[key]['values']
    return pars

data_dir = r"D:\Marcos\IFSC\Dados\Agilent\NaPO3-GaF3\spe_19F_20GaF_35KHz.fid"
# data_dir = r"D:\Marcos\IFSC\Dados\Agilent\NaPO3-GaF3\spe_19F_40GaF_35KHz_2FIDs.fid"
# data_dir = r'D:\Marcos\IFSC\Python\NMRproc\data\raw\agilent2pipe_1d\agilent_1d'



dic, aux = ng.agilent.read(data_dir)
procpar = ng.fileio.agilent.read_procpar(data_dir + r'\procpar')
pars = getpars(procpar)
udic = ng.agilent.guess_udic(dic, aux)

if aux.ndim == 3:
    fid = aux[0,0,:]+aux[0,1,:]
else:
    fid = aux   



lb = float(pars['lb'][0])
fn = int(pars['fn'][0])
np = int(pars['np'][0])/2
sw = float(pars['sw'][0])
sfrq = float(pars['sfrq'][0])
tn = pars['tn'][0]
rfl = float(pars['rfl'][0])
rfp = float(pars['rfp'][0])
rp = float(pars['lp'][0])
lp = float(pars['rp'][0])
ls = int(pars['lsfid'][0])

SR = sw/2-rfl+rfp
reffrq = (sfrq*1e6-SR)*1e-6
car = (reffrq-sfrq)*1e6



# Set the spectral parameters.
udic = ng.agilent.guess_udic(dic, fid)
udic[0]['size']     = np
udic[0]['complex']  = True
udic[0]['encoding'] = 'direct'
udic[0]['sw']       = sw
udic[0]['obs']      = sfrq
udic[0]['car']      = 99.0*125.681
udic[0]['label']    = tn


proc_fid = ng.proc_base.em(fid, lb = lb/sw) # Window function
proc_fid = ng.proc_base.ls(proc_fid, pts = ls) # leftshift
spc = ng.proc_base.fft(proc_fid) # FFT
# spc = ng.proc_autophase.autops(spc, fn = 'acme', p1 = lp) # Automatic phase correction

p0,p1 = ps.manual(spc) # manual phase correction


baseline_fitter = Baseline() # Baseline corrector object from pybaselines

base, params_3 = baseline_fitter.mor(spc.real, half_window=30) # Morfological baseline correcttion, half_window takes the size of the broader beak.
# spc = ng.proc_autophase.manual_ps(spc)

# Build xaxis

fig, axs = plt.subplots(3,1, figsize=(6, 5), num = 2)

axs[0].plot(proc_fid)
axs[0].plot(proc_fid)
axs[1].plot(spc)
axs[1].plot(base,color = 'r')
axs[2].plot(spc-base)
plt.xlim(600,400)
nmrplot.nmrstyle(axs[1])
nmrplot.nmrstyle(axs[2])


