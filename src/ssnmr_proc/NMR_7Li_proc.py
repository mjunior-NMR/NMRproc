# -*- coding: utf-8 -*-
"""
Created on Sun Sep 10 12:27:56 2023

Plot multiple pre-processe Bruker spectra

@author: Marcos de Oliveira Jr.
"""

import ssnmr_proc.bruker as br
import matplotlib.pyplot as plt


# from itertools import zip

#%% set paths and sample names and load data into spc dictionary
samples = {'FE593-1': r'Reference',
              'FE723-4': r'50 ALD cycles',
              'FE723-5': r'Reference',
              'FE723-6': r'400 ALD cycles'}


data_dir = r'D:\Marcos\IFSC\Orientacoes\Wellington\Dados\Agilent'


# mainpath = (r'G:/Outros computadores/PC-IFSC/IFSC/Dados/Bruker/Pablo_ZIFs/')


fullpath = dict() # Dictionary with the files path
spcs_MAS = dict() # Dictionary containing proc-spectra for each sample.
spcs_static = dict() # Dictionary containing proc-spectra for each sample.

#Load MAS data into spc_MAS
for key in samples.keys():
    fullpath[key] = data_dir + '\onepul_' + key + r'_6khz.fid'
    print('loading from: ',fullpath[key])
    spcs_MAS[key] = br.ProcData(fullpath[key])
    
#Load static data into spc_static
for key in samples.keys():
    fullpath[key] = data_dir + 'onepul_' + key + r'_estatico.fid'
    print('loading from: ',fullpath[key])
    spcs_static[key] = br.ProcData(fullpath[key])
    
# Save ascii data


fig, axs = plt.subplots(len(spcs_MAS),1,sharex=True,sharey=True,figsize=(6, 5))
fig.subplots_adjust(hspace=-0.01)

for i,(sample,spc) in enumerate(spcs_MAS.items()):
     spc.normalize()
     spc.save
     spc.plot(axis = axs[i], nucleus = '13C')
     
axs[0].set_xlim([160,5])
