# -*- coding: utf-8 -*-
"""
Created on Sun Sep 10 12:27:56 2023

Plot multiple pre-processe Bruker spectra

@author: Marcos de Oliveira Jr.
"""

import ssnmr_processing.bruker as br
import matplotlib.pyplot as plt
# from itertools import zip

#%% set paths and sample names and load data into spc dictionary
mIm = {'Zy2-m': r'mIm, acet. sin agitación',
       'Zy3-m': r'mIm, acet. con agitación',
       'Zx3-m': r'mIm, metanol, con agitación',
       'Zx2-m': r'mIm, metanol, sin agitación'}

mIm_PhIm = {'Y3-m':  r'mIm+PhIm, acet., con agitación', 
            'Y2-m':  r'mIm+PhIm, acet., sin agitación',
            'Y1-h':  r'mIm+PhIm, acet., solvotermal'}

mainpath = (r'd:/marcos/ifsc/dados/bruker/Pablo_Zifs/')
# mainpath = (r'G:/Outros computadores/PC-IFSC/IFSC/Dados/Bruker/Pablo_ZIFs/')


fullpath = dict() # Dictionary with the files path
spcs = dict() # Dictionary containing proc-spectra for each sample.

for key in mIm_PhIm.keys():
    fullpath[key] = mainpath + '2.5mm_' + key + r'\3\pdata\1'
    print('loading from: ',fullpath[key])
    spcs[key] = br.ProcData(fullpath[key])
    

#%% Normalize and plot the data

fig, axs = plt.subplots(len(spcs),1,sharex=True,sharey=True,figsize=(6, 5))
fig.subplots_adjust(hspace=-0.01)

for i,spc in enumerate(spcs.values()):
     spc.normalize(method = 'area', region = (160,5))
     spc.plot(axis = axs[i], nucleus = '13C')
     
axs[0].set_xlim([160,5])
