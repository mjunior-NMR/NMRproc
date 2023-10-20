# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 15:00:14 2023

@author: Marcos de Oliveira Jr.
"""

# import mrsimulator as mrs
# from mrsimulator import Simulator, SpinSystem, Site, Coupling
# from mrsimulator.spin_system.isotope import Isotope
# from mrsimulator.method.lib import BlochDecaySpectrum
# from mrsimulator.spin_system.tensors import SymmetricTensor
# import csdmpy as cp


import ssnmr_proc.bruker as br
# import nmrglue as ng
import matplotlib.pyplot as plt
import os
import git

'''=========== Just to get git-repository main directory ================='''
def get_git_root(path):  
    git_repo = git.Repo(path, search_parent_directories=True)
    git_root = git_repo.git.rev_parse("--show-toplevel")
    return git_root
''' ======================================================================'''
# from matplotlib.ticker import (MultipleLocator)
# import numpy as np

#%% set data_dir and create ProcData class objet


data_dir = get_git_root(os.getcwd()) + r'\data\raw\C13_spectrum\3\pdata\1'

# data_dir = r'D:\Marcos\IFSC\Dados\Bruker\camila\2.5mm_Cdbtc_1hBM_2022-11-04\3\pdata\1'

spc = br.ProcData(data_dir)  #Atribui a classe bruker.ProcData à variável spc
spc.normalize() #Normaliza o espectro
# spc.normalize(method = 'area', region = (190,110)) #Normaliza o espectro

#%% Plot single spectrum using plot method of class ProcData

fig, ax = plt.subplots(sharex=True,sharey=True,figsize=(6, 5))

spc.plot(axis = ax, nucleus = '13C',major_ticks_space=10,minor_ticks_space=2) 
ax.set_xlim([190,110])
ax.set_ylim([-0.05,1])

#%% Save data in multiple formats
spc.save(r'spectrum.dat') # Save spectrum in ascii file in xy format.
spc.save(r'spectrum.txt', file_type='dmfit') # Save a dmfit readable spectrum file
spc.save(r'spectrum.csdf', file_type='csdf') #save a csdm object in json format.
spc.save(r'procdata.pkl', file_type = 'pickle') # Save a binary file containing the spc object of class ProcData



# converter = ng.fileio.convert.converter()

# converter.from_bruker(dic=spc.dic, data=spc.rdata, udic=None, remove_digital_filter=False)


# spec = converter.to_csdm()

# spec = spec.real

# spec.x[0].to("ppm", "nmr_frequency_ratio")
# plt.figure(figsize=(4.25, 3.0))
# ax = plt.subplot(projection="csdm")
# ax.plot(spec, color="black", linewidth=0.5, label="Experiment")
# ax.set_xlim(150, -150)
# plt.grid()
# plt.tight_layout()
# plt.show()



