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
import csdmpy as cd
import matplotlib.pyplot as plt
import os
import git
from mrsimulator import signal_processor as sp


'''=========== Just to get git-repository main directory ================='''
def get_git_root(path):  
    git_repo = git.Repo(path, search_parent_directories=True)
    git_root = git_repo.git.rev_parse("--show-toplevel")
    return git_root
''' ======================================================================'''
# from matplotlib.ticker import (MultipleLocator)
# import numpy as np

#%% set data_dir and create ProcData class objet

root_dir = get_git_root(os.getcwd())
data_dir = root_dir + r'\data\raw\C13_spectrum\3\pdata\1'



spc = br.ProcData(data_dir)  #Atribui a classe bruker.ProcData à variável spc
spc.normalize() #Normaliza o espectro


#%% Convert to csdm object
csdm_spc = spc.to_csdm()

csdm_spc.save(r'my_file.csdf')



#%% Save xy ascii data

plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.plot(csdm_spc, color="black", linewidth=0.5, label="Experiment")
ax.set_xlim(190, 110)
plt.grid()
plt.tight_layout()
plt.show()



