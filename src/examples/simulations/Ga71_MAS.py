# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 16:47:00 2023

@author: Marcos
"""
from mrsimulator import Simulator, SpinSystem, Site, SpectralDimension
from mrsimulator.spin_system.isotope import Isotope
from mrsimulator.method.lib import BlochDecaySpectrum
from mrsimulator.spin_system.tensors import SymmetricTensor
from mrsimulator import signal_processor as sp
from mrsimulator.models import CzjzekDistribution
from mrsimulator.utils.collection import single_site_system_generator
import matplotlib.pyplot as plt
from ssnmr_proc import bruker as br
import numpy as np
from scipy.stats import multivariate_normal

## get experimental details

'''=========== Just to get git-repository main directory ================='''
import os
import git
def get_git_root(path):  
    git_repo = git.Repo(path, search_parent_directories=True)
    git_root = git_repo.git.rev_parse("--show-toplevel")
    return git_root
''' ======================================================================'''

'''================ Get Czjzek distribution ============================= '''
# Load experimental data into class ProcData
data_dir = get_git_root(os.getcwd()) + r'\data\raw\GaGeO_glasses\1.3mm_Younes_21_10_28_71Ga_50-20-30\101\pdata\1'
pd = br.ProcData(data_dir) # Create pd as object from ProcData class
pd.normalize()
pd.reffrq

# Create czjzek distribution object
cz =  CzjzekDistribution(sigma=2.0)

# Define ranges for Cq and eta
Cq_range = np.linspace(0,26,40)
eta_range = np.linspace(0,1,20)

#get Czjzek probability distribution function (pdf)
Cq_dist, eta_dist, cz_amp = cz.pdf(pos=[Cq_range,eta_range])

fig, ax = plt.subplots(figsize=(5, 5))
ax.contourf(Cq_dist, eta_dist, cz_amp, levels=10)
plt.xlabel("$C_q$ / MHz")
plt.ylabel("$\eta$")
plt.tight_layout()
plt.show()

''' ======================================================================'''

''' =================== Simulation step ================================='''
systems = single_site_system_generator(isotope="71Ga", 
                                       quadrupolar={
                                           "Cq": Cq_dist * 1e6, 
                                           "eta": eta_dist}, 
                                       abundance= cz_amp
)

MAS = BlochDecaySpectrum(
    channels=["71Ga"],
    rotor_frequency=60000,  # in Hz
    spectral_dimensions=[
        SpectralDimension(spectral_width=pd.udic[0]['sw'], reference_offset=pd.udic[0]['car'])  # values in Hz
    ],
)

sim = Simulator()
sim.spin_systems = systems  # add the spin systems
sim.methods = [MAS]
sim.run()

sim_spc = sim.methods[0].simulation.real.y[0]
sim_spc /= max(sim_spc)
sim_ppm_scale = sim.methods[0].simulation.x[0].coordinates #tenho que extrair o eixo x

fig, ax = plt.subplots(figsize=(4.25, 3.0))
ax.plot(sim_ppm_scale,sim_spc, color="black", linewidth=1)
ax.invert_xaxis()
plt.tight_layout()
plt.show()
