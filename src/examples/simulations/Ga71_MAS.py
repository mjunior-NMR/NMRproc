# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 16:47:00 2023

@author: Marcos
"""
import time
start_time = time.time()
from mrsimulator import Simulator, SpinSystem, Site, SpectralDimension
from mrsimulator.spin_system.isotope import Isotope
from mrsimulator.method.lib import BlochDecayCTSpectrum
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
pd.normalize(method = '01')
exp_spc = pd.to_csdm()

# Create czjzek distribution object
cz =  CzjzekDistribution(sigma=3.9)

# Define ranges for Cq and eta
Cq_range = np.linspace(0,30,100)
eta_range = np.linspace(0,1,20)

# get Czjzek probability distribution function (pdf)
Cq_dist, eta_dist, cz_amp = cz.pdf(pos=[Cq_range,eta_range])

# plt.close('all')
'''
fig, ax = plt.subplots(figsize=(5, 5))
ax.contourf(Cq_dist, eta_dist, cz_amp, levels=10)
plt.xlabel("$C_q$ / MHz")
plt.ylabel("$\eta$")
plt.tight_layout()
plt.show()
'''
''' ======================================================================'''

''' =================== Simulation step ================================='''

#%% Simulation
systems = single_site_system_generator(isotope="71Ga",
                                       isotropic_chemical_shift = 219.6,
                                       quadrupolar={
                                           "Cq": Cq_dist * 1e6, 
                                           "eta": eta_dist}, 
                                       abundance= cz_amp,
                                       shielding_symmetric = {'zeta': 100, 'eta': 1}
                                       )


MAS = BlochDecayCTSpectrum(
    channels=["71Ga"],
    rotor_frequency=60000,  # in Hz
    spectral_dimensions=[
        SpectralDimension(spectral_width=pd.udic[0]['sw'], reference_offset=pd.udic[0]['car'])  # values in Hz
    ],
    magnetic_flux_density = 14.1
)

sim = Simulator()
sim.config.number_of_sidebands = 16 # eight sidebands are sufficient for this example
sim.spin_systems = systems  # add the spin systems
sim.methods = [MAS]
sim.run()

#%% Post processing
sim_spc = sim.methods[0].simulation

simproc = sp.SignalProcessor(
operations=[
sp.IFFT(),
sp.apodization.Gaussian(FWHM="3000 Hz"),
sp.FFT(),
]
)

sim_spc = simproc.apply_operations(sim_spc)
sim_spc /= sim_spc.max()

#%% Plot the spectra
fig, ax = plt.subplots(figsize=(6, 3), subplot_kw={"projection": "csdm"})
ax.plot(exp_spc.real)
ax.plot(sim_spc.real)
plt.tight_layout()
# plt.xlim([1500,-2000])
plt.show()

print("--- %s seconds ---" % (time.time() - start_time))
