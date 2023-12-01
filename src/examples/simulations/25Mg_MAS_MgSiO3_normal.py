# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 16:47:00 2023

@author: Marcos
"""
import time
start_time = time.time()
from mrsimulator import Simulator, SpinSystem, Site, SpectralDimension
from mrsimulator.spin_system.isotope import Isotope
from mrsimulator.method.lib import BlochDecayCTSpectrum,BlochDecaySpectrum
from mrsimulator.spin_system.tensors import SymmetricTensor
from mrsimulator import signal_processor as sp
from mrsimulator.models import CzjzekDistribution
from mrsimulator.utils.collection import single_site_system_generator
import matplotlib.pyplot as plt
from ssnmr_proc import bruker as br
import numpy as np
from scipy.stats import multivariate_normal,moment

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
data_dir = get_git_root(os.getcwd()) + r'\data\raw\MgO_glasses\25Mg_MgSiO3_180W\pdata\1'
pd = br.ProcData(data_dir) # Create pd as object from ProcData class
pd.normalize()
exp_spc = pd.to_csdm()

# Create czjzek distribution object
cz =  CzjzekDistribution(sigma=2.22)

# Define ranges for Cq and eta
Cq_range = np.linspace(0,18,50)
eta_range = np.linspace(0,1,10)

# get Czjzek probability distribution function (pdf)
Cq_dist, eta_dist, cz_amp = cz.pdf(pos=[Cq_range,eta_range])


# get average Cq
cz_sum_over_eta = np.sum(cz_amp,axis=0)
# plt.plot(Cq_range,cz_sum_over_eta)
Cq_av = np.sum(cz_sum_over_eta*Cq_range)/np.sum(cz_sum_over_eta)
print(Cq_av)


# plt.close('all')

fig, ax = plt.subplots(figsize=(5, 6),num=1)
ax.contourf(Cq_dist, eta_dist, cz_amp, levels=10)
plt.xlabel("$C_q$ / MHz")
plt.ylabel("$\eta$")


''' ======================================================================'''

''' =================== Simulation step ================================='''

#%% Simulation
systems = single_site_system_generator(isotope="25Mg",
                                       isotropic_chemical_shift = 13,
                                       quadrupolar={
                                           "Cq": Cq_dist * 1e6, 
                                           "eta": eta_dist}, 
                                       abundance= cz_amp,
                                       shielding_symmetric = {'zeta': 150, 'eta': 0}
                                       )


MAS = BlochDecayCTSpectrum(
    channels=["25Mg"],
    rotor_frequency=20000,  # in Hz
    spectral_dimensions=[
        SpectralDimension(spectral_width=pd.udic[0]['sw'], reference_offset=pd.udic[0]['car'])  # values in Hz
    ],
    magnetic_flux_density = 14.1
)

sim = Simulator()
sim.config.number_of_sidebands = 16 # eight sidebands are sufficient for this example
sim.config.integration_volume = "octant"
sim.config.integration_density = 20 # (see https://mrsimulator.readthedocs.io/en/stable/user_guide/simulator/simulator.html?highlight=sphere#integration-volume)
sim.spin_systems = systems  # add the spin systems
sim.methods = [MAS]
sim.run()

sim_spc = sim.methods[0].simulation



#%% Post processing and plot

simproc = sp.SignalProcessor(
operations=[
sp.IFFT(),
sp.apodization.Gaussian(FWHM="5000 Hz"),
sp.FFT(),
]
)

sim_spc = simproc.apply_operations(sim_spc)
sim_spc /= sim_spc.max()



fig, ax = plt.subplots(figsize=(6, 5), subplot_kw={"projection": "csdm"}, num=2)
ax.plot(exp_spc.real)
ax.plot(sim_spc.real, color = 'k')
plt.xlim([3000,-4000])
plt.show()

print("--- %s seconds ---" % (time.time() - start_time))
