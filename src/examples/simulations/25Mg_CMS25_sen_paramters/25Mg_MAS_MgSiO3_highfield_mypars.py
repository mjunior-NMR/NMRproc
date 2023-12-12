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
import nmrglue as ng
import csdmpy as cp
## get experimental details

'''=========== Just to get git-repository main directory ================='''
import os
import git
def get_git_root(path):  
    git_repo = git.Repo(path, search_parent_directories=True)
    git_root = git_repo.git.rev_parse("--show-toplevel")
    return git_root
''' ======================================================================'''


# Load experimental data into class ProcData
data_dir = get_git_root(os.getcwd()) + r'src/examples/simulations/25Mg_CMS25_sen_paramters/simulation_Sen_lowfield.csdf'

exp_spc = cp.load(data_dir)
exp_spc.x[0].to("ppm", "nmr_frequency_ratio")
exp_spc /= exp_spc.real.max()

# simproc = sp.SignalProcessor(
# operations=[
# sp.baseline.Polynomial(polynomial_dictionary = {'c0':0.05, 'c2':-4e-9})
# ]
# )

# exp_spc = simproc.apply_operations(exp_spc.copy())




# Define ranges for Cq and eta
# Define ranges for Cq and eta
Cq_range = np.linspace(0,25,100)
eta_range = np.linspace(0,1,10)


diso1 = 13; sigma1 = 4.242/2


# get Czjzek probability distribution function (pdf)
cz =  CzjzekDistribution(sigma=sigma1)
Cq_dist, eta_dist, cz_amp = cz.pdf(pos=[Cq_range,eta_range])

# get average Cq
cz_sum_over_eta = np.sum(cz_amp,axis=0)
Cq_av = np.sum(cz_sum_over_eta*Cq_range)/np.sum(cz_sum_over_eta)
print("the average Cq for site_1 is  %s MHz" % (round(Cq_av*10)/10))

 # plt.close('all')

# fig, ax = plt.subplots(figsize=(5, 6),num=1)
# ax.contourf(Cq_dist, eta_dist, cz_amp, levels=10)
# plt.xlabel("$C_q$ / MHz")
# plt.ylabel("$\eta$")



############# Simulation step #################

#%% Simulation
system_1 = single_site_system_generator(isotope="25Mg",
                                       isotropic_chemical_shift = diso1,
                                       quadrupolar={
                                           "Cq": Cq_dist * 1e6, 
                                           "eta": eta_dist}, 
                                       abundance= cz_amp                                       
                                       )



MAS = BlochDecayCTSpectrum(
    channels=["25Mg"],
    rotor_frequency=20000,  # in Hz
    spectral_dimensions=[
        SpectralDimension(spectral_width=500000, reference_offset=0)  # values in Hz
    ],
    magnetic_flux_density = 35.2
)

sim = Simulator()
sim.config.number_of_sidebands = 64 # eight sidebands are sufficient for this example
sim.config.integration_volume = "octant"
sim.config.integration_density = 20 # (see https://mrsimulator.readthedocs.io/en/stable/user_guide/simulator/simulator.html?highlight=sphere#integration-volume)
sim.spin_systems = system_1  # add the spin systems
sim.methods = [MAS]
sim.run()

sim_spc = sim.methods[0].simulation
sim_spc.y[0].components = sim_spc.y[0].components/abs(np.trapz(sim_spc.y[0].components.real))


#%% Post processing and plot

simproc = sp.SignalProcessor(
operations=[
sp.IFFT(),
sp.apodization.Gaussian(FWHM="1000 Hz"),
sp.FFT(),
]
)

sim_spc = simproc.apply_operations(sim_spc)

sim_spc /= sim_spc.max()



fig, ax = plt.subplots(figsize=(6, 5), subplot_kw={"projection": "csdm"}, num=2)
ax.plot(exp_spc.real, color = 'k', label = 'experimental')
ax.plot(sim_spc.real, color = 'g', linestyle = (0,(3,1,1,1)), label = 'Sum')
plt.legend(loc='center right')

#Plot opptions
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_yticks([])
plt.xlim([300,-400])
plt.ylim([-0.05,1.2])
plt.xlabel("$^{25}$Mg-$\delta$ / ppm", fontsize = 14)
plt.ylabel("")
plt.savefig(r'Sim_highfield_mypars.png', dpi=300)
# plt.show()


