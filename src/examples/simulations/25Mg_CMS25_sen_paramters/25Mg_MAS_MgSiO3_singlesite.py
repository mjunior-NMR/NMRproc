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
from mrsimulator.models import ExtCzjzekDistribution
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
data_dir = get_git_root(os.getcwd()) + r'\data\raw\MgO_glasses\25Mg_CMS25_180W\pdata\1'
base_corrected = np.loadtxt(
    get_git_root(os.getcwd()) + 
    r'\data\raw\MgO_glasses\25Mg_CMS25_180W\pdata\1\25Mg_CMS25_baseline.txt')


pd = br.ProcData(data_dir) # Create pd as object from ProcData class
pd.data = base_corrected[:,1]+1j*base_corrected[:,2]
pd.flip()

pd.normalize()

exp_spc = pd.to_csdm()







# Define ranges for Cq and eta
Cq_range = np.linspace(0,25,100)
eta_range = np.linspace(0,1,10)


diso1 = 16; Cq1 = 7.8; sigma1 = 1.06



# Create czjzek distribution object for site 1
quad_tensor = {"Cq": Cq1, "eta": 0.01}  # Cq assumed in MHz
model_quad = ExtCzjzekDistribution(quad_tensor, eps= sigma1)
Cq_dist, eta_dist, cz_amp = model_quad.pdf(pos=[Cq_range, eta_range])

# get average Cq
cz_sum_over_eta = np.sum(cz_amp,axis=0)
# plt.plot(Cq_range,cz_sum_over_eta)
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
        SpectralDimension(spectral_width=pd.udic[0]['sw'], reference_offset=pd.udic[0]['car'])  # values in Hz
    ],
    magnetic_flux_density = 14.1
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
ax.plot(exp_spc.real,color = 'k')
ax.plot(sim_spc.real*0.8, color = 'g')

#Plot opptions
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_yticks([])
plt.xlim([2000,-2000])
plt.ylim([-0.05,1.2])
plt.xlabel("$^{25}$Mg-$\delta$ / ppm", fontsize = 14)
plt.ylabel("")
plt.savefig(r'Sim_lowfield_singlesite.png', dpi=300)
# plt.show()


#Save simulation to XRI and csdf formats.

sim_spc.y[0].components.astype('complex64')



sim_spc.save(r'simulation_Sen_lowfield.csdf') 




# R = sim_spc.y[0].components.real.transpose()
# I = sim_spc.y[0].components.imag.transpose()
# X = sim_spc.x[0].coordinates.value.transpose().reshape(sim_spc.size,1)
# # X = X*pd.reffrq


# # plt.plot(X,R)

# np.savetxt(r'simulation_Sen_lowfield.dat',np.concatenate((X,R,I),axis=1))




print("--- %s seconds ---" % (time.time() - start_time))
