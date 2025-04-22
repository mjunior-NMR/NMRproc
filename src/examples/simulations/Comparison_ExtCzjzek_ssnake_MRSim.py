# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 11:51:43 2024

@author: Marcos
"""

from mrsimulator import Simulator, SpinSystem, Site, SpectralDimension
from mrsimulator.spin_system.isotope import Isotope
from mrsimulator.method.lib import BlochDecayCTSpectrum,BlochDecaySpectrum
from mrsimulator.spin_system.tensors import SymmetricTensor
from mrsimulator import signal_processor as sp
from mrsimulator.models import ExtCzjzekDistribution
from mrsimulator.utils.collection import single_site_system_generator
import matplotlib.pyplot as plt
from ssnmr_proc import bruker as br
from ssnmr_proc import nmrplot
from scipy.stats import multivariate_normal,moment
import csdmpy as cp
from pathlib import Path, PureWindowsPath
import numpy as np


# Import ssnake simulation
ssnake_sim = cp.load('ssnake_sim_Cq7MHz7_sigma3.csdf')
ssnake_sim /= ssnake_sim.max()
ssnake_sim.dimensions[0].to("ppm", "nmr_frequency_ratio")

# Define ranges for Cq and eta
Cq_range = np.linspace(0,20,100)
eta_range = np.linspace(0,1,10)


diso = 30; Cq= 7.7; eps_val = 0.87



# Create czjzek distribution object for site 1
quad_tensor = {"Cq": Cq, "eta": 0.01}  # Cq assumed in MHz
model_quad = ExtCzjzekDistribution(quad_tensor, eps= eps_val)
Cq_dist, eta_dist, cz_amp = model_quad.pdf(pos=[Cq_range, eta_range])
print("ssNake ExtCzjzek sigma is %s MHz" % (round(model_quad.ssNake_sigma*100)/100)) #MRsim code was modified to get ssNake parameter (2*rho)

fig,czax = plt.subplots(1, num = 1)
czax.plot(Cq_range, np.sum(cz_amp, axis = 0))

rho_z = 0.32607*(1-np.exp(-2.097*eps_val**1.151))
sigma = rho_z*np.average(Cq_range, weights = np.sum(cz_amp, axis = 0))

system = single_site_system_generator(isotope="25Mg",
                                       isotropic_chemical_shift = diso,
                                       quadrupolar={
                                           "Cq": Cq_dist * 1e6, 
                                           "eta": eta_dist}, 
                                       abundance= cz_amp                                       
                                       )

MAS = BlochDecayCTSpectrum(
    channels=["25Mg"],
    rotor_frequency=20000,  # in Hz
    spectral_dimensions=[
        SpectralDimension(spectral_width=500000.0, reference_offset=-2770.569999995587)  # values in Hz
    ],
    magnetic_flux_density = 14.1
)

sim = Simulator()
sim.config.number_of_sidebands = 32 # eight sidebands are sufficient for this example
sim.config.integration_volume = "octant"
sim.config.integration_density = 40 # (see https://mrsimulator.readthedocs.io/en/stable/user_guide/simulator/simulator.html?highlight=sphere#integration-volume)
sim.spin_systems = system  # add the spin systems
sim.methods = [MAS]
sim.run()

sim_spc = sim.methods[0].simulation
#%%
simproc = sp.SignalProcessor(
operations=[
sp.IFFT(),
sp.apodization.Gaussian(FWHM="1000 Hz"),
sp.FFT(),
]
)

MR_sim = simproc.apply_operations(sim_spc)
MR_sim /= MR_sim.max()


fig,ax = plt.subplots(1, subplot_kw={'projection':'csdm'}, num = 2)
ax.cla()
ax.plot(MR_sim.real, label = 'MRSimulator')
ax.plot(ssnake_sim, label = "ssNake")
plt.legend()
plt.show()