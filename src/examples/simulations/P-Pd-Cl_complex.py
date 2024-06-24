# -*- coding: utf-8 -*-
"""
Created on Tue May 28 11:30:23 2024

@author: Marcos
"""

import ssnmr_proc.agilent as ag
import matplotlib.pyplot as plt
from ssnmr_proc.nmrplot import nmrstyle
from mrsimulator import Simulator, SpinSystem, Site, Coupling, SpectralDimension
from mrsimulator.spin_system.isotope import Isotope
from mrsimulator.method.lib import BlochDecaySpectrum
from mrsimulator.spin_system.tensors import SymmetricTensor
from mrsimulator import signal_processor as sp
import csdmpy as cp
import numpy as np

#%% Read experimental data and create csdm object
data_dir = r'D:\Marcos\IFSC\Dados\Agilent\Renato\A3campinas_31P_BPO4_10kHz.fid'
spc = ag.ssNakeProc(data_dir)
spc.normalize()

csdm_spec = spc.to_csdm()
csdm_spec.dimensions[0].to("ppm", "nmr_frequency_ratio") 

#%% Define system:

Cq = 100e6 # MHz
Cq_Cl = 1e6
Pd = Isotope(symbol='105Pd')
Cl35 = Isotope(symbol='35Cl')
Cl37 = Isotope(symbol='37Cl')
gamma_ratio = Cl37.gyromagnetic_ratio/Cl35.gyromagnetic_ratio
J_Pd = 200
J_Cl = 10
lb = 50 #Hz



P_site = Site(isotope="31P", isotropic_chemical_shift=147.6) #0
P_site_ = Site(isotope="31P", isotropic_chemical_shift=149.2) #1
# P_site = Site(isotope="31P", isotropic_chemical_shift=148.4) #0
# P_site_ = Site(isotope="31P", isotropic_chemical_shift=148.4) #1
Pd_site = Site(isotope="105Pd", quadrupolar=SymmetricTensor(Cq=Cq, eta=0)) #2
Cl35_site = Site(isotope="35Cl") #3
Cl37_site = Site(isotope="37Cl") #4

all_sites = [P_site, Pd_site, P_site_, Cl35_site, Cl37_site]

coupling_PPd = Coupling(site_index=[0, 1],isotropic_j=J_Pd)
coupling_PdP = Coupling(site_index=[1, 2],isotropic_j=J_Pd)
coupling_PCl35 = Coupling(site_index=[0, 3],isotropic_j=J_Cl)
coupling_PCl37 = Coupling(site_index=[0, 4],isotropic_j=J_Cl*gamma_ratio)
coupling_Cl35P = Coupling(site_index=[2, 3],isotropic_j=J_Cl)
coupling_Cl37P = Coupling(site_index=[2, 4],isotropic_j=J_Cl*gamma_ratio)

couplings_1 = [coupling_PPd, coupling_PdP, coupling_PCl35, coupling_Cl35P]
couplings_2 = [coupling_PPd, coupling_PdP, coupling_PCl37, coupling_Cl37P]
couplings_3 = [coupling_PCl35, coupling_Cl35P]
couplings_4 = [coupling_PCl37, coupling_Cl37P]


# coupling2 = Coupling(site_index=[0, 1],dipole)


# Make Site and SpinSystem objects


spin_system_1 = SpinSystem(sites=all_sites,
                           couplings = couplings_1,
                           name = '31P-105Pd-31P', abundance = Pd.natural_abundance*Cl35.natural_abundance/100
                           )
spin_system_2 = SpinSystem(sites=all_sites,
                           couplings = couplings_2,
                           name = '31P-105Pd-31P', abundance = Pd.natural_abundance*Cl37.natural_abundance/100
                           )
spin_system_3 = SpinSystem(sites=all_sites,
                           couplings = couplings_3,
                           name = '31P-105Pd-31P', abundance = (100-Pd.natural_abundance)*Cl35.natural_abundance/100
                           )
spin_system_4 = SpinSystem(sites=all_sites,
                           couplings = couplings_4,
                           name = '31P-105Pd-31P', abundance = (100-Pd.natural_abundance)*Cl37.natural_abundance/100
                           )




# Make static and MAS one-pulse acquire Method objects
mas = BlochDecaySpectrum(channels=["31P"],
                         rotor_frequency = 10000.0,
                         spectral_dimensions=[
                             SpectralDimension(count = csdm_spec.size/2, spectral_width=spc.json['sw'][0], reference_offset=spc.offset)  # values in Hz
                         ],
                         magnetic_flux_density = 5.64)  # in Hz

# Setup and run the Simulation object

systems=[spin_system_1, spin_system_2, spin_system_3, spin_system_4]



sim = Simulator(spin_systems=systems, methods=[mas])
sim.config.number_of_sidebands = 4 # eight sidebands are sufficient for this example
sim.config.integration_volume = "octant"
sim.config.integration_density = 15 
sim.config.decompose_spectrum = 'spin_system'
sim.run()


a = sim.methods[0].simulation
b = a[0]

# Create the SignalProcessor object
processor = sp.SignalProcessor(
    operations=[
        sp.IFFT(),
        sp.apodization.Exponential(FWHM=f"{lb} Hz"),                
        sp.FFT(),
    ]
)



# Apply the processor to the simulation dataset
proc_sim = processor.apply_operations(dataset=sim.methods[0].simulation)

soma = np.zeros((1,int(csdm_spec.size/2)), dtype = 'complex128')
for i in range(len(systems)):
    soma += proc_sim.y[i].components

new_dependent = cp.as_dependent_variable(soma, unit="", name = 'sum')

proc_sim.dependent_variables.append(new_dependent)


proc_sim /= proc_sim.max()[-1].real
proc_sim.dimensions[0].to("ppm", "nmr_frequency_ratio") 

plt.clf()
plt.figure(figsize=(5, 3.5), num = 1)
ax = plt.subplot(projection="csdm")
ax.plot(csdm_spec.real)
ax.plot(proc_sim.real)
plt.tight_layout()
plt.show()
plt.xlim(136,160)
nmrstyle(ax)

