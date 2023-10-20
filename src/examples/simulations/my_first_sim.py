# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 16:47:00 2023

@author: Marcos
"""
from mrsimulator import Simulator, SpinSystem, Site, Coupling
from mrsimulator.spin_system.isotope import Isotope
from mrsimulator.method.lib import BlochDecaySpectrum
from mrsimulator.spin_system.tensors import SymmetricTensor
from mrsimulator import signal_processor as sp
import matplotlib.pyplot as plt


Cq = 30e6
carbon = Isotope(symbol='13C')
Cu63 = Isotope(symbol='63Cu')
Cu65 = Isotope(symbol='65Cu')
J = 1500


coupling = Coupling(site_index=[0, 1],isotropic_j=J)
coupling2 = Coupling(site_index=[0, 1],isotropic_j=J*Cu65.gyromagnetic_ratio/Cu63.gyromagnetic_ratio)


# Make Site and SpinSystem objects
P_site = Site(isotope="31P", isotropic_chemical_shift=0)
Cu63_site = Site(isotope="63Cu", quadrupolar=SymmetricTensor(Cq=Cq, eta=0))
Cu65_site = Site(isotope="65Cu", quadrupolar=SymmetricTensor(Cq=Cq*Cu65.gyromagnetic_ratio/Cu63.gyromagnetic_ratio , eta=0))

spin_system_63 = SpinSystem(sites=[P_site, Cu63_site], abundance = Cu63.natural_abundance, couplings = [coupling])
spin_system_65 = SpinSystem(sites=[P_site, Cu65_site], abundance = Cu65.natural_abundance, couplings = [coupling2])


# Make static and MAS one-pulse acquire Method objects
static = BlochDecaySpectrum(channels=["31P"])
mas = BlochDecaySpectrum(channels=["31P"], rotor_frequency = 20000.0)  # in Hz

# Setup and run the Simulation object
sim = Simulator(spin_systems=[spin_system_63,spin_system_65], methods=[mas])
sim.run()

# Create the SignalProcessor object
processor = sp.SignalProcessor(
    operations=[
        sp.IFFT(),
        sp.apodization.Exponential(FWHM="150 Hz"),        
        sp.FFT(),
    ]
)

# Apply the processor to the simulation dataset
processed_simulation = processor.apply_operations(dataset=sim.methods[0].simulation)

# Plot the spectra
# plt.figure(figsize=(6, 3))
ax = plt.subplot(projection="csdm")
ax.plot(processed_simulation.real)
ax.set_title("Static")
plt.xlim(30,-30)
plt.show()


processed_simulation.real.save('simulation.csdf')

# a = sim.methods[0].simulation.real