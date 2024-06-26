# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 08:54:31 2024

@author: Marcos
"""

import ssnmr_proc.ssnake as ss
import ssnmr_proc.bruker as br
import matplotlib.pyplot as plt
from ssnmr_proc.nmrplot import nmrstyle
from mrsimulator import Simulator, SpinSystem, Site, Coupling, SpectralDimension
from mrsimulator.spin_system.isotope import Isotope
from mrsimulator.method.lib import BlochDecayCTSpectrum
from mrsimulator.spin_system.tensors import SymmetricTensor
from mrsimulator import signal_processor as sp
import csdmpy as cp
import numpy as np

#%% Read experimental data and create csdm object
data_dir = r'D:\Marcos\IFSC\Dados\Bruker\references\47-49Ti_TiO2_3.2mm\94\pdata\1'
exp = br.ProcData(data_dir)
spc = exp.to_csdm()
# data_dir = r'D:\Marcos\IFSC\Dados\Bruker\references\47-49Ti_TiO2_3.2mm\95\pdata\1\TiO2_HE.json'
# spc = ss.to_csdm(data_dir)
spc /= spc.max().real

glass_data_dir = r'D:\Marcos\IFSC\Dados\Bruker\2020\3.2mm_Philip_24_06_21_TiO2P2O5\91\pdata\1\RLS_glass_Ti.json'
glass = ss.to_csdm(glass_data_dir)


#%% Define system:

Cq49 = 4.8e6 # Hz
diso49 = -913 #ppm
daniso = -60 #ppm 
lb = 400 #Hz
lb2 = 300 #Hz
eta = 0.0
Ti49 = Isotope(symbol='49Ti')
Ti47 = Isotope(symbol='47Ti')

Cq47 = Cq49*Ti47.quadrupole_moment/Ti49.quadrupole_moment
diso47 = diso49-268

Ti49_site = Site(isotope="49Ti",
                 isotropic_chemical_shift = diso49,  # in ppm
                 shielding_symmetric = SymmetricTensor(zeta=daniso, eta=0),
                quadrupolar = SymmetricTensor(Cq=Cq49, eta= eta))

Ti47_site = Site(isotope="47Ti",
                 isotropic_chemical_shift = diso47,  # in ppm
                 shielding_symmetric = SymmetricTensor(zeta=daniso, eta=0),
                quadrupolar = SymmetricTensor(Cq=Cq47, eta= eta))

system_49 = SpinSystem(sites=[Ti49_site], name = '$^{49}$Ti', abundance = Ti49.natural_abundance)
system_47 = SpinSystem(sites=[Ti47_site], name = '$^{47}$Ti', abundance = Ti47.natural_abundance)

MAS49 = BlochDecayCTSpectrum(
    channels=["49Ti"],
    rotor_frequency= 22000,  # in Hz
    magnetic_flux_density= 14.1,  # in tesla
    spectral_dimensions=[
        SpectralDimension(
            count=spc.dimensions[0].count,
            spectral_width=spc.dimensions[0].increment.value*spc.dimensions[0].count,  # in Hz
            reference_offset=spc.dimensions[0].coordinates_offset.value,  # in Hz
        )
    ],
)

MAS47 = BlochDecayCTSpectrum(
    channels=["47Ti"],
    rotor_frequency= 22000,  # in Hz
    magnetic_flux_density= 14.1,  # in tesla
    spectral_dimensions=[
        SpectralDimension(
            count=spc.dimensions[0].count,
            spectral_width=spc.dimensions[0].increment.value*spc.dimensions[0].count,  # in Hz
            reference_offset=spc.dimensions[0].coordinates_offset.value,  # in Hz
        )
    ],
)

systems=[system_49, system_47]


sim = Simulator(spin_systems=systems, methods=[MAS49,MAS47])
sim.config.number_of_sidebands = 4 # eight sidebands are sufficient for this example
sim.config.integration_volume = "hemisphere"
sim.config.integration_density = 90 
sim.run()

processor = sp.SignalProcessor(
    operations=[
        sp.IFFT(),
        sp.apodization.Exponential(FWHM=f"{lb} Hz"),                
        sp.FFT(),
    ]
)


# Apply the processor to the simulation dataset
sim_spc = processor.apply_operations(dataset=sim.methods[0].simulation) 

processor = sp.SignalProcessor(
    operations=[
        sp.IFFT(),
        sp.apodization.Exponential(FWHM=f"{lb2} Hz"),                
        sp.FFT(),
    ]
)

sim2_spc = processor.apply_operations(dataset=sim.methods[1].simulation)

sim2_spc = sim2_spc/sim_spc.max().real
sim_spc /= sim_spc.max().real

plt.clf()
fig,ax = plt.subplots(3,1, sharex=True, figsize=(5, 3.5), num = 1, subplot_kw={"projection": "csdm"}, squeeze=True)
ax[1].plot(spc.real)
ax[1].text(1800, 0.2, r'TiO$_2$', fontsize=15)
ax[2].plot(sim_spc.real, label = '$^{49}$Ti')
ax[2].plot(sim2_spc.real, label = '$^{47}$Ti')
ax[2].text(1800, 0.2, r'Simulation (TiO$_2$)', fontsize=15)
ax[0].plot(glass)
ax[0].text(1800, 0.5, r'Glass (x = 0.71)', fontsize=15)
plt.legend(fontsize = 16)
plt.xlim(-4000,2000)
fig.subplots_adjust(hspace=0)
# plt.tight_layout()
# plt.show()
nmrstyle(ax[0])
nmrstyle(ax[1])
nmrstyle(ax[2])
ax[0].set_yticks([])
ax[0].set_ylabel('')

