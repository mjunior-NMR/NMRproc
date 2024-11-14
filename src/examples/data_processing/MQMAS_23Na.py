# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 17:22:47 2024

@author: mcabe
"""

import matplotlib.pyplot as plt

from mrsimulator import Simulator, SpinSystem, Site
from mrsimulator.method.lib import ThreeQ_VAS
from mrsimulator import signal_processor as sp
from mrsimulator.spin_system.tensors import SymmetricTensor
from mrsimulator.method import SpectralDimension
import ssnmr_proc.bruker as br

exp = br.ProcData(r'G:\Outros computadores\PC-IFSC\IFSC\Dados\Bruker\Yara\Glass Digital 1 - PS\10\pdata\1')

exp_csdm = exp.to_csdm()

Na_1 = Site(
    isotope="23Na",
    isotropic_chemical_shift=-27.4,  # in ppm
    quadrupolar=SymmetricTensor(Cq=1.68e6, eta=0.2),  # Cq is in Hz
)

sites = [Na_1]  # all sites

spin_systems = [SpinSystem(sites=[s]) for s in sites]


method = ThreeQ_VAS(
    channels=["23Na"],
    magnetic_flux_density=9.4,  # in T
    spectral_dimensions=[
        SpectralDimension(
            count=128,
            spectral_width=7e3,  # in Hz
            reference_offset=-7e3,  # in Hz
            label="Isotropic dimension",
        ),
        SpectralDimension(
            count=256,
            spectral_width=1e4,  # in Hz
            reference_offset=-4e3,  # in Hz
            label="MAS dimension",
        ),
    ],
)

sim = Simulator()
sim.spin_systems = spin_systems  # add the spin systems
sim.methods = [method]  # add the method.
sim.run()

dataset = sim.methods[0].simulation

plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
cb = ax.imshow(dataset.real / dataset.real.max(), aspect="auto", cmap="gist_ncar_r")
cb = ax.imshow(exp_csdm.real / exp_csdm.real.max(), aspect="auto", cmap="gist_ncar_r")
plt.colorbar(cb)
ax.invert_xaxis()
ax.invert_yaxis()
plt.tight_layout()
plt.show()