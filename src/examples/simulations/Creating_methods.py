# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 10:40:02 2023

@author: Marcos
"""

hahn_echo = Method(
channels=["2H"],
magnetic_flux_density=9.4, # in T
spectral_dimensions=[
{"count": 512,
"spectral_width": 2e4, # in Hz
"events": [
SpectralEvent(fraction=0.5, transition_query=[{"ch1": {"P": [1]}}]),
MixingEvent(mixing_query={"ch1": {"tip_angle": np.pi, "phase": 0}}),
SpectralEvent(fraction=0.5, transition_query=[{"ch1": {"P": [-1]}}])]}],
)
solid_echo = Method(
channels=["2H"],
magnetic_flux_density=9.4, # in T
spectral_dimensions=[
{"count": 512,
"spectral_width": 2e4, # in Hz
"events": [
SpectralEvent(fraction=0.5, transition_query=[{"ch1": {"P": [-1]}}]),
MixingEvent(mixing_query={"ch1": {"tip_angle": np.pi/2, "phase": 0}}),
SpectralEvent(fraction=0.5, transition_query=[{"ch1": {"P": [-1]}}])]}],
)
