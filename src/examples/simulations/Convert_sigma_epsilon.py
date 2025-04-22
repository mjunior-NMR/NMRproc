# -*- coding: utf-8 -*-
"""
Created on Wed Mar 26 09:51:36 2025

Convert epsilon MRSimulator values to Sigma ssNake parameter for extended Czjzek.
The ssNake sigma is twice the rho parameter from Caër et al. J.Phys.:Condens.Matter22(2010)065402

@author: Marcos
"""
import numpy as np
eps = 0.83
Cq = 5.3

rho = 2*eps*np.sqrt(1.5)*Cq/np.sqrt(30)

print(eps,rho)