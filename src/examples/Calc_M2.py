# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 16:57:52 2026

@author: mcabe
"""

import ssnmr_proc.tools as tools
import numpy as np


# filename = r'''G:\Outros computadores\PC-IFSC\IFSC\Dados\Crystal_structures\dia-Cd(Im)2_CdH.csv'
filename = r'D:\Marcos\IFSC\Dados\Crystal_structures\dia-Cd(Im)2_CdH.csv'

crys = tools.read_distances(filename)

r = crys['Cd']["distances"]
n = crys['Cd']["multiplicities"]

M2, D_eff_Hz = tools.M2("113Cd", "1H", r, n)

print(f"M2(Cd-H)           = {M2:.6e} rad²/s²")
print(f"D_eff        = {D_eff_Hz:.3f} Hz")


# filename2 = r'G:\Outros computadores\PC-IFSC\IFSC\Dados\Crystal_structures\dia-Cd(Im)2_HH.csv'
filename2 = r'D:\Marcos\IFSC\Dados\Crystal_structures\dia-Cd(Im)2_HH.csv'

crys_HH = tools.read_distances(filename2)

D_HH = []
M2_HH = []

for site, values in crys_HH.items():

    M2, D_eff_Hz = tools.M2(
        "1H",
        "1H",
        values["distances"],
        values["multiplicities"]
    )
    D_HH.append(D_eff_Hz)
    M2_HH.append(M2)

    print(f"\n{site}")
    print(f"M2(H-H)      = {M2:.6e} rad²/s²")
    print(f"D_eff   = {D_eff_Hz:.3f} Hz")
    
print("\n----------------------------------------")
print(f"Average M2(H-H) = {np.mean(M2_HH):.3f} rad²/s²")
print(f"Average D_eff(H-H) = {np.mean(D_HH):.3f} Hz")
