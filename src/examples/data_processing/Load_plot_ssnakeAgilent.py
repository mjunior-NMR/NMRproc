# -*- coding: utf-8 -*-
"""
Created on Mon May 27 17:15:25 2024

@author: Marcos
"""

import ssnmr_proc.agilent as ag
import matplotlib.pyplot as plt
from ssnmr_proc.nmrplot import nmrstyle

data_dir = r'D:\Marcos\IFSC\Dados\Agilent\Renato\A3campinas_31P_BPO4_10kHz.fid'
spc = ag.ssNakeProc(data_dir)


csdm_spec = spc.to_csdm()


plt.figure(figsize=(5, 3.5), num = 1)
ax = plt.subplot(projection="csdm")
ax.plot(csdm_spec)
plt.tight_layout()
plt.show()
nmrstyle(ax)
