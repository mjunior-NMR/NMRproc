# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 17:22:47 2024

@author: mcabe
"""

import matplotlib.pyplot as plt
import ssnmr_proc.bruker as br

import numpy as np


# data = np.random.rand(1024,512)
# data_dir = r'G:\Outros computadores\PC-IFSC\IFSC\Dados\Bruker\Renato\13C_1H_MOF808_2.5mm\8\pdata\1'
# data_dir = r'D:\Marcos\IFSC\Dados\Bruker\Renato\13C_1H_MOF808_2.5mm\8\pdata\1'
data_dir = r'D:\Marcos\IFSC\dados\Bruker\Renato\13C_1H_MOF808+EDTA_posFC_2.5mm\8\pdata\1'

spc = br.ProcData(data_dir)
spc.normalize()

spc.plot2D(f1lim=[200,0],f2lim=[12,-5],proj_type=['sum','skyline'])

# Salva o gráfico
plt.savefig('grafico_comparacao.png', dpi=600)
plt.show()