# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 15:00:14 2023

@author: mcabe
"""

import ssnmr_processing.bruker as br
# import nmrglue as ng
import matplotlib.pyplot as plt
# from matplotlib.ticker import (MultipleLocator)
# import numpy as np

path = r'G:\Outros computadores\PC-IFSC\IFSC\Dados\Bruker\camila\2.5mm_Cdbtc_1hBM_2022-11-04\3\pdata\1'

spc = br.ProcData(path)  #Atribui a classe bruker.ProcData à variável spc
spc.normalize() #Normaliza o espectro


#Preparar eixo para plot
fig, ax = plt.subplots(sharex=True,sharey=True,figsize=(5, 5))

''' plot() é uma função da classe br.ProcData que plota o 
espectro com atribtos definidos para apresentação de espectros de RMN'''
spc.plot(axis = ax, nucleus = '13C') 
ax.set_xlim([190,110])
ax.set_ylim([-0.05,1])

spc.save(r'spectrum.dat') # Utiliza o método da classe bruker.ProcData para salvar o espectro em um arquivo.