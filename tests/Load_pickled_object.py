# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 17:12:03 2023

@author: Marcos
"""

import pickle
# import ssnmr_proc.bruker
import matplotlib.pyplot as plt

with open('procdata.pkl' , 'rb') as f:
    pd = pickle.load(f)

fig, ax = plt.subplots(sharex=True,sharey=True,figsize=(6, 5))
pd.plot(axis = ax)

