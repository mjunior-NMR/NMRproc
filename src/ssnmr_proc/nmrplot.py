# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 15:50:03 2024

@author: Marcos
"""
def nmrstyle(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_yticks([])
    ax.set_ylabel('')
    ax.invert_xaxis()
    print('Plot has been converted to NMR style!')