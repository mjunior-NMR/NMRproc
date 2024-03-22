# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 15:50:03 2024

@author: Marcos
"""
def nmrstyle(ax, nucleus = None, label = r'Chemical shift /ppm'):
    import matplotlib.pyplot as plt
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_yticks([])
    ax.set_ylabel('')
    ax.invert_xaxis()
    if nucleus != None:
        import re
        split_str = re.split(r'(\d+)', nucleus) # Slipts string into numbers and text 
        isotope = split_str[1]
        atom = split_str[2]
        ax.set_xlabel(r'$^{'+ isotope + r'}$' + atom + r'-$\delta$ /ppm', fontsize = 16)
    else:
        ax.set_xlabel(label, fontsize = 16)
    plt.xticks(fontsize=14)
    # print('Plot has been converted to NMR style!')