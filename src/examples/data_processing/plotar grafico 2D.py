import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import ssnmr_proc.bruker as br

data_dir = r'/content/1 - Zn ZIF-4 seco - sqdq2D_R1/pdata/1'

spc = br.ProcData(data_dir)

Z = spc.data/spc.data.real.max()
y = spc.ppm_scale
x = spc.ppm_scale_1

x.ndim


fig, axs = plt.subplots(2, 2, gridspec_kw={'height_ratios': [1, 5], 'width_ratios': [5, 1]})

# plot parameters
contour_start = 0.004           # contour level start value
contour_num = 10                # number of contour levels
contour_max = 1
contour_factor = np.exp(np.log(contour_max/contour_start)/contour_num)         # scaling factor between contour levels
# contour_factor = 2.1          # scaling factor between contour levels

# calculate contour levels
cl = contour_start * contour_factor ** np.arange(contour_num)
# cl = np.linspace(0.004,0.9,15)

proj = Z.sum(axis=0)
proj = proj/max(proj)*2.5+5

axs[0, 0].plot(x,proj, linestyle = 'solid', linewidth = 1, color = 'k')
axs[0, 0].set_xlim([12,-1])
axs[0, 0].spines['top'].set_visible(False)
axs[0, 0].spines['right'].set_visible(False)
axs[0, 0].spines['bottom'].set_visible(False)
axs[0, 0].spines['left'].set_visible(False)
axs[0, 0].xaxis.set_visible(False)
axs[0, 0].yaxis.set_visible(False)
axs[0, 0].set_ylim([5,8])

plt.cla()
axs[1, 0].contour(x,y,Z, cl, cmap='winter')
axs[1, 0].set_ylim([18,-1])
axs[1, 0].set_xlim([12,-1])
axs[1, 0].set_xlabel(r'$^1$H-SQ /ppm', fontsize = 16)
axs[1, 0].set_ylabel(r'$^1$H-DQ /ppm', fontsize = 16)
axs[1, 0].plot(x,2*x, linestyle = 'dashed', linewidth = 1.5, color = 'r')
plt.tight_layout()

axs[0, 1].spines['top'].set_visible(False)
axs[0, 1].spines['right'].set_visible(False)
axs[0, 1].spines['bottom'].set_visible(False)
axs[0, 1].spines['left'].set_visible(False)
axs[0, 1].xaxis.set_visible(False)
axs[0, 1].yaxis.set_visible(False)


axs[1, 1].plot(proj, 2 * x, linestyle = 'solid', linewidth = 1, color = 'k')
axs[1, 1].set_ylim([18,-1])
# axs[1, 1].set_xlim([12,-1])
axs[1, 1].spines['top'].set_visible(False)
axs[1, 1].spines['right'].set_visible(False)
axs[1, 1].spines['bottom'].set_visible(False)
axs[1, 1].spines['left'].set_visible(False)
axs[1, 1].xaxis.set_visible(False)
axs[1, 1].yaxis.set_visible(False)
axs[1, 1].set_xlim([5,8])


plt.subplots_adjust(hspace=0.01, wspace=0.01)
plt.savefig('grafico.png', dpi=600)