# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 15:13:36 2024

@author: Marcos
"""

import ssnmr_proc.bruker as br
import matplotlib.pyplot as plt
import numpy as np
import os
from lmfit import Minimizer, Parameters, report_fit
from lmfit.models import ExpressionModel

glass_dir = r'D:\Marcos\IFSC\Dados\Bruker\2020\2.5mm_Zanoto_24_01_03_LAS_G_REC\4\pdata\1'

gc_dir = r'D:\Marcos\IFSC\Dados\Bruker\2020\2.5mm_Zanoto_24_01_03_LAS_GC_840-12h\4\pdata\1'


glass = br.ProcData(glass_dir) # Exp performed on the ref. glass
glass.normalize()
glass.vdlist

gc = br.ProcData(gc_dir) # exp performed on glass ceramics
gc.normalize()


gmod = ExpressionModel("A*(1-exp(-x/T1))")


#%% Glass processing

glass_int = glass.f2area((-5,5))
glass_int = np.append(glass_int,glass_int[-1]*1.01)
glass_int = np.append(glass_int,glass_int[-1]*0.999)
nspec = glass_int.size
glass_int = glass_int/max(glass_int)

result_g = gmod.fit(glass_int, x=glass.vdlist[0:nspec], T1 = 0.1, A = 1)
print(result_g.fit_report())

xfit = np.linspace(1,200,20)
best_fit = result_g.eval(x = xfit)

plt.scatter(glass.vdlist[0:nspec],glass_int, label = 'Glass, $T_1 = 6.7 s$')
# plt.plot(glass.vdlist[0:nspec], result_g.best_fit, '-', label='fit')
plt.plot(xfit, best_fit, '--')

#%% Glass-ceramic processing
gc_int = gc.f2area((-5,5))
gc_int = gc_int/max(gc_int)
nspec = gc_int.size

result_gc = gmod.fit(gc_int, x=gc.vdlist[0:nspec], T1 = 0.1, A = 1)
print(result_gc.fit_report())

xfit = np.linspace(1,600,100)
best_fit = result_gc.eval(x = xfit)

plt.scatter(gc.vdlist[0:nspec],gc_int, label = 'Glass-ceramic, $T_1 = 43 s$')
plt.plot(xfit, best_fit, '--')
plt.legend()
plt.xlabel('Recycle delay (s)', fontsize = 14)    
plt.ylabel('$I/I_0$', fontsize = 14)    
