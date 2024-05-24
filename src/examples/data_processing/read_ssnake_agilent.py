# -*- coding: utf-8 -*-
"""
Created on Fri May 24 16:03:50 2024

@author: Marcos
"""

import ssnmr_proc.agilent as ag

data_dir = r'D:\Marcos\IFSC\Dados\Agilent\Renato\A3campinas_31P_BPO4_10kHz.fid'

spec = ag.ssNakeProc(data_dir)
