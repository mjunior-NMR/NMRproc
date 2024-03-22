# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 15:27:29 2023

@author: Marcos de Oliveira Jr.
"""
import numpy as np
import nmrglue as ng
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator)

class ProcData(): 
    """ A Class to work with Bruker pre-processed data:
        spc = ProcData('data_dir') will return:
            spc.dic --> nmrglue dictionary with all acqu ans procs metadata
            spc.data.real --> Real spectrum
            spc.data.imag --> Imag spectrum
            spc.udic = nmrglue universal dictionary
            spc.ppm_scale = ppm axis
            spc.normalize(method = 'area' or 'intensity', range = 'full' or a tuple with integration region) Normalize spectrum by maximum or area
            spc.plot() --> plot spectrum with NMR-like style"""
    def __init__(self,data_dir):
        #Data dir must be a directory with processed Bruker 1r or 2rr files.
        self.dic, tmp = ng.bruker.read_pdata(data_dir, scale_data=True,all_components = True)
        self.data = tmp[0]+1j*tmp[1]        
        self.udic = ng.bruker.guess_udic(self.dic, self.data.real) # Universal nmrglue dictionary
        self.uc = ng.fileiobase.uc_from_udic(self.udic, dim = 1)        
        self.ppm_scale = self.uc.ppm_scale() # ppm axis        
        self.hz_scale = self.uc.hz_scale() # hz axis
        
        if self.data.ndim == 2:
            self.uc1 = ng.fileiobase.uc_from_udic(self.udic, dim=1)
            self.ppm_scale_1 = self.uc1.ppm_scale() # ppm axis        
            self.hz_scale_1 = self.uc1.hz_scale() # hz axis    
            self.data = self.data.transpose()
        
        self.fid = ng.process.proc_base.ifft(tmp[0]+1j*tmp[1])
        self.udic[0]['time'] = False
        self.udic[0]['freq'] = True
        self.reffrq = self.dic['acqus']['SFO1']

#%% Get area of an spectral region defined from a tupl variable
    def area(self,region = tuple()):
        if len(region) == 0:
            area = abs(np.trapz(self.data.real,x = self.ppm_scale))            
        else:
            p1 = self.uc(str(region[0])+' ppm') #convert from ppm_scale to points p1 and p2
            p2 = self.uc(str(region[1])+' ppm')
            reduced_rdata = self.data.real[min(p1,p2):max(p1,p2)]
            reduced_ppm_scale = self.ppm_scale[min(p1,p2):max(p1,p2)]
            area = abs(np.trapz(reduced_rdata,reduced_ppm_scale))
        return area


#%%  Normalize data by max intensity or area within region
    def normalize(self,method = 'intensity', region = tuple()):        
        if method == 'intensity':
            norm = 1/self.data.real.max()            
            # self.data.imag = self.data.imag*norm  # real and imaginary are normalized by intensity of real data
            self.data = self.data*norm
            # self.data.real = self.data.real*norm
            
        elif method == 'area':
            if region == ():
                raise ValueError(r'Region must be informed in Tuple format')
            norm = 1/self.area(region)
            self.data = self.data*norm
        elif method == '01': #Normaliza entre zero e um
            minimo = min(self.data.real)
            self.data = self.data-minimo
            
            norm = 1/max(self.data.real)
            self.data = self.data*norm                    
        else:
            raise ValueError(r'method must be either "intensity" or "area"!')
                
#%%    '''Plot single spectrum'''
    def plot(self,
             axis, # axis is an axes object, e.g. fig, ax = plt.subplots()
             nucleus=None,
             major_ticks_space = None,
             minor_ticks_space = None,
             line_color = 'k',
             font_name = 'Times new roman'):
        
        axis.plot(self.ppm_scale, self.data.real, color = line_color)
        axis.set_yticks([])
        axis.spines['top'].set_visible(False)
        axis.spines['right'].set_visible(False)
        axis.spines['left'].set_visible(False)
       
        if nucleus != None:
            import re
            split_str = re.split(r'(\d+)', nucleus) # Slipts string into numbers and text 
            isotope = split_str[1]
            atom = split_str[2]
            axis.set_xlabel(r'$^{'+ isotope + r'}$' + atom + r'-$\delta$ /ppm', fontsize = 16)
        else:
            axis.set_xlabel(r'Chemical shift /ppm', fontsize = 16)
            
        plt.xticks(fontsize=14, fontname = font_name)
        
        if major_ticks_space != None:
            axis.xaxis.set_major_locator(MultipleLocator(major_ticks_space))
        axis.xaxis.set_major_formatter('{x:.0f}')
        
        if major_ticks_space != None:
            axis.xaxis.set_minor_locator(MultipleLocator(minor_ticks_space))
    

#%%   '''Save proc_spec as csdm file'''
    def to_csdm(self):
        """
        Return a csdm object containing processed data.

        
        """

        try:
            import csdmpy as cp
            dv = cp.as_dependent_variable(self.data, unit="")
            dim = cp.LinearDimension(
                count = self.udic[0]['size'], 
                origin_offset = f'{self.udic[0]["obs"]*1e6} Hz',  
                coordinates_offset = f'{self.udic[0]["car"]} Hz',
                increment = f'{self.hz_scale[1]-self.hz_scale[0]} Hz',
                complex_fft=True,
                label="Frequency",
                reciprocal={'quantity_name': 'time', 'label': self.udic[0]['label']}
                )
            csdm_spec = cp.CSDM(dependent_variables=[dv], dimensions=[dim])
            csdm_spec.dimensions[0].to("ppm", "nmr_frequency_ratio")            
            return csdm_spec
        except:
            raise ImportError("csdmpy must be installed to use this function. Please install by typing 'pip install csdmpy' in the terminal.")




#%% Forward FFT (FFT)
    def FFT(self):
        self.data = ng.process.proc_base.fft(self.data)
        

#%% Inverse iFFT (iFFT)
    def iFFT(self):
        self.data = ng.process.proc_base.ifft(self.data)

#%% Shift_data
    def ls(self,pts=0):
        self.data = ng.process.proc_base.ls(self.data, pts)

#%% Flip data
    def flip(self):
        self.data = np.flip(self.data)
        
#%%    '''Save spectrum to file'''
    def save(self,filename = r'data.dat', file_type = r'ascii', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None):
        
        if file_type == 'ascii':
            self.xydata = np.column_stack([self.ppm_scale,self.data.real])
            np.savetxt(filename, 
                       self.xydata, 
                       fmt='%.18e', 
                       delimiter=' ', 
                       newline='\n', 
                       header='', 
                       footer='', 
                       comments='# ', 
                       encoding=None)
        elif file_type == 'dmfit': #Readable by dmfit software
            self.xydata = np.column_stack([self.hz_scale,self.data.real])
            np.savetxt(filename, 
                       self.xydata, 
                       fmt='%.18e', 
                       delimiter=' ', 
                       newline='\n', 
                       header='ti:\t exported from NMRproc \n##freq '+f'{self.reffrq}', 
                       footer='', 
                       comments='', 
                       encoding=None)
        elif file_type == 'csdf': # serialized csdm file format
            csdm_output = self.to_csdm()
            csdm_output.save(filename)
        elif file_type == 'pickle':
            import pickle
            with open(filename, 'wb') as file:
                pickle.dump(self, file)
                print(f'Object successfully saved to "{filename}"')
        else:
            print('Select a valid file_format: "ascii", "dmfit", "csdf", "pickle"')
            
'''======================== End of class ProcData ==================================='''


def load_pdata_obj(file_dir):
    
    import pickle
    
    with open(file_dir, 'rb') as f:
        procdata = pickle.load(f)
    return procdata




#%%
# def ppm_scale(dic,data):
    
#     # read in the data
#     data_dir = r"G:\Outros computadores\PC-IFSC\IFSC\Dados\Bruker\Pablo_ZIFs\2.5mm_Fiy3-m\3\pdata\1"

#     # From pre-procced data.
#     dic, data = ng.bruker.read_pdata(data_dir, scale_data=True)

#     udic = ng.bruker.guess_udic(dic, data)
#     uc = ng.fileiobase.uc_from_udic(udic)
#     ppm_scale = uc.ppm_scale()
    
#     return ppm_scale

# def plot_all_spectra(samples,mainpath,line_color = 'k',font_name = 'Times new roman',xlim = None,ylim=None,text_xpos = 0, text_ypos = 0,ticks_space = 20):
#     #Function to load and plot all 13C spectra creating the paths from the keys in samples dictionary.
#     spec = dict()
#     XAxis = dict()
#     for key in samples: 
#         Path = mainpath + '2.5mm_' + key + r'\3\pdata\1'
#         dic,data= ng.bruker.read_pdata(Path,all_components= 'true')
#         XAxis[key] = ppm_axis(dic) #Function created by me to reconstruct the x-axis
#         spec[key] = data[0]/max(data[0])
        
        
#     num_plots = len(samples) # Number of spectra to be ploted
    
#     #Create plot environment
#     fig, axs = plt.subplots(num_plots, 1, figsize=(5, 8),sharex=True,sharey=True,facecolor='none')
#     # Remove vertical space between axes
#     fig.subplots_adjust(hspace=-0.01)
    
    
    
#     for i, ax in enumerate(axs.flatten()):
#         # Remove eixo y e arestas dos plots
#         ax.set_yticks([])
#         ax.spines['top'].set_visible(False)
#         ax.spines['right'].set_visible(False)
#         ax.spines['left'].set_visible(False)
#         ax.text(text_xpos,text_ypos, samples[list(samples)[i]],fontname=font_name, fontsize = 12)
#         # if i != num_plots-1: # Mantem eixo x no plot inferior e configura x-label
#         #      ax.set_xticks([])
#         #     ax.spines['bottom'].set_visible(False)                
#         ax.xaxis.set_major_locator(MultipleLocator(ticks_space))
#         ax.xaxis.set_major_formatter('{x:.0f}')
#         ax.xaxis.set_minor_locator(MultipleLocator(5))
#         # user can choose plot scale
#         if xlim != None:
#             plt.xlim(xlim)
#         if ylim != None:
#             plt.ylim(ylim)
#         ax.plot(XAxis[list(spec)[i]],spec[list(spec)[i]],color = line_color)
#         ax.set_clip_on=False
    
#     plt.show()
#     ax.set_xlabel('$^{13}$C-$\delta$ /ppm', fontsize = 14)
#     return