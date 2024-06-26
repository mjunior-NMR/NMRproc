# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 15:27:29 2023

@author: Marcos de Oliveira Jr.
"""
import numpy as np
import json as js
import csdmpy as cp


def temp():
    return

def to_csdm(file_path): 
    """ A function to convert ssNake json files to csdm format. """
    f = open(file_path)
    json = js.load(f)
    data = np.array(json['dataReal'][0])+1j*np.array(json['dataImag'][0])
    npts = data.size
    hz_scale = json['xaxArray'][0]
    offset = json["freq"][0]-json["ref"][0]

    try:
        dv = cp.as_dependent_variable(data, unit="")
        dim = cp.LinearDimension(
        count = npts, 
        origin_offset = f'{json["freq"][0]} Hz',  
        #coordinates_offset = f'{spec["metaData"]["Offset [Hz]"]} Hz',
        coordinates_offset = f'{offset} Hz',
        increment = f'{hz_scale[1]-hz_scale[0]} Hz',
        complex_fft=True,
        label="Frequency",
        reciprocal={'quantity_name': 'time'}
        )
        csdm_spec = cp.CSDM(dependent_variables=[dv], dimensions=[dim])
        csdm_spec.dimensions[0].to("ppm", "nmr_frequency_ratio")         
    except:
        raise ImportError("csdmpy must be installed to use this function. Please install by typing 'pip install csdmpy' in the terminal.")    
    return csdm_spec