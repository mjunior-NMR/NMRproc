# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 16:19:35 2021

@author: Henrik Bradtmüller - mail@bradtmueller.net - https://hbrmn.github.io/

todo: make outputstring stating the nucleus-nucleus interactions
-make it auto detect delimiters
-make it detect coordination from shortest X-O distances
-include distance cut-off to be able to separately look into contributions from
the first few and higher coordination spheres

- Code the CIF to dist-list part
"""

import os
import numpy as np
from scipy import constants as cnst
# From ssnake
from loadIsotopes import getIsotopeInfo as gii

class DataSet():

    def __init__(self, data_path, nuc1, nuc2):

        isoPath = (os.path.dirname(os.path.realpath(__file__)) + os.path.sep
           + "IsotopeProperties")
        self.isotopes = gii(isoPath)
        self.data = np.genfromtxt(data_path,
                     dtype=str, delimiter=' ', encoding=None)
        self.data = np.split(self.data, self.find_distinct_indices(), axis=0)
        self.ind1 = self.find_isotope_index(nuc1)
        self.ind2 = self.find_isotope_index(nuc2)
        self.gamma1 = self.isotopes['gamma'][self.ind1] * 1e7
        self.gamma2 = self.isotopes['gamma'][self.ind2] * 1e7
        self.quant_number = self.isotopes['spin'][self.ind2]
        self.abundance = self.isotopes['abundance'][self.ind2]/100
        print(self.gamma1)
        print(cnst.hbar)

    def find_isotope_index(self, nuc):
        try:
            mass = float(nuc.rstrip('AaBbCcDdEeFfGgHhIiJjKkLlMmNn' +
                                       'OoPpQqRrSsTtUuVvWwXxYyZz'))
        except:
            raise ValueError('Nuclei have to be input in the format' +
                             '"mass number" + "atomic symbol",' +
                             ' e.g., 11B or 125Te.') from None
        ind = np.where(mass == np.array(self.isotopes['atomMass']))[0][0]
        return ind

    def find_distinct_indices(self):
        """ Splits data into arrays according to individual crystal sites"""
        strings = np.unique(self.data[:,0])
        indices = [np.where(self.data[:,0] == strings[x])[0][0]
                   for x in range(1, len(strings))]
        return indices

    def dist_sum(self, x):
        """Calculates the sum of inverse sixth power distances, multiplied
        by their corresponding multiplicities"""
        inv_sp_dist = ((np.array(self.data[x][:,3].astype(np.float))
                          * 1e-10)**6 )**-1
        multiplicity = np.array(self.data[x][:,2].astype(np.float))
        print(np.sum(inv_sp_dist*multiplicity))
        return np.sum(inv_sp_dist*multiplicity)

    def m2het(self, dist_sum):
        """Calculates heteronuclear second moment from given distance sum"""
        vac_perm = cnst.physical_constants['vacuum mag. permeability'][0]
        return (self.abundance * (4/15) * (vac_perm / (4 * cnst.pi))**2 *
                self.quant_number * (self.quant_number + 1)*
                self.gamma1**2 * self.gamma2**2 * cnst.hbar**2 * dist_sum)

    def m2hom(self, dist_sum):
        """Calculates homonuclear second moment from given distance sum"""
        vac_perm = cnst.physical_constants['vacuum mag. permeability'][0]
        return (self.abundance * 0.9562 * (vac_perm / (4 * cnst.pi))**2 *
                self.gamma1**2 * self.gamma2**2 * cnst.hbar**2 * dist_sum)

    def second_moment(self):
        """Calculates second moments for all "observe nuclei x" in dataset """
        if self.gamma1 == self.gamma2:
            result = [self.m2hom(self.dist_sum(x))/1e6
                      for x in range(len(self.data))]
        else:
            result = [self.m2het(self.dist_sum(x))/1e6
                      for x in range(len(self.data))]
        return result

##########################################################
a = DataSet(r'C:\Users\edwu5ea1\data_work\XRD Data\ScF3\ScF3.csv', '45Sc', '19F')
m2 = a.second_moment()
const_d = np.mean(m2)/(2*np.pi)