import numpy as np
import openmc
import openmc.capi as capi

class Polynomials(object):
    """Class for Polynomials"""
    def __init__(self, coeff):
        self.coeff = coeff
        self.n_coeff = len(coeff)

class Legendre(Polynomials):
    def __init__(self, coeff, domain_min = -1, domain_max = 1):
        super().__init__(coeff)
        self.order = self.n_coeff - 1
        self.domain_min = domain_min
        self.domain_max = domain_max
        self.norm_vec = (2 * np.arange(self.n_coeff) + 1) / (self.domain_max
            - self.domain_min)
        self.norm_coeff = np.multiply(self.norm_vec,self.coeff)

    def __call__(self, x):
        # x is the normalized position on [-1,1]
        pn = capi.calc_pn(self.order, x)
        return np.sum(np.multiply(self.norm_coeff, pn))

class ZernikeRadial(Polynomials):
    def __init__(self,coeff,radius = 1):
        super().__init__(coeff)
        self.order = 2 * (self.n_coeff - 1)
        self.radius = radius
        self.norm_vec = (2 * np.arange(self.n_coeff) + 1) / (np.pi * 
            radius ** 2)
        self.norm_coeff = np.multiply(self.norm_vec,self.coeff)
        
    def __call__(self, r):
        # r is the normalized position on radius [0,1]
        zn_rad = capi.calc_zn_rad(self.order, r)
        return np.sum(np.multiply(self.norm_coeff, zn_rad))