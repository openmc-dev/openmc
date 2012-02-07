#!/usr/bin/python
# Filename: eigenfunction_rms.py

import h5py
import numpy as np

class EigenFunction:

 	'''Represents the reference eigenfunction'''
 	reference = 0.

	def __init__(self,h5_file,cycle):

		'''Initializes the eigenfunction'''
		self.cycle = cycle
		f = h5py.File(h5_file,'r')
		group = '/cycle'+str(cycle)+'/flux'
		dataset = f[group]
		self.function = np.empty(dataset.shape,dataset.dtype)
		dataset.read_direct(self.function)
		self.iamref = 'F'

	def set_reference(self):

		'''Sets instance to be reference calc'''
		self.iamref = 'T'
		EigenFunction.reference = self.function

	def compute_rms(self):

		'''Computes RMS value'''
		Np = self.function.size
                tmp = (self.function - EigenFunction.reference)*(self.function - EigenFunction.reference)
                tmp2 = tmp.sum()
                self.rms = np.sqrt((1.0/float(Np))*tmp2)

# set up 1 million case                 
onemil = []
cycle_start = 11
cycle_end = 110
i = cycle_start
while i <= cycle_end:
	tmp = EigenFunction('output.h5',i)
	onemil.append(tmp)
	i += 1

# set reference
onemil[cycle_end-cycle_start].set_reference()

# loop through and compute rms
i = 0 
while i < len(onemil):
	onemil[i].compute_rms()
	print onemil[i].rms
	i += 1
