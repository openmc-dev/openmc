#!/usr/bin/python
# Filename: eigenfunction_rms.py

import h5py
import numpy as np
import cPickle
import matplotlib.pyplot as plt
import sys
#
class EigenFunction:
#
 	'''Represents the reference eigenfunction'''
 	reference = 0.
#
	def __init__(self):
#
		'''Initializes the eigenfunction'''
		self.function = 0.
#
        def read_hdf5(self,h5_file,cycle):
#
		'''Read data from HDF5 file'''
                self.cycle = cycle
                f = h5py.File(h5_file,'r')
                group = '/cycle'+str(cycle)+'/openmc_src'
                dataset = f[group]
                self.function = np.empty(dataset.shape,dataset.dtype)
                dataset.read_direct(self.function)
                self.iamref = 'F'
#
	def set_reference(self):
#
		'''Sets instance to be reference calc'''
		self.iamref = 'T'
		EigenFunction.reference = self.function
#
	def compute_rms(self):
#
		'''Computes RMS value'''
		Np = self.function.size
                tmp = (self.function - EigenFunction.reference)*(self.function - EigenFunction.reference)
                tmp2 = tmp.sum()
                self.rms = np.sqrt((1.0/float(Np))*tmp2)
#
if sys.argv[1] != 'restart':

	# calculate reference solution
	print 'Calculating Reference solution...'
	runpath = '/media/Backup/opr_runs/1mil/run' # the directory prefix
	hdfile  = '/output.h5'  # hdf5 file name
	cycle = 840    # cycle number to extract
	run_end = 25   # number of runs to look at
	tmp = EigenFunction()
	tmp.read_hdf5(runpath+str(1)+hdfile,cycle) # load first eigenfunction
	indices = tmp.function.shape # extent of all dimensions
	ref = np.zeros((run_end,indices[0],indices[1],indices[2],indices[3])) # initialize ref array
	ref[0] = tmp.function # set the first run in ref
	i = 2
	while i <= run_end: # begin loop around all runs
		tmp.read_hdf5(runpath+str(i)+hdfile,cycle)
		print 'Read in: '+runpath+str(i)+hdfile
       		ref[i-1] = tmp.function
		i += 1
	meanref = np.average(ref, axis=0) # compute average of all runs
	EigenFunction.reference = meanref # set to global space in EigenFunction instances

	# calculate rms for 1 million case
	print 'Calculating 1million History case...'
	onemil = []
	cycle_start = 201
	cycle_end = 840
	run_end = 10
	i = cycle_start
	while i <= cycle_end:
		meantmp = EigenFunction() # init mean object
		j = 1
		runs = np.zeros((run_end,indices[0],indices[1],indices[2],indices[3])) # init runs array
 		while j <= run_end:
			tmp.read_hdf5(runpath+str(j)+hdfile,i) # read hdf5 file
			runs[j-1] = tmp.function # put function into runs
			j += 1
		meantmp.function = np.average(runs, axis=0) # compute the mean
		onemil.append(meantmp)
		print 'Read in from path: '+runpath+' Cycle: '+str(i)
		i += 1

	# calculate rms array
	print 'Calculating rms...'
	rms = np.zeros((cycle_end - cycle_start + 1))
	i = 0
	while i < len(onemil):
		onemil[i].compute_rms()
		rms[i] = onemil[i].rms
		i += 1

	# write out numpy array to binary file
	print 'Writing output...'
	output = {}
	output.update({'1milrms':rms})
	output.update({'ref':meanref})
	fileout = open('rms.out','wb')
	cPickle.dump(output,fileout)
	fileout.close()
else:
	# load in data
	print 'Loading input...'
	filein = open('rms.out','r')
	output = cPickle.load(filein)
	filein.close()
	rms = output['1milrms']
	meanref = output['ref']
	EigenFunction.reference = meanref

# plot rms
print 'Generating plot...'
ax = plt.subplot(111)
x = np.linspace(1,640,640)*1e6
y = rms[0]/1e-3*x**(-0.5)
plt.loglog(x,rms*100,'b--')
plt.loglog(x,y*100,'g-')
ax.xaxis.grid(True,'minor')
ax.yaxis.grid(True,'minor')
ax.xaxis.grid(True,'major',linewidth=2)
ax.yaxis.grid(True,'major',linewidth=2)
plt.xlabel('# of Total Neutron Histories (active cycles)')
plt.ylabel('RMS Error [%]')
plt.legend(('1 million (10 runs)','Ideal Error 1 mil'))

# plot mean source distribution
plt.figure()
X = np.linspace(1,272,272) 
Y = np.linspace(1,272,272)
Y,X = np.meshgrid(Y,X)
plt.contourf(X,Y,meanref[0,:,:,0],100)
plt.colorbar()
plt.xlabel('Mesh Cell in x-direction')
plt.ylabel('Mesh Cell in y-direction')
plt.title('OPR Converged Fission Source Distribution')
plt.show()
