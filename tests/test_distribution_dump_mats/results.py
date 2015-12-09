#!/usr/bin/env python

import os
import sys

import h5py

# set up output string
outstr = ''

# get material file contents
matfile = 'materials-out.h5'
file = h5py.File(matfile, 'r')
for i in file.keys():
	group = file[i]
	outstr += '%s\n' % group.name
	outstr += '%s\n' % group['n_nuclides'].value
	outstr += '%s\n' % group['n_instances'].value
	for comp in group['comps']:
		outstr += ' '.join([str(comp)]) + '\n'

# write results to file
with open('results_test.dat','w') as fh:
    fh.write(outstr)
