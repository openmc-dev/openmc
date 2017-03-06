#!/usr/bin/env python

import os
from glob import glob
import sys
import h5py as h5

n = len(os.listdir('.'))
ndirs = 0
for i in range(0,n):
  if os.path.isdir(os.listdir('.')[i]):
    ndirs += 1

k = []
sig = []
for i in range(1,ndirs+1):
  os.chdir('./case-'+str(i))
  statefiles = glob('statepoint*')
  errfiles = glob('error.log')
  nfiles = len(statefiles) + len(errfiles)
  if nfiles == 1:
    if len(statefiles) == 1:
      hf = h5.File(statefiles[0], 'r')
      k.append(hf['k_combined'][0])
      sig.append(hf['k_combined'][1])
    else:
      file = open(errfiles[0], 'r')
      for line in file:
        k.append(line)
        break
      file.close()
  elif nfiles > 1:
    sys.exit('Multiple output files for case-'+str(i))
  elif nfiles == 0:
    sys.exit('No output files for case-'+str(i))
  os.chdir('..')

outfile = open('referenceSystemResults.dat', 'w')
for kval in k:
  if isinstance(kval, basestring):
    outfile.write(kval)
  else:
    outfile.write('%.17f'%float(kval)+'\n')
outfile.close()
