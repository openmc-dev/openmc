#!/usr/bin/env python

import os
from glob import glob
import sys
from subprocess import call
import h5py as h5

n = len(os.listdir('.'))
ndirs = 0
for i in range(0,n):
  if os.path.isdir(os.listdir('.')[i]):
    ndirs += 1

ktest = []
sigtest = []
for i in range(1,ndirs+1):
  os.chdir('./case-'+str(i))
  statefiles = glob('statepoint*')
  errfiles = glob('error.log')
  nfiles = len(statefiles) + len(errfiles)
  if nfiles == 1:
    if len(statefiles) == 1:
      hf = h5.File(statefiles[0], 'r')
      ktest.append(hf['k_combined'][0])
      sigtest.append(hf['k_combined'][1])
    else:
      file = open(errfiles[0], 'r')
      for line in file:
        ktest.append(line)
        break
      file.close()
  elif nfiles == 0:
    sys.exit('No output files for case-'+str(i))
  elif nfiles > 1:
    sys.exit('Multiple output files for case-'+str(i))
  os.chdir('..')

testfile = open('testSystemResults.dat', 'w')
for kval in ktest:
  if isinstance(kval, basestring):
    testfile.write(kval)
  else:
    testfile.write('%.17f'%float(kval)+'\n')
testfile.close()

refk = []
testk = []
reffile = open('referenceSystemResults.dat', 'r')
testfile = open('testSystemResults.dat', 'r')
for line in reffile:
  refk.append(line)
for line in testfile:
  testk.append(line)
reffile.close()
testfile.close()

cnt = 0
for kval in refk:
  if kval == testk[cnt]:
    call(['echo', 'case-'+str(cnt+1)+' passed'])
  else:
    call(['echo', 'case-'+str(cnt+1)+' failed'])
  cnt += 1
