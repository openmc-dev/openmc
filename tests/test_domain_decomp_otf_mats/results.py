#!/usr/bin/env python

import os
import sys

# import dumpmat
sys.path.insert(0, '../../src/utils')
import dumpmat

# set up output string
outstr = ''

# get list of files in order
matfiles = []
matfile = 'material.m1.binary'
if not os.path.exists(matfile): matfile = 'material.m1.h5'
matfiles.append(matfile)
matfile = 'material.m2.binary'
if not os.path.exists(matfile): matfile = 'material.m2.h5'
matfiles.append(matfile)
matfile = 'material.m3.binary'
if not os.path.exists(matfile): matfile = 'material.m3.h5'
matfiles.append(matfile)

# get material file contents
for matfile in matfiles:
  outstr += '%s\n' % dumpmat.get_n_nucs(matfile)
  outstr += '%s\n' % dumpmat.get_n_inst(matfile)
  for comp in dumpmat.get_comps(matfile):
    outstr += ' '.join([str(c) for c in comp]) + '\n'

# write results to file
with open('results_test.dat','w') as fh:
    fh.write(outstr)
