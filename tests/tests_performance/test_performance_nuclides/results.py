#!/usr/bin/env python

import os

with open('stdout.log') as fh:
  lines = fh.readlines()

outstr = ''
for line in lines:
  if 'Total time for initialization' in line:
    words = line.split()
    outstr += 'initialization {0}\n'.format(words[5])
  if 'Total time in simulation' in line:
    words = line.split()
    outstr += 'simulation {0}\n'.format(words[5])

# write results to file
with open('results_test.dat','w') as fh:
    fh.write(outstr)
