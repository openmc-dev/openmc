#!/usr/bin/env python

import sys
import numpy as np

sys.path.insert(0, '../..')
from openmc.statepoint import StatePoint

# read in statepoint file
if len(sys.argv) > 1:
    sp = StatePoint(sys.argv[1])
else:
    sp = StatePoint('statepoint.07.binary')

sp.read_results()

# extract tally results and convert to vector
tally5 = sp._tallies[5]
results5 = np.zeros((tally5._sum.size + tally5._sum.size, ))
results5[0::2] = tally5._sum.ravel()
results5[1::2] = tally5._sum_sq.ravel()

tally10 = sp._tallies[10]
results10 = np.zeros((tally10._sum.size + tally10._sum.size, ))
results10[0::2] = tally10._sum.ravel()
results10[1::2] = tally10._sum_sq.ravel()

# set up output string
outstr = ''

# write out k-combined
outstr += 'k-combined:\n'
outstr += "{0:12.6E} {1:12.6E}\n".format(sp._k_combined[0], sp._k_combined[1])

# write out tally results
outstr += 'tally 1:\n'
for item in results10:
  outstr += "{0:12.6E}\n".format(item)
outstr += 'tally 2:\n'
for item in results5:
  outstr += "{0:12.6E}\n".format(item)

# write results to file
with open('results_test.dat','w') as fh:
    fh.write(outstr)
