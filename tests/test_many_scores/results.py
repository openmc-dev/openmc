#!/usr/bin/env python

import sys
import numpy as np

sys.path.insert(0, '../..')
from openmc.statepoint import StatePoint

# read in statepoint file
if len(sys.argv) > 1:
    sp = StatePoint(sys.argv[1])
else:
    sp = StatePoint('statepoint.5.binary')

sp.read_results()

# extract tally results and convert to vector
tally1 = sp._tallies[1]
results1 = np.zeros((tally1._sum.size + tally1._sum.size, ))
results1[0::2] = tally1._sum.ravel()
results1[1::2] = tally1._sum_sq.ravel()

# set up output string
outstr = ''

# write out k-combined
outstr += 'k-combined:\n'
outstr += "{0:12.6E} {1:12.6E}\n".format(sp._k_combined[0], sp._k_combined[1])

# write out tally results
outstr += 'tally 1:\n'
for item in results1:
  outstr += "{0:12.6E}\n".format(item)

# write results to file
with open('results_test.dat','w') as fh:
    fh.write(outstr)
