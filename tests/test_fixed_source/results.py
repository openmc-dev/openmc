#!/usr/bin/env python

import sys
import numpy as np

sys.path.insert(0, '../..')
from openmc.statepoint import StatePoint

# read in statepoint file
if len(sys.argv) > 1:
    sp = StatePoint(sys.argv[1])
else:
    sp = StatePoint('statepoint.10.binary')

sp.read_results()

# extract tally results and convert to vector
tally = sp._tallies[1]
results = np.zeros((tally._sum.size + tally._sum.size, ))
results[0::2] = tally._sum.ravel()
results[1::2] = tally._sum_sq.ravel()

# set up output string
outstr = ''

# write out tally results
outstr += 'tallies:\n'
for item in results:
  outstr += "{0:12.6E}\n".format(item)

# write results to file
with open('results_test.dat','w') as fh:
    fh.write(outstr)
