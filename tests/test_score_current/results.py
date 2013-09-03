#!/usr/bin/env python

import sys
import numpy as np

# import statepoint
sys.path.append('../../src/utils')
import statepoint

# read in statepoint file
if len(sys.argv) > 1:
    sp = statepoint.StatePoint(sys.argv[1])
else:
    sp = statepoint.StatePoint('statepoint.10.binary')
sp.read_results()

# extract tally results and convert to vector
results1 = sp.tallies[0].results
shape1 = results1.shape
size1 = (np.product(shape1))
results1 = np.reshape(results1, size1)
results2 = sp.tallies[1].results
shape2 = results2.shape
size2 = (np.product(shape2))
results2 = np.reshape(results2, size2)

# set up output string
outstr = ''
 
# write out k-combined
outstr += 'k-combined:\n'
outstr += "{0:12.6E} {1:12.6E}\n".format(sp.k_combined[0], sp.k_combined[1])

# write out tally results
outstr += 'tally 1:\n'
for item in results1:
  outstr += "{0:12.6E}\n".format(item)
outstr += 'tally 2:\n'
for item in results2:
  outstr += "{0:12.6E}\n".format(item)

# write results to file
with open('results_test.dat','w') as fh:
    fh.write(outstr)
