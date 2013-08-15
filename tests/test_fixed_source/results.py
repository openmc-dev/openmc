#!/usr/bin/env python

import sys
import numpy as np

# import statepoint
sys.path.append('../../src/utils')
import statepoint

# read in statepoint file
sp = statepoint.StatePoint('statepoint.10.binary')
sp.read_results()

# extract tally results and convert to vector
results = sp.tallies[0].results
shape = results.shape
size = (np.product(shape))
results = np.reshape(results, size)

# set up output string
outstr = ''
 
# write out tally results
outstr += 'tallies:\n'
for item in results:
  outstr += "{0:10.8f}\n".format(item)

# write results to file
with open('results_test.dat','w') as fh:
    fh.write(outstr)
