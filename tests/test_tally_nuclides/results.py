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
results0 = sp.tallies[0].results
results1 = sp.tallies[1].results
results = np.concatenate((results0.flatten(), results1.flatten()))

# set up output string
outstr = ''
 
# write out k-combined
outstr += 'k-combined:\n'
outstr += "{0:12.6E} {1:12.6E}\n".format(sp.k_combined[0], sp.k_combined[1])

# write out tally results
outstr += 'tallies:\n'
for item in results:
  outstr += "{0:12.6E}\n".format(item)

# write results to file
with open('results_test.dat','w') as fh:
    fh.write(outstr)
