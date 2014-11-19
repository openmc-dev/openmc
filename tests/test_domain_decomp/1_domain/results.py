#!/usr/bin/env python

import os
import sys
import numpy as np

# import statepoint
sys.path.insert(0, '../../../src/utils')
import statepoint

def order_by(arr, ordering):
    ordered = np.zeros(arr.shape)
    for i,val in enumerate(arr):
        ordered[ordering[i]-1] = val
    return ordered 

spfile = 'statepoint.20.domain_1.binary'
if not os.path.exists(spfile): spfile = 'statepoint.20.domain_1.h5'
sp = statepoint.StatePoint(spfile)
sp.read_results()

# extract tally results (means only) and convert to vector
results = sp.tallies[0].results[:,:,0]
results = order_by(results, sp.tallies[0].otf_filter_bin_map)
shape = results.shape
size = (np.product(shape))
results = np.reshape(results, size)

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
