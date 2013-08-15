#!/usr/bin/env python

import sys

# import statepoint
sys.path.append('../../src/utils')
import statepoint

# read in statepoint file
sp = statepoint.StatePoint('statepoint.10.binary')
sp.read_results()

# set up output string
outstr = ''
 
# write out k-combined
outstr += 'k-combined:\n'
outstr += "{0:10.8f} {1:10.8f}\n".format(sp.k_combined[0], sp.k_combined[1])

# write results to file
with open('results_test.dat','w') as fh:
    fh.write(outstr)
