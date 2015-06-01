#!/usr/bin/env python

import sys

sys.path.insert(0, '../..')
from openmc.statepoint import StatePoint

# read in statepoint file
if len(sys.argv) > 1:
    sp = StatePoint(sys.argv[1])
else:
    sp = StatePoint('statepoint.7.binary')

sp.read_results()

# set up output string
outstr = ''

# write out k-combined
outstr += 'k-combined:\n'
outstr += "{0:12.6E} {1:12.6E}\n".format(sp._k_combined[0], sp._k_combined[1])

# write results to file
with open('results_test.dat','w') as fh:
    fh.write(outstr)
