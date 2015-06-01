#!/usr/bin/env python

import sys

sys.path.insert(0, '../..')
from openmc.statepoint import StatePoint

# read in statepoint file
if len(sys.argv) > 1:
    sp = StatePoint(sys.argv[1])
else:
    sp = StatePoint('statepoint.08.binary')

sp.read_results()
sp.read_source()

# set up output string
outstr = ''

# write out k-combined
outstr += 'k-combined:\n'
outstr += "{0:12.6E} {1:12.6E}\n".format(sp._k_combined[0], sp._k_combined[1])

# write out xyz
xyz = sp._source[0]._xyz
for i in xyz:
    outstr += "{0:12.6E} ".format(i)
outstr += "\n"

# write results to file
with open('results_test.dat','w') as fh:
    fh.write(outstr)
