#!/usr/bin/env python

import sys

# import statepoint
sys.path.append('../../src/utils')
import statepoint

# read in statepoint file
if len(sys.argv) > 1:
    sp = statepoint.StatePoint(sys.argv[1])
else:
    sp = statepoint.StatePoint('statepoint.8.binary')
sp.read_results()
sp.read_source()

# set up output string
outstr = ''
 
# write out k-combined
outstr += 'k-combined:\n'
outstr += "{0:12.6E} {1:12.6E}\n".format(sp.k_combined[0], sp.k_combined[1])

# write out xyz
xyz = sp.source[0].xyz
for i in xyz:
    outstr += "{0:12.6E} ".format(i)
outstr += "\n"

# write results to file
with open('results_test.dat','w') as fh:
    fh.write(outstr)
