#!/usr/bin/env python

import sys
import numpy as np

# import statepoint
sys.path.insert(0, '../../src/utils')
import statepoint

# read in statepoint file
if len(sys.argv) > 1:
    sp = statepoint.StatePoint(sys.argv[1])
else:
    sp = statepoint.StatePoint('statepoint.20.binary')
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
results3 = sp.tallies[2].results
shape3 = results3.shape
size3 = (np.product(shape3))
results3 = np.reshape(results3, size3)
results4 = sp.tallies[3].results
shape4 = results4.shape
size4 = (np.product(shape4))
results4 = np.reshape(results4, size4)

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
outstr += 'tally 3:\n'
for item in results3:
    outstr += "{0:12.6E}\n".format(item)
outstr += 'tally 4:\n'
for item in results4:
    outstr += "{0:12.6E}\n".format(item)

# write out cmfd answers
outstr += 'cmfd indices\n'
for item in sp.cmfd_indices:
    outstr += "{0:12.6E}\n".format(item)
outstr += 'k cmfd\n'
for item in sp.k_cmfd:
    outstr += "{0:12.6E}\n".format(item)
outstr += 'cmfd entropy\n'
for item in sp.cmfd_entropy:
    outstr += "{0:12.6E}\n".format(item)
outstr += 'cmfd balance\n'
for item in sp.cmfd_balance:
    outstr += "{0:12.6E}\n".format(item)
outstr += 'cmfd dominance ratio\n'
for item in sp.cmfd_dominance:
    outstr += "{0:10.3E}\n".format(item)
outstr += 'cmfd openmc source comparison\n'
for item in sp.cmfd_srccmp:
    outstr += "{0:12.6E}\n".format(item)
outstr += 'cmfd source\n'
cmfdsrc = np.reshape(sp.cmfd_src, np.product(sp.cmfd_indices), order='F')
for item in cmfdsrc:
    outstr += "{0:12.6E}\n".format(item)

# write results to file
with open('results_test.dat', 'w') as fh:
    fh.write(outstr)
