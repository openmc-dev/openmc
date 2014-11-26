#!/usr/bin/env python

import sys
import numpy as np

# import statepoint
from openmc.statepoint import StatePoint

# read in statepoint file
if len(sys.argv) > 1:
    sp = StatePoint(sys.argv[1])
else:
    sp = StatePoint('statepoint.20.binary')

sp.read_results()

# extract tally results and convert to vector
tally1 = sp._tallies[1]
results1 = np.zeros((tally1._sum.size + tally1._sum.size, ))
results1[0::2] = tally1._sum.ravel()
results1[1::2] = tally1._sum_sq.ravel()

for tally_id in sp._tallies:
  if 'CMFD flux, total, scatter-1' in sp._tallies[tally_id]._label:
    tally2 = sp._tallies[tally_id]
  elif 'CMFD neutron production' in sp._tallies[tally_id]._label:
    tally3 = sp._tallies[tally_id]
  elif 'CMFD surface currents' in sp._tallies[tally_id]._label:
    tally4 = sp._tallies[tally_id]

results2 = np.zeros((tally2._sum.size + tally2._sum.size, ))
results2[0::2] = tally2._sum.ravel()
results2[1::2] = tally2._sum_sq.ravel()

results3 = np.zeros((tally3._sum.size + tally3._sum.size, ))
results3[0::2] = tally3._sum.ravel()
results3[1::2] = tally3._sum_sq.ravel()

results4 = np.zeros((tally4._sum.size + tally4._sum.size, ))
results4[0::2] = tally4._sum.ravel()
results4[1::2] = tally4._sum_sq.ravel()

# set up output string
outstr = ''

# write out k-combined
outstr += 'k-combined:\n'
outstr += "{0:12.6E} {1:12.6E}\n".format(sp._k_combined[0], sp._k_combined[1])

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
for item in sp._cmfd_indices:
    outstr += "{0:12.6E}\n".format(item)
outstr += 'k cmfd\n'
for item in sp._k_cmfd:
    outstr += "{0:12.6E}\n".format(item)
outstr += 'cmfd entropy\n'
for item in sp._cmfd_entropy:
    outstr += "{0:12.6E}\n".format(item)
outstr += 'cmfd balance\n'
for item in sp._cmfd_balance:
    outstr += "{0:12.6E}\n".format(item)
outstr += 'cmfd dominance ratio\n'
for item in sp._cmfd_dominance:
    outstr += "{0:10.3E}\n".format(item)
outstr += 'cmfd openmc source comparison\n'
for item in sp._cmfd_srccmp:
    outstr += "{0:12.6E}\n".format(item)
outstr += 'cmfd source\n'
cmfdsrc = np.reshape(sp._cmfd_src, np.product(sp._cmfd_indices), order='F')
for item in cmfdsrc:
    outstr += "{0:12.6E}\n".format(item)

# write results to file
with open('results_test.dat', 'w') as fh:
    fh.write(outstr)
