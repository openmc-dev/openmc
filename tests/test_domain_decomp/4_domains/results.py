#!/usr/bin/env python

import os
import sys
import numpy as np

# import statepoint
sys.path.insert(0, '../../../src/utils')
import statepoint

# read in statepoint files

spfile = 'statepoint.20.domain_1.binary'
if not os.path.exists(spfile): spfile = 'statepoint.20.domain_1.h5'
sp1 = statepoint.StatePoint(spfile)
spfile = 'statepoint.20.domain_2.binary'
if not os.path.exists(spfile): spfile = 'statepoint.20.domain_2.h5'
sp2 = statepoint.StatePoint(spfile)
spfile = 'statepoint.20.domain_3.binary'
if not os.path.exists(spfile): spfile = 'statepoint.20.domain_3.h5'
sp3 = statepoint.StatePoint(spfile)
spfile = 'statepoint.20.domain_4.binary'
if not os.path.exists(spfile): spfile = 'statepoint.20.domain_4.h5'
sp4 = statepoint.StatePoint(spfile)
sp1.read_results()
sp2.read_results()
sp3.read_results()
sp4.read_results()

# extract tally results, means only since sum_sq won't match the 1_domain case
results = [sp1.tallies[0].results[:,:,0],
           sp2.tallies[0].results[:,:,0],
           sp3.tallies[0].results[:,:,0],
           sp4.tallies[0].results[:,:,0]]

# combine results for bins that were on more than one domain, in real_bin order
maps = np.array([sp1.tallies[0].otf_filter_bin_map,
                 sp2.tallies[0].otf_filter_bin_map,
                 sp3.tallies[0].otf_filter_bin_map,
                 sp4.tallies[0].otf_filter_bin_map])
maxbin = max(maps.flatten())
tmp_results = np.zeros((maxbin, 1))
for i in range(maxbin):
    for j, map_ in enumerate(maps):
      if i+1 in map_:
          bin_ = np.nonzero(map_ == i+1)[0][0]
          tmp_results[i] += results[j][bin_]
results = np.reshape(tmp_results, (np.product(tmp_results.shape)))

# set up output string
outstr = ''
 
# write out k-combined
outstr += 'k-combined:\n'
outstr += "{0:12.6E} {1:12.6E}\n".format(sp1.k_combined[0], sp1.k_combined[1])

# write out tally results
outstr += 'tallies:\n'
for item in results:
  outstr += "{0:12.6E}\n".format(item)

# write results to file
with open('results_test.dat','w') as fh:
    fh.write(outstr)
