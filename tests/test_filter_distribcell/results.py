#!/usr/bin/env python

import sys
import numpy as np

sys.path.insert(0, '../..')
from openmc.statepoint import StatePoint
from openmc import Filter

# read in statepoint file
sp3 = StatePoint(sys.argv[1])
sp2 = StatePoint(sys.argv[2])
sp1 = StatePoint(sys.argv[3])
sp1.read_results()
sp2.read_results()
sp3.read_results()
sp1.compute_stdev()
sp2.compute_stdev()
sp3.compute_stdev()

# analyze sp1
# compare distrib sum to cell
sp1_t1 = sp1.tallies[1]
sp1_t2 = sp1.tallies[2]

sp1_t1_d1 = sp1_t1.get_values(scores=['total'], filters=['distribcell'],
                              filter_bins=[(0,)], value='mean')
sp1_t1_d2 = sp1_t1.get_values(scores=['total'], filters=['distribcell'],
                              filter_bins=[(1,)], value='mean')
sp1_t1_d3 = sp1_t1.get_values(scores=['total'], filters=['distribcell'],
                              filter_bins=[(2,)], value='mean')
sp1_t1_d4 = sp1_t1.get_values(scores=['total'], filters=['distribcell'],
                              filter_bins=[(3,)], value='mean')
sp1_t1_sum = sp1_t1_d1 + sp1_t1_d2 + sp1_t1_d3 + sp1_t1_d4

cell_filter = sp1_t2.find_filter(filter_type='cell')
sp1_t2_c201 = sp1_t2.get_values(scores=['total'], value='mean')


# analyze sp2
sp2_t1 = sp2.tallies[1]
sp2_t1_c201 = sp2_t1.get_values(scores=['total'], filters=['cell'],
                                filter_bins=[(2,)], value='mean')
sp2_t1_c203 = sp2_t1.get_values(scores=['total'], filters=['cell'],
                                filter_bins=[(4,)], value='mean')
sp2_t1_c205 = sp2_t1.get_values(scores=['total'], filters=['cell'],
                                filter_bins=[(6,)], value='mean')
sp2_t1_c207 = sp2_t1.get_values(scores=['total'], filters=['cell'],
                                filter_bins=[(8,)], value='mean')


# analyze sp3
sp3_t1 = sp3.tallies[1]
sp3_t2 = sp3.tallies[2]
sp3_t3 = sp3.tallies[3]
sp3_t4 = sp3.tallies[4]
sp3_t5 = sp3.tallies[5]
sp3_t6 = sp3.tallies[6]

sp3_t1_c1 = sp3_t1.get_values(scores=['total'], filters=['cell'],
                              filter_bins=[(1,)], value='mean')
sp3_t2_d1 = sp3_t2.get_values(scores=['total'], filters=['distribcell'],
                              filter_bins=[(0,)], value='mean')

sp3_t3_c60 = sp3_t3.get_values(scores=['total'], filters=['cell'],
                               filter_bins=[(26,)], value='mean')
sp3_t4_d = sum(sp3_t4.get_values(scores=['total'], value='mean'))

sp3_t5_c27 = sp3_t5.get_values(scores=['total'], filters=['cell'],
                               filter_bins=[(19,)], value='mean')
sp3_t6_d = sum(sp3_t6.get_values(scores=['total'], value='mean'))

# set up output string
outstr = ''

outstr += "{0:12.6E}\n".format(sp1_t1_sum)
outstr += "{0:12.6E}\n".format(sp1_t2_c201[()])
outstr += "{0:12.6E}\n".format(sp2_t1_c201[()])
outstr += "{0:12.6E}\n".format(sp1_t1_d4[()])
outstr += "{0:12.6E}\n".format(sp2_t1_c203[()])
outstr += "{0:12.6E}\n".format(sp1_t1_d1[()])
outstr += "{0:12.6E}\n".format(sp2_t1_c205[()])
outstr += "{0:12.6E}\n".format(sp1_t1_d3[()])
outstr += "{0:12.6E}\n".format(sp2_t1_c207[()])
outstr += "{0:12.6E}\n".format(sp1_t1_d2[()])
outstr += "{0:12.6E}\n".format(sp3_t1_c1[()])
outstr += "{0:12.6E}\n".format(sp3_t2_d1[()])
outstr += "{0:12.6E}\n".format(sp3_t3_c60[()])
outstr += "{0:12.6E}\n".format(sp3_t4_d)
outstr += "{0:12.6E}\n".format(sp3_t5_c27[()])
outstr += "{0:12.6E}\n".format(sp3_t6_d)

# write results to file
with open('results_test.dat','w') as fh:
    fh.write(outstr)
