#!/usr/bin/env python

# import statepoint
from openmc.statepoint import StatePoint
from openmc import Filter
import numpy as np
import sys

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
sp1_t1 = sp1._tallies[1]
sp1_t2 = sp1._tallies[2]

distribcell_filter = sp1_t1.find_filter(filter_type='distribcell', bins=[2])
sp1_t1_d1 = sp1_t1.get_value(score='total', filters=[distribcell_filter],
                            filter_bins=[0], value='mean')
sp1_t1_d2 = sp1_t1.get_value(score='total', filters=[distribcell_filter],
                            filter_bins=[1], value='mean')
sp1_t1_d3 = sp1_t1.get_value(score='total', filters=[distribcell_filter],
                            filter_bins=[2], value='mean')
sp1_t1_d4 = sp1_t1.get_value(score='total', filters=[distribcell_filter],
                            filter_bins=[3], value='mean')
sp1_t1_sum = sp1_t1_d1 + sp1_t1_d2 + sp1_t1_d3 + sp1_t1_d4

cell_filter = sp1_t2.find_filter(filter_type='cell', bins=[2])
sp1_t2_c201 = sp1_t2.get_value(score='total', filters=[cell_filter],
                              filter_bins=[2], value='mean')


# analyze sp2
sp2_t1 = sp2._tallies[1]
cell_filter = sp2_t1.find_filter(filter_type='cell', bins=[2,4,6,8])
sp2_t1_c201 = sp2_t1.get_value(score='total', filters=[cell_filter],
                               filter_bins=[2], value='mean')
sp2_t1_c203 = sp2_t1.get_value(score='total', filters=[cell_filter],
                               filter_bins=[4], value='mean')
sp2_t1_c205 = sp2_t1.get_value(score='total', filters=[cell_filter],
                               filter_bins=[6], value='mean')
sp2_t1_c207 = sp2_t1.get_value(score='total', filters=[cell_filter],
                               filter_bins=[8], value='mean')


# analyze sp3
sp3_t1 = sp3._tallies[1]
sp3_t2 = sp3._tallies[2]
sp3_t3 = sp3._tallies[3]
sp3_t4 = sp3._tallies[4]
sp3_t5 = sp3._tallies[5]
sp3_t6 = sp3._tallies[6]

cell_filter = sp3_t1.find_filter(filter_type='cell', bins=[1])
distribcell_filter = sp3_t2.find_filter(filter_type='distribcell', bins=[1])
sp3_t1_c1 = sp3_t1.get_value(score='total', filters=[cell_filter],
                            filter_bins=[1], value='mean')
sp3_t2_d1 = sp3_t2.get_value(score='total', filters=[distribcell_filter],
                            filter_bins=[0], value='mean')

cell_filter = sp3_t3.find_filter(filter_type='cell', bins=[26])
distribcell_filter = sp3_t4.find_filter(filter_type='distribcell', bins=[26])
sp3_t3_c60 = sp3_t3.get_value(score='total', filters=[cell_filter],
                              filter_bins=[26], value='mean')
sp3_t4_d = 0
for i in range(241):
  sp3_t4_d += sp3_t4.get_value(score='total', filters=[distribcell_filter],
                              filter_bins=[i], value='mean')

cell_filter = sp3_t5.find_filter(filter_type='cell', bins=[19])
distribcell_filter = sp3_t6.find_filter(filter_type='distribcell', bins=[19])
sp3_t5_c27 = sp3_t5.get_value(score='total', filters=[cell_filter],
                              filter_bins=[19], value='mean')
sp3_t6_d = 0
for i in range(63624):
  sp3_t6_d += sp3_t6.get_value(score='total', filters=[distribcell_filter],
                              filter_bins=[i], value='mean')

# set up output string
outstr = ''

outstr += "{0:12.6E}\n".format(sp1_t1_sum)
outstr += "{0:12.6E}\n".format(sp1_t2_c201)
outstr += "{0:12.6E}\n".format(sp2_t1_c201)
outstr += "{0:12.6E}\n".format(sp1_t1_d4)
outstr += "{0:12.6E}\n".format(sp2_t1_c203)
outstr += "{0:12.6E}\n".format(sp1_t1_d1)
outstr += "{0:12.6E}\n".format(sp2_t1_c205)
outstr += "{0:12.6E}\n".format(sp1_t1_d3)
outstr += "{0:12.6E}\n".format(sp2_t1_c207)
outstr += "{0:12.6E}\n".format(sp1_t1_d2)
outstr += "{0:12.6E}\n".format(sp3_t1_c1)
outstr += "{0:12.6E}\n".format(sp3_t2_d1)
outstr += "{0:12.6E}\n".format(sp3_t3_c60)
outstr += "{0:12.6E}\n".format(sp3_t4_d)
outstr += "{0:12.6E}\n".format(sp3_t5_c27)
outstr += "{0:12.6E}\n".format(sp3_t6_d)

# write results to file
with open('results_test.dat','w') as fh:
    fh.write(outstr)
