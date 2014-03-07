#!/usr/bin/env python

import sys

# import statepoint
sys.path.append('../../src/utils')
import statepoint

# read in statepoint file
sp3 = statepoint.StatePoint(sys.argv[1])
sp2 = statepoint.StatePoint(sys.argv[2])
sp1 = statepoint.StatePoint(sys.argv[3])
sp1.read_results()
sp2.read_results()
sp3.read_results()

# analyze sp1
# compare distrib sum to cell
sp1_t1_c1 = sp1.get_value(0,[('distribcell',0)],0)[0]
sp1_t1_c2 = sp1.get_value(0,[('distribcell',1)],0)[0]
sp1_t1_c3 = sp1.get_value(0,[('distribcell',2)],0)[0]
sp1_t1_c4 = sp1.get_value(0,[('distribcell',3)],0)[0]
sp1_t1 = sp1_t1_c1 + sp1_t1_c2 + sp1_t1_c3 + sp1_t1_c4
sp1_t2 =  sp1.get_value(1,[('cell',0)],0)[0]
   
# analyze sp2
sp2_t1_c1 =  sp2.get_value(0,[('cell',0)],0)[0]
sp2_t1_c2 =  sp2.get_value(0,[('cell',1)],0)[0]
sp2_t1_c3 =  sp2.get_value(0,[('cell',2)],0)[0]
sp2_t1_c4 =  sp2.get_value(0,[('cell',3)],0)[0]

# analyze sp3
sp3_t1 = sp3.get_value(0,[('cell',0)],0)[0]
sp3_t2 = sp3.get_value(1,[('distribcell',0)],0)[0]

sp3_t3 = sp3.get_value(2,[('cell',0)],0)[0]
sp3_t4 = 0
for i in range(241):
  sp3_t4 += sp3.get_value(3,[('distribcell',i)],0)[0]

sp3_t5 = sp3.get_value(4,[('cell',0)],0)[0]
sp3_t6 = 0
for i in range(63624):
  sp3_t6 += sp3.get_value(5,[('distribcell',i)],0)[0]



# What if we don't know the index in the array?
# Well, we know how to get there. Explained:
# Universe 0 -> Cell 1 -> Lattice 200 (x,y,z = 4,6,1) -> Universe 6 -> Cell 60
path = [0,1,(200,4,6,1),6,60]

# So we have the path to our cell now.
# Lets get that specific value
sp3_t4_c60 = sp3.get_value(3,[('distribcell',path)],0)[0]

# set up output string
outstr = ''
 
outstr += "{0:12.6E}\n".format(sp1_t1)
outstr += "{0:12.6E}\n".format(sp1_t2)
outstr += "{0:12.6E}\n".format(sp2_t1_c1)
outstr += "{0:12.6E}\n".format(sp1_t1_c4)
outstr += "{0:12.6E}\n".format(sp2_t1_c2)
outstr += "{0:12.6E}\n".format(sp1_t1_c1)
outstr += "{0:12.6E}\n".format(sp2_t1_c3)
outstr += "{0:12.6E}\n".format(sp1_t1_c3)
outstr += "{0:12.6E}\n".format(sp2_t1_c4)
outstr += "{0:12.6E}\n".format(sp1_t1_c2)
outstr += "{0:12.6E}\n".format(sp3_t1)
outstr += "{0:12.6E}\n".format(sp3_t2)
outstr += "{0:12.6E}\n".format(sp3_t3)
outstr += "{0:12.6E}\n".format(sp3_t4)
outstr += "{0:12.6E}\n".format(sp3_t5)
outstr += "{0:12.6E}\n".format(sp3_t6)
outstr += "{0:12.6E}\n".format(sp3_t4_c60)

# write results to file
with open('results_test.dat','w') as fh:
    fh.write(outstr)
