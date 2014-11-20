#!/usr/bin/env python

import os
import sys

# import dumpmat
sys.path.insert(0, '../../src/utils')
import statepoint

spfile = 'statepoint.20.domain_1.binary'
if not os.path.exists(spfile): spfile = 'statepoint.20.domain_1.h5'
sp = statepoint.StatePoint(spfile)
sp.read_domain_tally_metadata()
sp.read_results()

# set up output string
outstr = ''

spec_list = [('distribcell', [0, 1, (6,1,1,1), 5, 2, (4,1,2,1), 1, 101])]
outstr += '%s\n' % sp.get_value(1, spec_list, 0)[0]
spec_list = [('distribcell', [0, 1, (6,1,2,1), 5, 2, (4,1,2,1), 1, 101])]
outstr += '%s\n' % sp.get_value(1, spec_list, 0)[0]
spec_list = [('distribcell', [0, 1, (6,2,1,1), 5, 2, (4,1,2,1), 1, 101])]
outstr += '%s\n' % sp.get_value(1, spec_list, 0)[0]
spec_list = [('distribcell', [0, 1, (6,2,2,1), 5, 2, (4,1,2,1), 1, 101])]
outstr += '%s\n' % sp.get_value(1, spec_list, 0)[0]

spec_list = [('distribcell', [0, 1, (6,1,1,1), 5, 2, (4,1,1,1), 2, 201])]
outstr += '%s\n' % sp.get_value(2, spec_list, 0)[0]
spec_list = [('distribcell', [0, 1, (6,1,1,1), 5, 2, (4,2,2,1), 2, 201])]
outstr += '%s\n' % sp.get_value(2, spec_list, 0)[0]
spec_list = [('distribcell', [0, 1, (6,2,1,1), 5, 2, (4,1,1,1), 2, 201])]
outstr += '%s\n' % sp.get_value(2, spec_list, 0)[0]
spec_list = [('distribcell', [0, 1, (6,2,1,1), 5, 2, (4,2,2,1), 2, 201])]
outstr += '%s\n' % sp.get_value(2, spec_list, 0)[0]
spec_list = [('distribcell', [0, 1, (6,1,2,1), 5, 2, (4,1,1,1), 2, 201])]
outstr += '%s\n' % sp.get_value(2, spec_list, 0)[0]
spec_list = [('distribcell', [0, 1, (6,1,2,1), 5, 2, (4,2,2,1), 2, 201])]
outstr += '%s\n' % sp.get_value(2, spec_list, 0)[0]
spec_list = [('distribcell', [0, 1, (6,2,2,1), 5, 2, (4,1,1,1), 2, 201])]
outstr += '%s\n' % sp.get_value(2, spec_list, 0)[0]
spec_list = [('distribcell', [0, 1, (6,2,2,1), 5, 2, (4,2,2,1), 2, 201])]
outstr += '%s\n' % sp.get_value(2, spec_list, 0)[0]

spec_list = [('distribcell', [0, 1, (6,1,1,1), 5, 2, (4,2,1,1), 3, 301])]
outstr += '%s\n' % sp.get_value(3, spec_list, 0)[0]
spec_list = [('distribcell', [0, 1, (6,1,2,1), 5, 2, (4,2,1,1), 3, 301])]
outstr += '%s\n' % sp.get_value(3, spec_list, 0)[0]
spec_list = [('distribcell', [0, 1, (6,2,1,1), 5, 2, (4,2,1,1), 3, 301])]
outstr += '%s\n' % sp.get_value(3, spec_list, 0)[0]
spec_list = [('distribcell', [0, 1, (6,2,2,1), 5, 2, (4,2,1,1), 3, 301])]
outstr += '%s\n' % sp.get_value(3, spec_list, 0)[0]

# write results to file
with open('results_test.dat','w') as fh:
    fh.write(outstr)
