#!/usr/bin/env python

import sys

# import particle restart 
sys.path.append('../../src/utils')
import particle_restart as pr

# read in particle restart file
p = pr.Particle('particle_10_638.binary')

# set up output string
outstr = ''
 
# write out properties 
outstr += 'current batch:\n'
outstr += "{0:10.8f}\n".format(p.current_batch)
outstr += 'current gen:\n'
outstr += "{0:10.8f}\n".format(p.current_gen)
outstr += 'particle id:\n'
outstr += "{0:10.8f}\n".format(p.id)
outstr += 'particle weight:\n'
outstr += "{0:10.8f}\n".format(p.weight)
outstr += 'particle energy:\n'
outstr += "{0:10.8f}\n".format(p.energy)
outstr += 'particle xyz:\n'
outstr += "{0:10.8f} {1:10.8f} {2:10.8f}\n".format(p.xyz[0],p.xyz[1],p.xyz[2])
outstr += 'particle uvw:\n'
outstr += "{0:10.8f} {1:10.8f} {2:10.8f}\n".format(p.uvw[0],p.uvw[1],p.uvw[2])

# write results to file
with open('results_test.dat','w') as fh:
    fh.write(outstr)
