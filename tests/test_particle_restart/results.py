#!/usr/bin/env python

import sys

# import particle restart 
sys.path.append('../../src/utils')
import particle_restart as pr

# read in particle restart file
if len(sys.argv) > 1:
    p = pr.Particle(sys.argv[1])
else:
    p = pr.Particle('particle_12_192.binary')

# set up output string
outstr = ''
 
# write out properties 
outstr += 'current batch:\n'
outstr += "{0:12.6E}\n".format(p.current_batch)
outstr += 'current gen:\n'
outstr += "{0:12.6E}\n".format(p.current_gen)
outstr += 'particle id:\n'
outstr += "{0:12.6E}\n".format(p.id)
outstr += 'particle weight:\n'
outstr += "{0:12.6E}\n".format(p.weight)
outstr += 'particle energy:\n'
outstr += "{0:12.6E}\n".format(p.energy)
outstr += 'particle xyz:\n'
outstr += "{0:12.6E} {1:12.6E} {2:12.6E}\n".format(p.xyz[0],p.xyz[1],p.xyz[2])
outstr += 'particle uvw:\n'
outstr += "{0:12.6E} {1:12.6E} {2:12.6E}\n".format(p.uvw[0],p.uvw[1],p.uvw[2])

# write results to file
with open('results_test.dat','w') as fh:
    fh.write(outstr)
