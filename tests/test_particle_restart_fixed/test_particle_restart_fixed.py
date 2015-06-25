#!/usr/bin/env python

import sys
sys.path.insert(0, '..')
from testing_harness import *

import openmc.particle_restart as pr


class ParticleRestartTestHarness(TestHarness):
    def _test_output_created(self):
        """Make sure the restart file has been created."""
        particle = glob.glob(os.path.join(os.getcwd(), self._sp_name))
        assert len(particle) == 1, 'Either multiple or no particle restart ' \
             'files exist.'
        assert particle[0].endswith('binary') \
             or particle[0].endswith('h5'), \
             'Particle restart file is not a binary or hdf5 file.'


    def _get_results(self):
        """Digest info in the restart file and create a simpler ASCII file."""
        # Read the particle restart file.
        particle = glob.glob(os.path.join(os.getcwd(), self._sp_name))[0]
        p = pr.Particle(particle)

        # Write out the properties.
        outstr = ''
        outstr += 'current batch:\n'
        outstr += "{0:12.6E}\n".format(p.current_batch)
        outstr += 'current gen:\n'
        outstr += "{0:12.6E}\n".format(p.current_gen)
        outstr += 'particle id:\n'
        outstr += "{0:12.6E}\n".format(p.id)
        outstr += 'run mode:\n'
        outstr += "{0:12.6E}\n".format(p.run_mode)
        outstr += 'particle weight:\n'
        outstr += "{0:12.6E}\n".format(p.weight)
        outstr += 'particle energy:\n'
        outstr += "{0:12.6E}\n".format(p.energy)
        outstr += 'particle xyz:\n'
        outstr += "{0:12.6E} {1:12.6E} {2:12.6E}\n".format(p.xyz[0],p.xyz[1],
                                                           p.xyz[2])
        outstr += 'particle uvw:\n'
        outstr += "{0:12.6E} {1:12.6E} {2:12.6E}\n".format(p.uvw[0],p.uvw[1],
                                                           p.uvw[2])

        # Write results to a file.
        with open('results_test.dat','w') as fh:
            fh.write(outstr)


if __name__ == '__main__':
    harness = ParticleRestartTestHarness('particle_7_6144.*')
    harness.execute_test()
