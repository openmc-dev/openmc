#!/usr/bin/env python

import os
import sys
import glob
import hashlib
sys.path.insert(0, os.pardir)
from testing_harness import PyAPITestHarness
import openmc


class TallyArithmeticTestHarness(PyAPITestHarness):
    def _build_inputs(self):

        # The summary.h5 file needs to be created to read in the tallies
        self._input_set.settings.output = {'summary': True}

        # Initialize the tallies file
        tallies_file = openmc.Tallies()

        # Initialize the nuclides
        u235 = openmc.Nuclide('U-235')
        u238 = openmc.Nuclide('U-238')
        pu239 = openmc.Nuclide('Pu-239')

        # Initialize Mesh
        mesh = openmc.Mesh(mesh_id=1)
        mesh.type = 'regular'
        mesh.dimension = [2, 2, 2]
        mesh.lower_left = [-160.0, -160.0, -183.0]
        mesh.upper_right = [160.0, 160.0, 183.0]

        # Initialize the filters
        energy_filter = openmc.Filter(type='energy', bins=(0.0, 0.253e-6,
                                                           1.0e-3, 1.0, 20.0))
        material_filter  = openmc.Filter(type='material', bins=(1, 3))
        distrib_filter   = openmc.Filter(type='distribcell', bins=(60))
        mesh_filter      = openmc.Filter(type='mesh')
        mesh_filter.mesh = mesh

        # Initialized the tallies
        tally = openmc.Tally(name='tally 1')
        tally.filters = [material_filter, energy_filter, distrib_filter]
        tally.scores = ['nu-fission', 'total']
        tally.nuclides = [u235, pu239]
        tallies_file.append(tally)

        tally = openmc.Tally(name='tally 2')
        tally.filters = [energy_filter, mesh_filter]
        tally.scores = ['total', 'fission']
        tally.nuclides = [u238, u235]
        tallies_file.append(tally)

        # Export tallies to file
        self._input_set.tallies = tallies_file
        super(TallyArithmeticTestHarness, self)._build_inputs()

    def _get_results(self, hash_output=False):
        """Digest info in the statepoint and return as a string."""

        # Read the statepoint file.
        statepoint = glob.glob(os.path.join(os.getcwd(), self._sp_name))[0]
        sp = openmc.StatePoint(statepoint)

        # Load the tallies
        tally_1 = sp.get_tally(name='tally 1')
        tally_2 = sp.get_tally(name='tally 2')

        # Perform all the tally arithmetic operations and output results
        outstr = ''
        tally_3 = tally_1 * tally_2
        outstr += str(tally_3.mean)

        tally_3 = tally_1.hybrid_product(tally_2, '*', 'entrywise', 'tensor',
                                         'tensor')
        outstr += str(tally_3.mean)

        tally_3 = tally_1.hybrid_product(tally_2, '*', 'entrywise', 'entrywise',
                                         'tensor')
        outstr += str(tally_3.mean)

        tally_3 = tally_1.hybrid_product(tally_2, '*', 'entrywise', 'tensor',
                                         'entrywise')
        outstr += str(tally_3.mean)

        tally_3 = tally_1.hybrid_product(tally_2, '*', 'entrywise', 'entrywise',
                                         'entrywise')
        outstr += str(tally_3.mean)

        # Hash the results if necessary
        if hash_output:
            sha512 = hashlib.sha512()
            sha512.update(outstr.encode('utf-8'))
            outstr = sha512.hexdigest()

        return outstr

    def _cleanup(self):
        super(TallyArithmeticTestHarness, self)._cleanup()
        f = os.path.join(os.getcwd(), 'tallies.xml')
        if os.path.exists(f): os.remove(f)

if __name__ == '__main__':
    harness = TallyArithmeticTestHarness('statepoint.10.h5', True)
    harness.main()
