#!/usr/bin/env python

import os
import sys
import glob
import hashlib
sys.path.insert(0, os.pardir)
from testing_harness import PyAPITestHarness
import openmc


class TallyArithmeticTestHarness(PyAPITestHarness):
    def __init__(self, *args, **kwargs):
        super(TallyArithmeticTestHarness, self).__init__(*args, **kwargs)

        # Initialize Mesh
        mesh = openmc.Mesh(mesh_id=1)
        mesh.type = 'regular'
        mesh.dimension = [2, 2, 2]
        mesh.lower_left = [-160.0, -160.0, -183.0]
        mesh.upper_right = [160.0, 160.0, 183.0]

        # Initialize the filters
        energy_filter = openmc.EnergyFilter((0.0, 0.253e-6, 1.0e-3, 1.0, 20.0))
        material_filter  = openmc.MaterialFilter((1, 3))
        distrib_filter   = openmc.DistribcellFilter(60)
        mesh_filter      = openmc.MeshFilter(mesh)

        # Initialized the tallies
        tally = openmc.Tally(name='tally 1')
        tally.filters = [material_filter, energy_filter, distrib_filter]
        tally.scores = ['nu-fission', 'total']
        tally.nuclides = ['U234', 'U235']
        self._model.tallies.append(tally)

        tally = openmc.Tally(name='tally 2')
        tally.filters = [energy_filter, mesh_filter]
        tally.scores = ['total', 'fission']
        tally.nuclides = ['U238', 'U235']
        self._model.tallies.append(tally)

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


if __name__ == '__main__':
    harness = TallyArithmeticTestHarness('statepoint.10.h5')
    harness.main()
