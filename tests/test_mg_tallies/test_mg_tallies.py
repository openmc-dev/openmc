#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import TestHarness, PyAPITestHarness
import openmc


class MGTalliesTestHarness(PyAPITestHarness):
    def _build_inputs(self):
        # Instantiate a tally mesh
        mesh = openmc.Mesh(mesh_id=1)
        mesh.type = 'regular'
        mesh.dimension = [17, 17, 1]
        mesh.lower_left = [0.0, 0.0, 0.0]
        mesh.upper_right = [21.42, 21.42, 100.0]

        # Instantiate some tally filters
        energy_filter = openmc.Filter(type='energy',
                                      bins=[0.0, 20.0])
        energyout_filter = openmc.Filter(type='energyout',
                                         bins=[0.0, 20.0])
        mesh_filter = openmc.Filter()
        mesh_filter.mesh = mesh

        mat_filter = openmc.Filter(type='material', bins=[1,2,3])

        tally1 = openmc.Tally(tally_id=1)
        tally1.filters = [mesh_filter]
        tally1.scores = ['total', 'absorption', 'flux',
                         'fission', 'nu-fission']

        tally2 = openmc.Tally(tally_id=2)
        tally2.filters = [mat_filter, energy_filter, energyout_filter]
        tally2.scores = ['scatter', 'nu-scatter']

        self._input_set.tallies = openmc.Tallies([tally1, tally2])

        super(MGTalliesTestHarness, self)._build_inputs()

    def _cleanup(self):
        super(MGTalliesTestHarness, self)._cleanup()
        f = os.path.join(os.getcwd(), 'tallies.xml')
        if os.path.exists(f): os.remove(f)


if __name__ == '__main__':
    harness = MGTalliesTestHarness('statepoint.10.*', True, mg=True)
    harness.main()
