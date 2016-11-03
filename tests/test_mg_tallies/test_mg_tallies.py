#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import HashedPyAPITestHarness
import openmc


class MGTalliesTestHarness(HashedPyAPITestHarness):
    def _build_inputs(self):
        """Write input XML files."""
        self._input_set.build_default_materials_and_geometry(as_macro=False)
        self._input_set.build_default_settings()

        # Instantiate a tally mesh
        mesh = openmc.Mesh(mesh_id=1)
        mesh.type = 'regular'
        mesh.dimension = [1, 1, 10]
        mesh.lower_left = [0.0, 0.0, 0.0]
        mesh.upper_right = [10, 10, 5]

        # Instantiate some tally filters
        energy_filter = openmc.EnergyFilter([0.0, 20.0e6])
        energyout_filter = openmc.EnergyoutFilter([0.0, 20.0e6])
        mesh_filter = openmc.MeshFilter(mesh)

        mat_ids = [mat.id for mat in self._input_set.materials]
        mat_filter = openmc.MaterialFilter(mat_ids)

        nuclides = [xs.name for xs in self._input_set.xs_data]

        scores= {False: ['total', 'absorption', 'flux', 'fission', 'nu-fission'],
                 True: ['total', 'absorption', 'fission', 'nu-fission']}

        tallies = []
        for do_nuclides in [False, True]:
            tallies.append(openmc.Tally())
            tallies[-1].filters = [mesh_filter]
            tallies[-1].estimator = 'analog'
            tallies[-1].scores = scores[do_nuclides]
            if do_nuclides:
                tallies[-1].nuclides = nuclides

            tallies.append(openmc.Tally())
            tallies[-1].filters = [mesh_filter]
            tallies[-1].estimator = 'tracklength'
            tallies[-1].scores = scores[do_nuclides]
            if do_nuclides:
                tallies[-1].nuclides = nuclides

            tallies.append(openmc.Tally())
            tallies[-1].filters = [mat_filter, energy_filter]
            tallies[-1].estimator = 'analog'
            tallies[-1].scores = scores[do_nuclides] + ['scatter',
                                                        'nu-scatter']
            if do_nuclides:
                tallies[-1].nuclides = nuclides

            tallies.append(openmc.Tally())
            tallies[-1].filters = [mat_filter, energy_filter]
            tallies[-1].estimator = 'collision'
            tallies[-1].scores = scores[do_nuclides]
            if do_nuclides:
                tallies[-1].nuclides = nuclides

            tallies.append(openmc.Tally())
            tallies[-1].filters = [mat_filter, energy_filter]
            tallies[-1].estimator = 'tracklength'
            tallies[-1].scores = scores[do_nuclides]
            if do_nuclides:
                tallies[-1].nuclides = nuclides

            tallies.append(openmc.Tally())
            tallies[-1].filters = [mat_filter, energy_filter,
                                   energyout_filter]
            tallies[-1].scores = ['scatter', 'nu-scatter', 'nu-fission']
            if do_nuclides:
                tallies[-1].nuclides = nuclides

        self._input_set.tallies = openmc.Tallies(tallies)

        self._input_set.export()

    def _cleanup(self):
        super(MGTalliesTestHarness, self)._cleanup()
        f = os.path.join(os.getcwd(), 'tallies.xml')
        if os.path.exists(f): os.remove(f)


if __name__ == '__main__':
    harness = MGTalliesTestHarness('statepoint.10.*', True, mg=True)
    harness.main()
