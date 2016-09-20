#!/usr/bin/env python

import os
import sys
import glob
import hashlib
sys.path.insert(0, os.pardir)
from testing_harness import HashedPyAPITestHarness
import openmc


class FilterMeshTestHarness(HashedPyAPITestHarness):
    def _build_inputs(self):

        # The summary.h5 file needs to be created to read in the tallies
        self._input_set.settings.output = {'summary': True}

        # Initialize the tallies file
        tallies_file = openmc.Tallies()

        # Initialize Meshes
        mesh_1d = openmc.Mesh(mesh_id=1)
        mesh_1d.type = 'regular'
        mesh_1d.dimension = [17]
        mesh_1d.lower_left = [-182.07]
        mesh_1d.upper_right = [182.07]

        mesh_2d = openmc.Mesh(mesh_id=2)
        mesh_2d.type = 'regular'
        mesh_2d.dimension = [17, 17]
        mesh_2d.lower_left = [-182.07, -182.07]
        mesh_2d.upper_right = [182.07, 182.07]

        mesh_3d = openmc.Mesh(mesh_id=3)
        mesh_3d.type = 'regular'
        mesh_3d.dimension = [17, 17, 17]
        mesh_3d.lower_left = [-182.07, -182.07, -183.00]
        mesh_3d.upper_right = [182.07, 182.07, 183.00]

        # Initialize the filters
        mesh_1d_filter      = openmc.Filter(type='mesh')
        mesh_2d_filter      = openmc.Filter(type='mesh')
        mesh_3d_filter      = openmc.Filter(type='mesh')
        mesh_1d_filter.mesh = mesh_1d
        mesh_2d_filter.mesh = mesh_2d
        mesh_3d_filter.mesh = mesh_3d

        # Initialized the tallies
        tally = openmc.Tally(name='tally 1')
        tally.filters = [mesh_1d_filter]
        tally.scores = ['total']
        tallies_file.append(tally)

        tally = openmc.Tally(name='tally 2')
        tally.filters = [mesh_1d_filter]
        tally.scores = ['current']
        tallies_file.append(tally)

        tally = openmc.Tally(name='tally 3')
        tally.filters = [mesh_2d_filter]
        tally.scores = ['total']
        tallies_file.append(tally)

        tally = openmc.Tally(name='tally 4')
        tally.filters = [mesh_2d_filter]
        tally.scores = ['current']
        tallies_file.append(tally)

        tally = openmc.Tally(name='tally 5')
        tally.filters = [mesh_3d_filter]
        tally.scores = ['total']
        tallies_file.append(tally)

        tally = openmc.Tally(name='tally 6')
        tally.filters = [mesh_3d_filter]
        tally.scores = ['current']
        tallies_file.append(tally)

        # Export tallies to file
        self._input_set.tallies = tallies_file
        super(FilterMeshTestHarness, self)._build_inputs()

    def _cleanup(self):
        super(FilterMeshTestHarness, self)._cleanup()
        f = os.path.join(os.getcwd(), 'tallies.xml')
        if os.path.exists(f): os.remove(f)


if __name__ == '__main__':
    harness = FilterMeshTestHarness('statepoint.10.*', True)
    harness.main()
