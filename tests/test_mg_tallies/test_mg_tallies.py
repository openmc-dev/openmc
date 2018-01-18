#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import HashedPyAPITestHarness
import openmc
from openmc.examples import slab_mg


if __name__ == '__main__':
    model = slab_mg(as_macro=False)

    # Instantiate a tally mesh
    mesh = openmc.Mesh(mesh_id=1)
    mesh.type = 'regular'
    mesh.dimension = [1, 1, 10]
    mesh.lower_left = [0.0, 0.0, 0.0]
    mesh.upper_right = [10, 10, 5]

    # Instantiate some tally filters
    energy_filter = openmc.EnergyFilter([0.0, 20.0e6])
    energyout_filter = openmc.EnergyoutFilter([0.0, 20.0e6])
    energies = [1e-5, 0.0635, 10.0, 1.0e2, 1.0e3, 0.5e6, 1.0e6, 20.0e6]
    matching_energy_filter = openmc.EnergyFilter(energies)
    matching_eout_filter = openmc.EnergyoutFilter(energies)
    mesh_filter = openmc.MeshFilter(mesh)

    mat_filter = openmc.MaterialFilter(model.materials)

    nuclides = model.xs_data

    scores = {False: ['total', 'absorption', 'flux', 'fission', 'nu-fission'],
              True: ['total', 'absorption', 'fission', 'nu-fission']}

    for do_nuclides in [False, True]:
        t = openmc.Tally()
        t.filters = [mesh_filter]
        t.estimator = 'analog'
        t.scores = scores[do_nuclides]
        if do_nuclides:
            t.nuclides = nuclides
        model.tallies.append(t)

        t = openmc.Tally()
        t.filters = [mesh_filter]
        t.estimator = 'tracklength'
        t.scores = scores[do_nuclides]
        if do_nuclides:
            t.nuclides = nuclides
        model.tallies.append(t)

        # Impose energy bins that dont match the MG structure and those
        # that do
        for match_energy_bins in [False, True]:
            if match_energy_bins:
                e_filter = matching_energy_filter
                eout_filter = matching_eout_filter
            else:
                e_filter = energy_filter
                eout_filter = energyout_filter

            t = openmc.Tally()
            t.filters = [mat_filter, e_filter]
            t.estimator = 'analog'
            t.scores = scores[do_nuclides] + ['scatter', 'nu-scatter']
            if do_nuclides:
                t.nuclides = nuclides
            model.tallies.append(t)

            t = openmc.Tally()
            t.filters = [mat_filter, e_filter]
            t.estimator = 'collision'
            t.scores = scores[do_nuclides]
            if do_nuclides:
                t.nuclides = nuclides
            model.tallies.append(t)

            t = openmc.Tally()
            t.filters = [mat_filter, e_filter]
            t.estimator = 'tracklength'
            t.scores = scores[do_nuclides]
            if do_nuclides:
                t.nuclides = nuclides
            model.tallies.append(t)

            t = openmc.Tally()
            t.filters = [mat_filter, e_filter, eout_filter]
            t.scores = ['scatter', 'nu-scatter', 'nu-fission']
            if do_nuclides:
                t.nuclides = nuclides
            model.tallies.append(t)

    harness = HashedPyAPITestHarness('statepoint.10.h5', model)
    harness.main()
