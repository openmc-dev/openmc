#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import PyAPITestHarness
from openmc import Filter, Mesh, Tally, TalliesFile
from openmc.source import Source
from openmc.stats import Box

class TalliesTestHarness(PyAPITestHarness):
    def _build_inputs(self):
        # Build default materials/geometry
        self._input_set.build_default_materials_and_geometry()

        # Set settings explicitly
        self._input_set.settings.batches = 5
        self._input_set.settings.inactive = 0
        self._input_set.settings.particles = 400
        self._input_set.settings.source = Source(space=Box(
            [-160, -160, -183], [160, 160, 183]))

        azimuthal_bins = (-3.1416, -1.8850, -0.6283, 0.6283, 1.8850, 3.1416)
        azimuthal_filter1 = Filter(type='azimuthal', bins=azimuthal_bins)
        azimuthal_tally1 = Tally()
        azimuthal_tally1.add_filter(azimuthal_filter1)
        azimuthal_tally1.add_score('flux')
        azimuthal_tally1.estimator = 'tracklength'

        azimuthal_tally2 = Tally()
        azimuthal_tally2.add_filter(azimuthal_filter1)
        azimuthal_tally2.add_score('flux')
        azimuthal_tally2.estimator = 'analog'

        azimuthal_filter2 = Filter(type='azimuthal', bins=(5,))
        azimuthal_tally3 = Tally()
        azimuthal_tally3.add_filter(azimuthal_filter2)
        azimuthal_tally3.add_score('flux')
        azimuthal_tally3.estimator = 'tracklength'

        mesh_2x2 = Mesh(mesh_id=1)
        mesh_2x2.lower_left  = [-182.07, -182.07]
        mesh_2x2.upper_right = [182.07,  182.07]
        mesh_2x2.dimension = [2, 2]
        mesh_filter = Filter(type='mesh', bins=(1,))
        azimuthal_tally4 = Tally()
        azimuthal_tally4.add_filter(azimuthal_filter2)
        azimuthal_tally4.add_filter(mesh_filter)
        azimuthal_tally4.add_score('flux')
        azimuthal_tally4.estimator = 'tracklength'

        cellborn_tally = Tally()
        cellborn_tally.add_filter(Filter(type='cellborn', bins=(10, 21, 22, 23)))
        cellborn_tally.add_score('total')

        dg_tally = Tally()
        dg_tally.add_filter(Filter(type='delayedgroup', bins=(1, 2, 3, 4, 5, 6)))
        dg_tally.add_score('delayed-nu-fission')

        four_groups = (0.0, 0.253e-6, 1.0e-3, 1.0, 20.0)
        energy_filter = Filter(type='energy', bins=four_groups)
        energy_tally = Tally()
        energy_tally.add_filter(energy_filter)
        energy_tally.add_score('total')

        energyout_filter = Filter(type='energyout', bins=four_groups)
        energyout_tally = Tally()
        energyout_tally.add_filter(energyout_filter)
        energyout_tally.add_score('scatter')

        transfer_tally = Tally()
        transfer_tally.add_filter(energy_filter)
        transfer_tally.add_filter(energyout_filter)
        transfer_tally.add_score('scatter')
        transfer_tally.add_score('nu-fission')

        material_tally = Tally()
        material_tally.add_filter(Filter(type='material', bins=(1, 2, 3, 4)))
        material_tally.add_score('total')

        mu_tally1 = Tally()
        mu_tally1.add_filter(Filter(type='mu', bins=(-1.0, -0.5, 0.0, 0.5, 1.0)))
        mu_tally1.add_score('scatter')
        mu_tally1.add_score('nu-scatter')

        mu_filter = Filter(type='mu', bins=(5,))
        mu_tally2 = Tally()
        mu_tally2.add_filter(mu_filter)
        mu_tally2.add_score('scatter')
        mu_tally2.add_score('nu-scatter')

        mu_tally3 = Tally()
        mu_tally3.add_filter(mu_filter)
        mu_tally3.add_filter(mesh_filter)
        mu_tally3.add_score('scatter')
        mu_tally3.add_score('nu-scatter')

        polar_bins = (0.0, 0.6283, 1.2566, 1.8850, 2.5132, 3.1416)
        polar_filter = Filter(type='polar', bins=polar_bins)
        polar_tally1 = Tally()
        polar_tally1.add_filter(polar_filter)
        polar_tally1.add_score('flux')
        polar_tally1.estimator = 'tracklength'

        polar_tally2 = Tally()
        polar_tally2.add_filter(polar_filter)
        polar_tally2.add_score('flux')
        polar_tally2.estimator = 'analog'

        polar_filter2 = Filter(type='polar', bins=(5,))
        polar_tally3 = Tally()
        polar_tally3.add_filter(polar_filter2)
        polar_tally3.add_score('flux')
        polar_tally3.estimator = 'tracklength'

        polar_tally4 = Tally()
        polar_tally4.add_filter(polar_filter2)
        polar_tally4.add_filter(mesh_filter)
        polar_tally4.add_score('flux')
        polar_tally4.estimator = 'tracklength'

        universe_tally = Tally()
        universe_tally.add_filter(Filter(type='universe', bins=(1, 2, 3, 4)))
        universe_tally.add_score('total')

        cell_filter = Filter(type='cell', bins=(10, 21, 22, 23))
        score_tallies = [Tally(), Tally(), Tally()]
        for t in score_tallies:
            t.add_filter(cell_filter)
            t.add_score('absorption')
            t.add_score('delayed-nu-fission')
            t.add_score('events')
            t.add_score('fission')
            t.add_score('inverse-velocity')
            t.add_score('kappa-fission')
            t.add_score('(n,2n)')
            t.add_score('(n,n1)')
            t.add_score('(n,gamma)')
            t.add_score('nu-fission')
            t.add_score('scatter')
            t.add_score('elastic')
            t.add_score('total')
        score_tallies[0].estimator = 'tracklength'
        score_tallies[1].estimator = 'analog'
        score_tallies[2].estimator = 'collision'

        cell_filter2 = Filter(type='cell', bins=(21, 22, 23, 27, 28, 29))
        flux_tallies = [Tally() for i in range(4)]
        [t.add_filter(cell_filter2) for t in flux_tallies]
        flux_tallies[0].add_score('flux')
        [t.add_score('flux-y5') for t in flux_tallies[1:]]
        flux_tallies[1].estimator = 'tracklength'
        flux_tallies[2].estimator = 'analog'
        flux_tallies[3].estimator = 'collision'

        scatter_tally1 = Tally()
        scatter_tally1.add_filter(cell_filter)
        scatter_tally1.add_score('scatter')
        scatter_tally1.add_score('scatter-1')
        scatter_tally1.add_score('scatter-2')
        scatter_tally1.add_score('scatter-3')
        scatter_tally1.add_score('scatter-4')
        scatter_tally1.add_score('nu-scatter')
        scatter_tally1.add_score('nu-scatter-1')
        scatter_tally1.add_score('nu-scatter-2')
        scatter_tally1.add_score('nu-scatter-3')
        scatter_tally1.add_score('nu-scatter-4')

        scatter_tally2 = Tally()
        scatter_tally2.add_filter(cell_filter)
        scatter_tally2.add_score('scatter-p4')
        scatter_tally2.add_score('scatter-y4')
        scatter_tally2.add_score('nu-scatter-p4')
        scatter_tally2.add_score('nu-scatter-y3')

        total_tallies = [Tally() for i in range(4)]
        [t.add_filter(cell_filter) for t in total_tallies]
        total_tallies[0].add_score('total')
        [t.add_score('total-y4') for t in total_tallies[1:]]
        [t.add_nuclide('U-235') for t in total_tallies[1:]]
        [t.add_nuclide('total') for t in total_tallies[1:]]
        total_tallies[1].estimator = 'tracklength'
        total_tallies[2].estimator = 'analog'
        total_tallies[3].estimator = 'collision'

        questionable_tally = Tally()
        questionable_tally.add_score('transport')
        questionable_tally.add_score('n1n')

        all_nuclide_tallies = [Tally(), Tally()]
        for t in all_nuclide_tallies:
            t.add_filter(cell_filter)
            t.add_nuclide('all')
            t.add_score('total')
        all_nuclide_tallies[0].estimator = 'tracklength'
        all_nuclide_tallies[0].estimator = 'collision'

        self._input_set.tallies = TalliesFile()
        self._input_set.tallies.add_tally(azimuthal_tally1)
        self._input_set.tallies.add_tally(azimuthal_tally2)
        self._input_set.tallies.add_tally(azimuthal_tally3)
        self._input_set.tallies.add_tally(azimuthal_tally4)
        self._input_set.tallies.add_tally(cellborn_tally)
        self._input_set.tallies.add_tally(dg_tally)
        self._input_set.tallies.add_tally(energy_tally)
        self._input_set.tallies.add_tally(energyout_tally)
        self._input_set.tallies.add_tally(transfer_tally)
        self._input_set.tallies.add_tally(material_tally)
        self._input_set.tallies.add_tally(mu_tally1)
        self._input_set.tallies.add_tally(mu_tally2)
        self._input_set.tallies.add_tally(mu_tally3)
        self._input_set.tallies.add_tally(polar_tally1)
        self._input_set.tallies.add_tally(polar_tally2)
        self._input_set.tallies.add_tally(polar_tally3)
        self._input_set.tallies.add_tally(polar_tally4)
        self._input_set.tallies.add_tally(universe_tally)
        [self._input_set.tallies.add_tally(t) for t in score_tallies]
        [self._input_set.tallies.add_tally(t) for t in flux_tallies]
        self._input_set.tallies.add_tally(scatter_tally1)
        self._input_set.tallies.add_tally(scatter_tally2)
        [self._input_set.tallies.add_tally(t) for t in total_tallies]
        self._input_set.tallies.add_tally(questionable_tally)
        [self._input_set.tallies.add_tally(t) for t in all_nuclide_tallies]
        self._input_set.tallies.add_mesh(mesh_2x2)

        self._input_set.export()

    def _get_results(self):
        return super(TalliesTestHarness, self)._get_results(hash_output=True)

    def _cleanup(self):
        super(TalliesTestHarness, self)._cleanup()
        f = os.path.join(os.getcwd(), 'tallies.xml')
        if os.path.exists(f): os.remove(f)


if __name__ == '__main__':
    harness = TalliesTestHarness('statepoint.5.*', True)
    harness.main()
