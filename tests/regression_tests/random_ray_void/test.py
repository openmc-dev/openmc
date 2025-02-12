import os

import openmc
from openmc.utility_funcs import change_directory
from openmc.examples import random_ray_three_region_cube
import pytest

from tests.testing_harness import TolerantPyAPITestHarness


class MGXSTestHarness(TolerantPyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = 'mgxs.h5'
        if os.path.exists(f):
            os.remove(f)


@pytest.mark.parametrize("shape", ["flat", "linear"])
def test_random_ray_void(shape):
    with change_directory(shape):
        openmc.reset_auto_ids()
        model = random_ray_three_region_cube()

        # There is some different logic for void depending on linear
        # vs. flat, so we test both
        model.settings.random_ray['source_shape'] = shape

        # As we are testing linear sources, need to have more than
        # 10 inactive batches so the moments start getting computed
        model.settings.inactive = 20
        model.settings.batches = 40

        # Begin by getting handles to the cells, and setting the
        # source and void areas to have no fill. We leave the absorber
        # as solid.
        absorber_cell = model.geometry.get_cells_by_name(
            'infinite absorber region', matching=True)[0]
        void_cell = model.geometry.get_cells_by_name(
            'infinite void region', matching=True)[0]
        source_cell = model.geometry.get_cells_by_name(
            'infinite source region', matching=True)[0]

        void_cell.fill = None
        source_cell.fill = None

        # We also need to redefine all three tallies to use cell
        # filters instead of material ones
        estimator = 'tracklength'
        absorber_filter = openmc.CellFilter(absorber_cell)
        absorber_tally = openmc.Tally(name="Absorber Tally")
        absorber_tally.filters = [absorber_filter]
        absorber_tally.scores = ['flux']
        absorber_tally.estimator = estimator

        void_filter = openmc.CellFilter(void_cell)
        void_tally = openmc.Tally(name="Void Tally")
        void_tally.filters = [void_filter]
        void_tally.scores = ['flux']
        void_tally.estimator = estimator

        source_filter = openmc.CellFilter(source_cell)
        source_tally = openmc.Tally(name="Source Tally")
        source_tally.filters = [source_filter]
        source_tally.scores = ['flux']
        source_tally.estimator = estimator

        tallies = openmc.Tallies([source_tally, void_tally, absorber_tally])
        model.tallies = tallies

        harness = MGXSTestHarness('statepoint.40.h5', model)
        harness.main()
