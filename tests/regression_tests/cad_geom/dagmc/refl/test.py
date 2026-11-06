import openmc
import openmc.lib
from openmc.stats import Box

import pytest
from tests.testing_harness import PyAPITestHarness

class UWUWTest(PyAPITestHarness):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # settings
        self._model.settings.batches = 5
        self._model.settings.inactive = 0
        self._model.settings.particles = 100

        source = openmc.IndependentSource(space=Box([-4, -4, -4],
                                                    [ 4,  4,  4]))
        self._model.settings.source = source

        # geometry
        dag_univ = openmc.DAGMCUniverse("dagmc.h5m", auto_geom_ids=True)
        self._model.geometry = openmc.Geometry(dag_univ)

        # materials
        u235 = openmc.Material(name="fuel")
        u235.add_nuclide('U235', 1.0, 'ao')
        u235.set_density('g/cc', 11)
        u235.id = 40

        water = openmc.Material(name="moderator")
        water.add_nuclide('H1', 2.0, 'ao')
        water.add_nuclide('O16', 1.0, 'ao')
        water.set_density('g/cc', 1.0)
        water.add_s_alpha_beta('c_H_in_H2O')
        water.id = 41

        self._model.materials = openmc.Materials([u235, water])

        # tally
        tally = openmc.Tally()
        tally.scores = ['total']
        tally.filters = [openmc.CellFilter(2)]
        self._model.tallies = [tally]


def test_refl():
    harness = UWUWTest('statepoint.5.h5', model=openmc.Model())
    harness.main()
