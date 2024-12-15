import openmc
import openmc.lib
from openmc.stats import Box

import pytest
from tests.testing_harness import PyAPITestHarness

pytestmark = pytest.mark.skipif(
    not openmc.lib._uwuw_enabled(), reason="UWUW is not enabled."
)


class UWUWTest(PyAPITestHarness):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # settings
        self._model.settings.batches = 5
        self._model.settings.inactive = 0
        self._model.settings.particles = 100

        source = openmc.IndependentSource(space=Box([-4, -4, -4], [4, 4, 4]))
        self._model.settings.source = source

        # geometry
        dag_univ = openmc.DAGMCUniverse("dagmc.h5m", auto_geom_ids=True)
        self._model.geometry = openmc.Geometry(dag_univ)

        # tally
        tally = openmc.Tally()
        tally.scores = ["total"]
        tally.filters = [openmc.CellFilter(2)]
        self._model.tallies = [tally]


def test_refl():
    harness = UWUWTest("statepoint.5.h5", model=openmc.Model())
    harness.main()
