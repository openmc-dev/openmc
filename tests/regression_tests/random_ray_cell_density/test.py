import os

from openmc.examples import random_ray_lattice

from tests.testing_harness import TolerantPyAPITestHarness


class MGXSTestHarness(TolerantPyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = 'mgxs.h5'
        if os.path.exists(f):
            os.remove(f)


def test_random_ray_basic():
    model = random_ray_lattice()
    # Double the densities of all the fuel cells.
    for id, cell in model.geometry.get_all_cells().items():
        if cell.fill.name == "UO2 fuel":
            cell.density = 2.0

    # Gold file was generated with manually scaled fuel cross sections.
    harness = MGXSTestHarness('statepoint.10.h5', model)
    harness.main()
