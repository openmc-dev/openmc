import openmc
import pytest
from openmc.utility_funcs import change_directory
from openmc.examples import delta_tracking_lattice

from tests.testing_harness import PyAPITestHarness


@pytest.mark.parametrize("photon", [False, True])
def test_lattice_checkerboard(photon):
    with change_directory(str(photon)):
        model = delta_tracking_lattice(photon, densities=[10.0, 20.0, 10.0, 20.0])
        harness = PyAPITestHarness('statepoint.10.h5', model)
        harness.main()
