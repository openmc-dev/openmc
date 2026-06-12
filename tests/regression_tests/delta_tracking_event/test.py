import openmc
import pytest
from openmc.utility_funcs import change_directory
from openmc.examples import delta_tracking_lattice

from tests.testing_harness import PyAPITestHarness


@pytest.mark.parametrize("photon", [False, True])
def test_lattice(photon):
    with change_directory(str(photon)):
        model = delta_tracking_lattice(photon)
        model.settings.event_based = True
        harness = PyAPITestHarness('statepoint.10.h5', model)
        harness.main()
