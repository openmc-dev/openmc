import openmc
import pytest
from openmc.utility_funcs import change_directory
from openmc.examples import delta_tracking_lattice

from tests.testing_harness import PyAPITestHarness


@pytest.mark.parametrize("boundary", ['vacuum', 'reflective', 'periodic', 'white'])
def test_lattice(boundary):
    with change_directory(boundary):
        model = delta_tracking_lattice(False, boundary)
        harness = PyAPITestHarness('statepoint.10.h5', model)
        harness.main()
