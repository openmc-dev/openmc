import openmc
import pytest
from openmc.utility_funcs import change_directory
from openmc.examples import delta_tracking_lattice

from tests.testing_harness import PyAPITestHarness


@pytest.mark.parametrize("boundary", ['vacuum', 'reflective', 'periodic', 'white'])
def test_lattice(boundary):
    delta_tracking_settings = {
        'enable' : True,
        'hybrid_type' : 'cross_section',
        'xs_threshold' : 1.0
    }
    with change_directory(boundary):
        model = delta_tracking_lattice(delta_tracking_settings, run_photon=False, boundary_type=boundary)
        harness = PyAPITestHarness('statepoint.10.h5', model)
        harness.main()
