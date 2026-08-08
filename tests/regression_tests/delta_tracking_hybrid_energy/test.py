import openmc
import pytest
from openmc.utility_funcs import change_directory
from openmc.examples import delta_tracking_lattice

from tests.testing_harness import PyAPITestHarness


@pytest.mark.parametrize("photon", [False, True])
@pytest.mark.parametrize("event", [False, True])
def test_lattice(photon, event):
    delta_tracking_settings = {
        'enable' : True,
        'hybrid_type' : 'energy',
        'neutron_energy_threshold' : 5e4,
        'photon_energy_threshold' : 1e6
    }
    with change_directory("photon_"+str(photon)+"_event_"+str(event)):
        model = delta_tracking_lattice(delta_tracking_settings, run_photon=photon)
        model.settings.event_based = event
        harness = PyAPITestHarness('statepoint.10.h5', model)
        harness.main()
