import openmc
import openmc.lib

import pytest
from tests.testing_harness import PyAPITestHarness

pytestmark = pytest.mark.skipif(
    not openmc.lib.feature_enabled('dagmc'),
    reason="DAGMC CAD geometry is not enabled.")


test_libs = []
if openmc.lib.feature_enabled('dagmc'):
    test_libs.append('dagmc')
if openmc.lib.feature_enabled('xdg'):
    test_libs.append('xdg')


@pytest.mark.parametrize("interface", test_libs)
def test_univ(lattice_model, interface):
    model = lattice_model(interface)
    input_file = f'inputs_{interface}.dat'
    harness = PyAPITestHarness('statepoint.10.h5', model, input_file)
    harness.main()
