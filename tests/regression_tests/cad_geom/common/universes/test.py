import pytest
from tests.testing_harness import PyAPITestHarness


@pytest.mark.parametrize("interface", ["dagmc", "xdg"])
def test_univ(lattice_model, interface):
    model = lattice_model(interface)
    input_file = f'inputs_{interface}.dat'
    harness = PyAPITestHarness('statepoint.10.h5', model, input_file)
    harness.main()
