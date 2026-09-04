import pytest

from tests.regression_tests.temperature_field.models import lattice_nested_cubes_vacuum
from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def temperature_field_model():
    return lattice_nested_cubes_vacuum()


def test_temperature_field(temperature_field_model):
    harness = PyAPITestHarness('statepoint.20.h5', temperature_field_model)
    harness.main()
