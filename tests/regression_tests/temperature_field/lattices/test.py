import pytest

from tests.regression_tests.temperature_field.models import lattice_nested_cubes_vacuum
from tests.regression_tests.temperature_field.harness import TemperatureFieldTestHarness


@pytest.fixture
def temperature_field_model():
    return lattice_nested_cubes_vacuum()


def test_temperature_field(temperature_field_model):
    harness = TemperatureFieldTestHarness(
        temperature_field_model, lattice_nested_cubes_vacuum(reference=True))
    harness.main()
