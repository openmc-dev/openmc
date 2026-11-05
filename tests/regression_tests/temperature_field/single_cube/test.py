import pytest

from tests.regression_tests.temperature_field.models import (
    csg_single_cube_reflective,
)
from tests.regression_tests.temperature_field.harness import TemperatureFieldTestHarness


@pytest.fixture
def temperature_field_model():
    return csg_single_cube_reflective()


def test_temperature_field(temperature_field_model):
    harness = TemperatureFieldTestHarness(
        temperature_field_model, csg_single_cube_reflective(reference=True))
    harness.main()
