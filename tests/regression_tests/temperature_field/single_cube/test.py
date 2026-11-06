import pytest

from tests.temperature_field import (
    csg_single_cube_reflective,
)
from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def temperature_field_model():
    return csg_single_cube_reflective()


def test_temperature_field(temperature_field_model):
    harness = PyAPITestHarness('statepoint.20.h5', temperature_field_model)
    harness.main()
