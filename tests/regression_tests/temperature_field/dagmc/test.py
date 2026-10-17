import openmc
import pytest

from tests.regression_tests.temperature_field.models import (
    dagmc_nested_cubes_reflective,
)
from tests.testing_harness import PyAPITestHarness


pytestmark = pytest.mark.skipif(
    not openmc.lib.feature_enabled('dagmc'),
    reason="DAGMC CAD geometry is not enabled.")


@pytest.fixture
def temperature_field_model():
    return dagmc_nested_cubes_reflective("nested_cubes_cell_reflecting.h5m")


def test_temperature_field(temperature_field_model):
    harness = PyAPITestHarness('statepoint.20.h5', temperature_field_model)
    harness.main()
