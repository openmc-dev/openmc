from pathlib import Path

import openmc
import pytest

from tests.regression_tests.temperature_field.models import (
    dagmc_nested_cubes_reflective,
)
from tests.regression_tests.temperature_field.harness import TemperatureFieldTestHarness


pytestmark = pytest.mark.skipif(
    not openmc.lib.feature_enabled('dagmc'),
    reason="DAGMC CAD geometry is not enabled.")


@pytest.fixture
def temperature_field_model():
    return dagmc_nested_cubes_reflective("nested_cubes_cell_reflecting.h5m")


def test_temperature_field(temperature_field_model):
    reference = dagmc_nested_cubes_reflective(
        Path('nested_cubes_multicell_reflecting.h5m').resolve(), reference=True)
    harness = TemperatureFieldTestHarness(temperature_field_model, reference)
    harness.main()
