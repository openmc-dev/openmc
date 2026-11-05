from pathlib import Path

import openmc
import pytest

from tests.regression_tests.temperature_field.harness import TemperatureFieldTestHarness
from tests.regression_tests.temperature_field.models import (
    csg_single_cube_reflective,
    csg_nested_cubes_periodic,
    lattice_nested_cubes_vacuum,
)
from tests.testing_harness import PyAPITestHarness


@pytest.mark.parametrize('builder', [
    csg_single_cube_reflective, csg_nested_cubes_periodic,
    lattice_nested_cubes_vacuum,
])
def test_reference_temperature_assignment(builder):
    model = builder(reference=True)
    assert model.settings.temperature_field is None
    for z in range(2):
        for y in range(2):
            for x in range(2):
                point = (1.25 + 2.5*x, 1.25 + 2.5*y, 1.25 + 2.5*z)
                cell = model.geometry.find(point)[-1]
                assert cell.temperature == 294.0 + 100.0*(x + 2*y + 4*z)
    if builder is not csg_single_cube_reflective:
        assert model.geometry.find((-2.0, -2.0, -2.0))[-1].temperature == 250.0


@pytest.mark.parametrize('matches', [False, True])
def test_regold_requires_equivalence(run_in_tmpdir, monkeypatch, matches):
    """A rejected comparison must preserve both preexisting baseline files."""
    model = csg_single_cube_reflective()
    harness = TemperatureFieldTestHarness(model, model)
    Path('inputs_true.dat').write_text('original inputs')
    Path('results_true.dat').write_text('original results')
    monkeypatch.setattr(PyAPITestHarness, '_run_openmc', lambda self: None)
    monkeypatch.setattr(PyAPITestHarness, '_get_results', lambda self: 'reference')
    monkeypatch.setattr(harness, '_get_results',
                        lambda: 'reference' if matches else 'wrong temperature')
    monkeypatch.setattr(harness, '_test_output_created', lambda: None)
    if matches:
        harness.update_results()
        assert Path('inputs_true.dat').read_text() != 'original inputs'
        assert Path('results_true.dat').read_text() == 'reference'
    else:
        with pytest.raises(AssertionError, match='explicit cell temperatures'):
            harness.update_results()
        assert Path('inputs_true.dat').read_text() == 'original inputs'
        assert Path('results_true.dat').read_text() == 'original results'
