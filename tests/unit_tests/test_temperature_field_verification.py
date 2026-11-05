from pathlib import Path

import openmc
import openmc.lib
import numpy as np
import pytest

from tests.regression_tests import config
from tests.temperature_field import (
    csg_single_cube_reflective,
    csg_nested_cubes_periodic,
    lattice_nested_cubes_vacuum,
    dagmc_nested_cubes_reflective,
)


def _run(model, directory):
    kwargs = {'openmc_exec': config['exe'], 'event_based': config['event']}
    if config['mpi']:
        kwargs['mpi_args'] = [config['mpiexec'], '-n', config['mpi_np']]
    statepoint = model.run(cwd=directory, **kwargs)
    with openmc.StatePoint(statepoint) as sp:
        return [np.array([sp.keff.n, sp.keff.s])] + [
            np.array([tally.sum, tally.sum_sq]) for tally in sp.tallies.values()]


def _compare_models(model, reference, directory):
    actual = _run(model, directory / 'field')
    expected = _run(reference, directory / 'reference')
    assert len(actual) == len(expected)
    for field_value, reference_value in zip(actual, expected):
        # Both configurations stop at the same interfaces. Permit only
        # accumulated floating-point rounding, not statistical discrepancies.
        np.testing.assert_allclose(field_value, reference_value, rtol=1.e-12,
                                   atol=0, err_msg='Field/cell temperature mismatch')


@pytest.mark.parametrize('builder', [
    csg_single_cube_reflective, csg_nested_cubes_periodic,
    lattice_nested_cubes_vacuum,
])
def test_transport_equivalence(builder, tmp_path):
    """Compare field transport directly with independent cell temperatures."""
    _compare_models(builder(), builder(reference=True), tmp_path)


@pytest.mark.skipif(not openmc.lib.feature_enabled('dagmc'), reason='DAGMC is disabled')
def test_dagmc_transport_equivalence(tmp_path):
    assets = Path(__file__).with_name('temperature_field')
    field_mesh = (Path(__file__).parents[1] / 'regression_tests' /
                  'temperature_field/dagmc/nested_cubes_cell_reflecting.h5m')
    model = dagmc_nested_cubes_reflective(field_mesh)
    reference = dagmc_nested_cubes_reflective(
        assets / 'nested_cubes_multicell_reflecting.h5m', reference=True)
    _compare_models(model, reference, tmp_path)


def test_transport_comparison_rejects_wrong_temperatures(tmp_path):
    model = csg_single_cube_reflective()
    model.settings.temperature_field.values = [294.0] * 8
    reference = csg_single_cube_reflective(reference=True)
    with pytest.raises(AssertionError, match='Field/cell temperature mismatch'):
        _compare_models(model, reference, tmp_path)


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
