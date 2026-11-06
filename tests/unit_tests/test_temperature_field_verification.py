from copy import deepcopy

import openmc
import numpy as np
import pytest

from tests.regression_tests import config
from tests.temperature_field import (
    csg_single_cube_reflective,
    csg_nested_cubes_periodic,
    lattice_nested_cubes_vacuum,
)


def _cell_temperature_model(model):
    """Express the test's eight temperatures as cells, without a field."""
    reference = deepcopy(model)
    settings = reference.settings.to_xml_element()
    settings.remove(settings.find('temperature_field'))
    settings.remove(settings.find('mesh'))
    reference.settings = openmc.Settings.from_xml_element(settings)
    # Indexed by z, y, x, with each axis increasing in physical coordinates.
    temperatures = np.array([
        [[294.0, 394.0], [494.0, 594.0]],
        [[694.0, 794.0], [894.0, 994.0]],
    ])
    lattices = reference.geometry.get_all_lattices()
    if lattices:
        lattice = next(iter(lattices.values()))
        universes = np.empty((2, 2, 2), dtype=object)
        for z, y, x in np.ndindex(temperatures.shape):
            # RectLattice stores y rows in descending order.
            universe = lattice.universes[z, 1-y, x].clone(
                clone_materials=False, clone_regions=False)
            for cell in universe.cells.values():
                cell.temperature = temperatures[z, y, x]
            universes[z, 1-y, x] = universe
        lattice.universes = universes
    else:
        root = reference.geometry.root_universe
        field_cell = reference.geometry.find((1.25, 1.25, 1.25))[-1]
        root.remove_cell(field_cell)
        midplanes = [openmc.XPlane(2.5), openmc.YPlane(2.5), openmc.ZPlane(2.5)]
        for z, y, x in np.ndindex(temperatures.shape):
            region = openmc.Intersection([field_cell.region] + [
                +plane if upper else -plane
                for plane, upper in zip(midplanes, (x, y, z))])
            cell = openmc.Cell(fill=field_cell.fill, region=region)
            cell.temperature = temperatures[z, y, x]
            root.add_cell(cell)
    return reference


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
        assert np.isfinite(field_value).all()
        assert np.isfinite(reference_value).all()
        # Both configurations stop at the same interfaces. Permit only
        # accumulated floating-point rounding, not statistical discrepancies.
        np.testing.assert_allclose(
            field_value, reference_value, rtol=1.e-12, atol=0,
            err_msg='Field/cell temperature mismatch')


@pytest.mark.parametrize('builder', [
    csg_single_cube_reflective, csg_nested_cubes_periodic,
    lattice_nested_cubes_vacuum,
])
def test_transport_equivalence(builder, tmp_path):
    """Compare field transport directly with independent cell temperatures."""
    model = builder()
    reference = _cell_temperature_model(model)
    _compare_models(model, reference, tmp_path)


def test_transport_comparison_rejects_wrong_temperatures(tmp_path):
    model = csg_single_cube_reflective()
    reference = _cell_temperature_model(model)
    model.settings.temperature_field.values = [294.0] * 8
    with pytest.raises(AssertionError, match='Field/cell temperature mismatch'):
        _compare_models(model, reference, tmp_path)


@pytest.mark.parametrize('builder', [
    csg_single_cube_reflective, csg_nested_cubes_periodic,
    lattice_nested_cubes_vacuum,
])
def test_reference_temperature_assignment(builder):
    model = _cell_temperature_model(builder())
    assert model.settings.temperature_field is None
    for z in range(2):
        for y in range(2):
            for x in range(2):
                point = (1.25 + 2.5*x, 1.25 + 2.5*y, 1.25 + 2.5*z)
                cell = model.geometry.find(point)[-1]
                assert cell.temperature == 294.0 + 100.0*(x + 2*y + 4*z)
    if builder is not csg_single_cube_reflective:
        assert model.geometry.find((-2.0, -2.0, -2.0))[-1].temperature == 250.0
