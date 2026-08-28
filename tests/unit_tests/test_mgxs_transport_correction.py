"""Tests for user-supplied transport correction ratios in openmc.mgxs.Library."""

import numpy as np
import pytest

import openmc
import openmc.mgxs


@pytest.fixture
def simple_geometry():
    openmc.reset_auto_ids()
    mat = openmc.Material(material_id=1)
    mat.add_nuclide('U235', 1.0)
    mat.set_density('g/cm3', 10.0)
    sph = openmc.Sphere(r=1.0, boundary_type='vacuum')
    cell = openmc.Cell(fill=mat, region=-sph)
    return openmc.Geometry([cell]), mat


@pytest.fixture
def library(simple_geometry):
    geometry, _ = simple_geometry
    groups = openmc.mgxs.EnergyGroups(group_edges=[0.0, 0.625, 2.0e7])
    lib = openmc.mgxs.Library(geometry)
    lib.energy_groups = groups
    lib.domain_type = 'material'
    lib.correction = None
    lib.scatter_format = 'legendre'
    return lib


def test_setter_accepts_and_normalizes(library):
    library.transport_correction_ratios = {'material': {1: [0.9, 0.8]}}
    stored = library.transport_correction_ratios
    assert set(stored) == {'material'}
    assert list(stored['material']) == [1]
    assert isinstance(stored['material'][1], np.ndarray)
    np.testing.assert_allclose(stored['material'][1], [0.9, 0.8])


def test_setter_none_clears(library):
    library.transport_correction_ratios = {'material': {1: [0.9, 0.8]}}
    library.transport_correction_ratios = None
    assert library.transport_correction_ratios is None


def test_setter_validation(library):
    # Not a mapping
    with pytest.raises(TypeError):
        library.transport_correction_ratios = [0.9, 0.8]

    # Invalid domain type
    with pytest.raises(ValueError):
        library.transport_correction_ratios = {'banana': {1: [0.9, 0.8]}}

    # Non-integer domain ID
    with pytest.raises(TypeError):
        library.transport_correction_ratios = {'material': {'1': [0.9, 0.8]}}

    # Wrong number of groups
    with pytest.raises(ValueError):
        library.transport_correction_ratios = {'material': {1: [0.9, 0.8, 0.7]}}

    # Non-positive ratio
    with pytest.raises(ValueError):
        library.transport_correction_ratios = {'material': {1: [0.9, -0.1]}}


def _make_xsdata(groups, sigma_t, scatter, representation='isotropic',
                 num_polar=1, num_azimuthal=1):
    xsdata = openmc.XSdata('set1', groups, representation=representation)
    xsdata.order = 0
    if representation == 'angle':
        xsdata.num_polar = num_polar
        xsdata.num_azimuthal = num_azimuthal
    xsdata.set_total(sigma_t)
    xsdata.set_scatter_matrix(scatter)
    return xsdata


def test_apply_isotropic(library, simple_geometry):
    _, mat = simple_geometry
    library.transport_correction_ratios = {'material': {1: [0.9, 0.8]}}

    sigma_t = np.array([2.0, 3.0])
    scatter = np.array([[[0.5], [0.3]], [[0.1], [1.2]]])
    scatter_orig = scatter.copy()
    absorption = sigma_t - scatter_orig[:, :, 0].sum(axis=1)

    xsdata = _make_xsdata(library.energy_groups, sigma_t, scatter)
    library._apply_transport_correction_ratios(xsdata, mat, 294.0)

    ratios = np.array([0.9, 0.8])
    delta = (1.0 - ratios) * sigma_t

    # Total is transport-corrected
    np.testing.assert_allclose(xsdata._total[0], ratios * sigma_t)

    # In-group P0 diagonal reduced by delta, off-diagonal unchanged
    sm = xsdata._scatter_matrix[0]
    np.testing.assert_allclose(sm[0, 0, 0], scatter_orig[0, 0, 0] - delta[0])
    np.testing.assert_allclose(sm[1, 1, 0], scatter_orig[1, 1, 0] - delta[1])
    np.testing.assert_allclose(sm[0, 1, 0], scatter_orig[0, 1, 0])
    np.testing.assert_allclose(sm[1, 0, 0], scatter_orig[1, 0, 0])

    # Absorption balance (total - out-scatter) is unchanged
    new_absorption = xsdata._total[0] - sm[:, :, 0].sum(axis=1)
    np.testing.assert_allclose(new_absorption, absorption)


def test_apply_angle(simple_geometry):
    geometry, mat = simple_geometry
    groups = openmc.mgxs.EnergyGroups(group_edges=[0.0, 0.625, 2.0e7])
    lib = openmc.mgxs.Library(geometry)
    lib.energy_groups = groups
    lib.domain_type = 'material'
    lib.correction = None
    lib.num_polar = 2
    lib.num_azimuthal = 2
    lib.transport_correction_ratios = {'material': {1: [0.9, 0.8]}}

    sigma_t = np.array([2.0, 3.0])
    scatter2d = np.array([[[0.5], [0.3]], [[0.1], [1.2]]])
    total = np.empty((2, 2, 2))
    total[...] = sigma_t
    scatter = np.zeros((2, 2, 2, 2, 1))
    scatter[...] = scatter2d
    scatter_orig = scatter2d.copy()

    xsdata = _make_xsdata(groups, total, scatter, representation='angle',
                          num_polar=2, num_azimuthal=2)
    lib._apply_transport_correction_ratios(xsdata, mat, 294.0)

    ratios = np.array([0.9, 0.8])
    delta = (1.0 - ratios) * sigma_t
    sm = xsdata._scatter_matrix[0]
    np.testing.assert_allclose(xsdata._total[0], np.broadcast_to(
        ratios * sigma_t, (2, 2, 2)))
    np.testing.assert_allclose(sm[0, 0, 0, 0, 0], scatter_orig[0, 0, 0] - delta[0])
    np.testing.assert_allclose(sm[1, 1, 1, 1, 0], scatter_orig[1, 1, 0] - delta[1])
    np.testing.assert_allclose(sm[:, :, 0, 1, 0], scatter_orig[0, 1, 0])


def test_apply_noop_without_entry(library, simple_geometry):
    _, mat = simple_geometry
    library.transport_correction_ratios = {'material': {999: [0.9, 0.8]}}

    sigma_t = np.array([2.0, 3.0])
    scatter = np.array([[[0.5], [0.3]], [[0.1], [1.2]]])
    xsdata = _make_xsdata(library.energy_groups, sigma_t, scatter)
    library._apply_transport_correction_ratios(xsdata, mat, 294.0)

    # Domain 1 has no entry, so nothing is changed
    np.testing.assert_allclose(xsdata._total[0], sigma_t)


def test_apply_requires_correction_none(library, simple_geometry):
    _, mat = simple_geometry
    library.correction = 'P0'
    library.transport_correction_ratios = {'material': {1: [0.9, 0.8]}}

    sigma_t = np.array([2.0, 3.0])
    scatter = np.array([[[0.5], [0.3]], [[0.1], [1.2]]])
    xsdata = _make_xsdata(library.energy_groups, sigma_t, scatter)
    with pytest.raises(ValueError, match='correction'):
        library._apply_transport_correction_ratios(xsdata, mat, 294.0)


def test_check_library_rejects_p0(library):
    library.mgxs_types = ['total', 'absorption', 'nu-scatter matrix',
                          'scatter matrix']
    library.correction = 'P0'
    library.transport_correction_ratios = {'material': {1: [0.9, 0.8]}}
    with pytest.raises(ValueError, match='Invalid MGXS configuration'):
        library.check_library_for_openmc_mgxs()


def test_check_library_warns_missing_domain_type(library):
    library.mgxs_types = ['total', 'absorption', 'nu-scatter matrix',
                          'scatter matrix']
    # Ratios provided for a domain type the library does not use
    library.transport_correction_ratios = {'cell': {1: [0.9, 0.8]}}
    with pytest.warns(UserWarning, match='do not contain any entries'):
        library.check_library_for_openmc_mgxs()
