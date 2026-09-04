"""Tests for transport correction ratios in openmc.mgxs.Library.

These cover both user-supplied ratios and the ratios that
:meth:`openmc.mgxs.Library.load_from_statepoint` computes automatically from
tallied transport data.
"""

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


def test_apply_uses_stored_ratio_with_p0(library, simple_geometry, monkeypatch):
    # With correction='P0' the dataset arrives with the tally-based transport
    # total, and the plain total is recovered from the transport MGXS. The
    # stored ratio (as if edited by the user) then determines the correction.
    _, mat = simple_geometry
    library.correction = 'P0'
    library.mgxs_types = ['transport', 'nu-scatter matrix', 'scatter matrix']
    library.transport_correction_ratios = {'material': {1: [0.9, 0.8]}}

    # sigma_t recovered from the transport MGXS
    sigma_t = np.array([2.0, 3.0])
    monkeypatch.setattr(library, 'get_mgxs',
                        lambda domain, mgxs_type: object())
    monkeypatch.setattr(
        library, '_uncorrected_total_xs',
        lambda tm, nuclides, xs_type, subdomains: sigma_t)

    # The tally-based transport total in the dataset differs from
    # ratios * sigma_t, so the stored ratio must visibly take effect.
    corrected_total = np.array([1.5, 2.0])
    scatter = np.array([[[0.5], [0.3]], [[0.1], [1.2]]])
    scatter_orig = scatter.copy()
    xsdata = _make_xsdata(library.energy_groups, corrected_total, scatter)

    library._apply_transport_correction_ratios(xsdata, mat, 294.0)

    ratios = np.array([0.9, 0.8])
    target_total = ratios * sigma_t
    delta = target_total - corrected_total

    np.testing.assert_allclose(xsdata._total[0], target_total)
    sm = xsdata._scatter_matrix[0]
    np.testing.assert_allclose(sm[0, 0, 0], scatter_orig[0, 0, 0] + delta[0])
    np.testing.assert_allclose(sm[1, 1, 0], scatter_orig[1, 1, 0] + delta[1])
    np.testing.assert_allclose(sm[0, 1, 0], scatter_orig[0, 1, 0])
    np.testing.assert_allclose(sm[1, 0, 0], scatter_orig[1, 0, 0])


def test_apply_p0_ratio_matches_tally_is_noop(library, simple_geometry,
                                              monkeypatch):
    # When the stored ratio equals the tally-derived sigma_tr / sigma_t (as it
    # is after automatic population), re-applying it reproduces the tally-based
    # total and leaves the scattering matrix unchanged (no double correction).
    _, mat = simple_geometry
    library.correction = 'P0'
    library.mgxs_types = ['transport', 'nu-scatter matrix', 'scatter matrix']

    sigma_t = np.array([2.0, 3.0])
    corrected_total = np.array([1.8, 2.4])  # sigma_tr from the tallies
    ratios = corrected_total / sigma_t
    library.transport_correction_ratios = {'material': {1: list(ratios)}}

    monkeypatch.setattr(library, 'get_mgxs',
                        lambda domain, mgxs_type: object())
    monkeypatch.setattr(
        library, '_uncorrected_total_xs',
        lambda tm, nuclides, xs_type, subdomains: sigma_t)

    scatter = np.array([[[0.5], [0.3]], [[0.1], [1.2]]])
    scatter_orig = scatter.copy()
    xsdata = _make_xsdata(library.energy_groups, corrected_total, scatter)
    library._apply_transport_correction_ratios(xsdata, mat, 294.0)

    np.testing.assert_allclose(xsdata._total[0], corrected_total)
    np.testing.assert_allclose(xsdata._scatter_matrix[0][:, :, 0],
                               scatter_orig[:, :, 0])


def test_check_library_allows_ratios_with_p0(library):
    # With correction='P0' the stored ratios are applied by recovering the
    # plain total from the transport MGXS, so a valid P0 configuration with a
    # 'transport' type must pass validation without error.
    library.mgxs_types = ['transport', 'absorption', 'nu-scatter matrix',
                          'scatter matrix']
    library.correction = 'P0'
    library.transport_correction_ratios = {'material': {1: [0.9, 0.8]}}
    library.check_library_for_openmc_mgxs()


def test_check_library_p0_requires_transport(library):
    # With correction='P0' and no 'transport'/'nu-transport' type there is no
    # way to recover the plain total, so validation must fail.
    library.mgxs_types = ['absorption', 'nu-scatter matrix', 'scatter matrix']
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


def test_store_ratios_noop_without_transport(library):
    # Without a 'transport' or 'nu-transport' MGXS type there is nothing to
    # compute, so no ratios are stored (keeping the standard total-based
    # workflow unaffected).
    library.mgxs_types = ['total', 'absorption', 'scatter matrix']
    library._store_computed_transport_correction_ratios()
    assert library.transport_correction_ratios is None


def test_store_ratios_skips_angle(library):
    # Angle-dependent data does not yield a single ratio per group, so the
    # helper returns before touching any tallies.
    library.mgxs_types = ['transport', 'absorption', 'scatter matrix']
    library.num_polar = 2
    library.num_azimuthal = 2
    library._store_computed_transport_correction_ratios()
    assert library.transport_correction_ratios is None


def test_store_ratios_preserves_user_entries(library):
    # A user-provided entry must never be overwritten by the automatic
    # computation, even when a transport MGXS type is present.
    library.mgxs_types = ['transport', 'absorption', 'scatter matrix']
    library.num_polar = 2  # skip the tally-based computation for this test
    library.num_azimuthal = 2
    library.transport_correction_ratios = {'material': {1: [0.9, 0.8]}}
    library._store_computed_transport_correction_ratios()
    np.testing.assert_allclose(
        library.transport_correction_ratios['material'][1], [0.9, 0.8])
