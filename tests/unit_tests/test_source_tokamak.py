import numpy as np
import pytest

import openmc
import openmc.stats

from tests.unit_tests import assert_sample_mean


def make_source(**kwargs):
    """Build a valid TokamakSource, overriding defaults via kwargs."""
    r_over_a = np.linspace(0.0, 1.0, 10)
    params = dict(
        major_radius=620.0,
        minor_radius=200.0,
        elongation=1.8,
        triangularity=0.45,
        shafranov_shift=10.0,
        r_over_a=r_over_a,
        emission_density=(1.0 - r_over_a**2),
        energy=openmc.stats.muir(e0=14.08e6, m_rat=5.0, kt=2.0e4),
    )
    params.update(kwargs)
    return openmc.TokamakSource(**params)


def test_tokamak_source_roundtrip():
    src = make_source(
        phi_start=0.1, phi_extent=np.pi, n_alpha=51, vertical_shift=5.0,
        strength=2.0, time=openmc.stats.Uniform(0.0, 1e-6))

    elem = src.to_xml_element()
    assert elem.get('type') == 'tokamak'

    new = openmc.SourceBase.from_xml_element(elem)
    assert isinstance(new, openmc.TokamakSource)
    assert new.major_radius == src.major_radius
    assert new.minor_radius == src.minor_radius
    assert new.elongation == src.elongation
    assert new.triangularity == src.triangularity
    assert new.shafranov_shift == src.shafranov_shift
    assert new.phi_start == src.phi_start
    assert new.phi_extent == src.phi_extent
    assert new.n_alpha == src.n_alpha
    assert new.vertical_shift == src.vertical_shift
    assert new.strength == src.strength
    np.testing.assert_allclose(new.r_over_a, src.r_over_a)
    np.testing.assert_allclose(new.emission_density, src.emission_density)
    assert len(new.energy) == 1
    assert isinstance(new.time, openmc.stats.Uniform)
    assert new.time.a == src.time.a
    assert new.time.b == src.time.b


def test_tokamak_source_default_time():
    src = make_source()
    assert src.time is None

    new = openmc.SourceBase.from_xml_element(src.to_xml_element())
    assert new.time is None

    with pytest.raises(TypeError):
        make_source(time=1.0)


def test_tokamak_source_multiple_energies():
    r_over_a = np.linspace(0.0, 1.0, 5)
    energies = [openmc.stats.muir(e0=14.08e6, m_rat=5.0, kt=kt)
                for kt in (1.0e4, 1.5e4, 2.0e4, 2.5e4, 3.0e4)]
    src = make_source(r_over_a=r_over_a,
                      emission_density=np.ones_like(r_over_a),
                      energy=energies)
    assert len(src.energy) == len(r_over_a)

    new = openmc.SourceBase.from_xml_element(src.to_xml_element())
    assert len(new.energy) == len(r_over_a)


@pytest.mark.parametrize("kwargs, match", [
    (dict(minor_radius=700.0), "smaller than major_radius"),
    (dict(shafranov_shift=150.0), "half the minor_radius"),
    (dict(emission_density=np.ones(5)), "same length as r_over_a"),
    (dict(energy=[openmc.stats.muir(14.08e6, 5.0, 2.0e4)] * 2),
     "Number of energy distributions"),
    (dict(r_over_a=np.linspace(0.1, 1.0, 10)), "must start at 0"),
    (dict(r_over_a=np.linspace(0.0, 0.9, 10)), "must end at 1"),
    (dict(emission_density=-np.linspace(0.0, 1.0, 10)), "cannot be negative"),
    (dict(emission_density=np.zeros(10)), "must contain a positive value"),
])
def test_tokamak_source_invalid(kwargs, match):
    with pytest.raises(ValueError, match=match):
        make_source(**kwargs)


@pytest.mark.parametrize("value", [-1.5, 1.5])
def test_tokamak_source_invalid_triangularity(value):
    with pytest.raises(ValueError):
        make_source(triangularity=value)


def test_tokamak_source_invalid_n_alpha():
    with pytest.raises(ValueError):
        make_source(n_alpha=2)


@pytest.mark.parametrize(("attribute", "value", "match"), [
    ("minor_radius", 700.0, "smaller than major_radius"),
    ("shafranov_shift", 150.0, "half the minor_radius"),
    ("emission_density", np.ones(5), "same length as r_over_a"),
    ("energy", [openmc.stats.delta_function(1.0)] * 2,
     "Number of energy distributions"),
])
def test_tokamak_source_mutation_validation(attribute, value, match):
    src = make_source()
    setattr(src, attribute, value)
    with pytest.raises(ValueError, match=match):
        src.to_xml_element()


@pytest.mark.parametrize(("emission_density", "expected_mean"), [
    ([1.0, 1.0], 2.0 / 3.0),
    ([1.0, 0.0], 0.5),
])
def test_tokamak_source_radial_sampling(
    run_in_tmpdir, emission_density, expected_mean
):
    """Check radial sampling for profiles on the coarsest valid grid."""
    major_radius = 620.0
    minor_radius = 200.0
    src = make_source(
        major_radius=major_radius,
        minor_radius=minor_radius,
        elongation=1.0,
        triangularity=0.0,
        shafranov_shift=0.0,
        r_over_a=[0.0, 1.0],
        emission_density=emission_density,
        energy=openmc.stats.delta_function(14.07e6),
    )

    sphere = openmc.Sphere(r=2000.0, boundary_type='vacuum')
    model = openmc.Model(
        geometry=openmc.Geometry([openmc.Cell(region=-sphere)]),
        settings=openmc.Settings(
            particles=100, batches=1, run_mode='fixed source', source=src),
    )

    sites = model.sample_external_source(20_000)
    xyz = np.array([site.r for site in sites])
    major_r = np.hypot(xyz[:, 0], xyz[:, 1])
    r_over_a = np.hypot(major_r - major_radius, xyz[:, 2]) / minor_radius
    assert_sample_mean(r_over_a, expected_mean)


def test_tokamak_source_poloidal_sampling(run_in_tmpdir):
    """Check linear-linear poloidal sampling on a coarse internal grid."""
    major_radius = 620.0
    minor_radius = 200.0
    with pytest.warns(UserWarning, match="below 51"):
        src = make_source(
            major_radius=major_radius,
            minor_radius=minor_radius,
            elongation=1.0,
            triangularity=0.0,
            shafranov_shift=0.0,
            emission_density=np.ones(10),
            n_alpha=3,
            energy=openmc.stats.delta_function(14.07e6),
        )

    sphere = openmc.Sphere(r=2000.0, boundary_type='vacuum')
    model = openmc.Model(
        geometry=openmc.Geometry([openmc.Cell(region=-sphere)]),
        settings=openmc.Settings(
            particles=100, batches=1, run_mode='fixed source', source=src),
    )

    sites = model.sample_external_source(100_000)
    xyz = np.array([site.r for site in sites])
    major_r = np.hypot(xyz[:, 0], xyz[:, 1])

    # With three alpha points, linear interpolation of cos(alpha) on [0, pi]
    # gives 1 - 2*alpha/pi. Integrating the resulting density gives this mean.
    expected_R = (
        major_radius
        + 2.0 * minor_radius**2 / (np.pi**2 * major_radius)
    )
    assert_sample_mean(major_r, expected_R)


@pytest.mark.flaky(reruns=1)
def test_tokamak_source_sampling(run_in_tmpdir):
    """Exercise the compiled C++ sampling path and check invariants.

    Sampled moments are compared against direct numerical quadrature of the
    exact source density S(r)*R*|J|, where the Jacobian J of the flux-surface
    map is computed from analytic partial derivatives. This check is
    independent of the Bernstein-mixture factorization used by the
    implementation.
    """
    R0, a, kappa, delta = 620.0, 200.0, 1.8, -0.5
    shafranov, zshift = 40.0, 25.0
    phi_start, phi_extent = 0.5, np.pi / 2

    # Fine grids to make discretization error negligible relative to
    # statistical uncertainty
    r_over_a = np.linspace(0.0, 1.0, 200)
    src = make_source(
        major_radius=R0, minor_radius=a, elongation=kappa,
        triangularity=delta, shafranov_shift=shafranov, vertical_shift=zshift,
        r_over_a=r_over_a, emission_density=1.0 - r_over_a**2,
        phi_start=phi_start, phi_extent=phi_extent, n_alpha=201,
        energy=openmc.stats.delta_function(14.07e6))

    sphere = openmc.Sphere(r=2000.0, boundary_type='vacuum')
    cell = openmc.Cell(region=-sphere)
    settings = openmc.Settings(
        particles=100, batches=1, run_mode='fixed source', source=src)
    model = openmc.Model(geometry=openmc.Geometry([cell]), settings=settings)

    n_samples = 20_000
    sites = model.sample_external_source(n_samples)

    xyz = np.array([s.r for s in sites])
    R = np.hypot(xyz[:, 0], xyz[:, 1])
    z = xyz[:, 2]

    # Energy, weight, and time invariants
    assert np.all([s.E == 14.07e6 for s in sites])
    assert np.all([s.wgt == 1.0 for s in sites])
    assert np.all([s.time == 0.0 for s in sites])

    # Toroidal angle within the requested sector
    phi = np.arctan2(xyz[:, 1], xyz[:, 0])
    assert phi.min() >= phi_start
    assert phi.max() <= phi_start + phi_extent

    # Positions bounded by the last closed flux surface
    assert R.min() >= R0 - a
    assert R.max() <= R0 + a + shafranov
    assert np.abs(z - zshift).max() <= kappa * a

    # Up-down symmetry about the vertical shift
    assert_sample_mean(z, zshift)

    # Reference moments by 2D quadrature of the exact density
    r = np.linspace(0.0, 1.0, 1001)[:, np.newaxis]  # r/a
    alpha = np.linspace(0.0, 2 * np.pi, 2001)[np.newaxis, :]
    psi = alpha + delta * np.sin(alpha)
    R_map = R0 + a * r * np.cos(psi) + shafranov * (1.0 - r**2)
    Z_map = kappa * a * r * np.sin(alpha)
    dR_dr = a * np.cos(psi) - 2.0 * shafranov * r
    dR_da = -a * r * np.sin(psi) * (1.0 + delta * np.cos(alpha))
    dZ_dr = kappa * a * np.sin(alpha) * np.ones_like(psi)
    dZ_da = kappa * a * r * np.cos(alpha)
    jac = np.abs(dR_dr * dZ_da - dR_da * dZ_dr)
    dens = (1.0 - r**2) * R_map * jac
    norm = dens.sum()
    expected_R = (R_map * dens).sum() / norm
    expected_z2 = (Z_map**2 * dens).sum() / norm

    assert_sample_mean(R, expected_R)
    assert_sample_mean((z - zshift)**2, expected_z2)
