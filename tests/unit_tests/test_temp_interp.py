from math import isnan
import os
from pathlib import Path

import numpy as np
import openmc.data
from openmc.data import K_BOLTZMANN
from openmc.stats import Uniform
import pytest


def make_fake_cross_section():
    """Create fake U235 nuclide

    This nuclide is designed to have k_inf=1 at 300 K, k_inf=2 at 600 K, and
    k_inf=1 at 900 K. The absorption cross section is also constant with
    temperature so as to make the true k-effective go linear with temperature.
    """

    def isotropic_angle(E_min, E_max):
        return openmc.data.AngleDistribution(
            [E_min, E_max],
            [Uniform(-1., 1.), Uniform(-1., 1.)]
        )

    def cross_section(value):
        return openmc.data.Tabulated1D(
            energy,
            value*np.ones_like(energy)
        )

    temperatures = (300, 600, 900)

    u235_fake = openmc.data.IncidentNeutron(
        'U235', 92, 235, 0, 233.0248, [T*K_BOLTZMANN for T in temperatures]
    )

    # Create energy grids
    E_min, E_max = 1e-5, 20.0e6
    energy = np.logspace(np.log10(E_min), np.log10(E_max))
    for T in temperatures:
        u235_fake.energy['{}K'.format(T)] = energy

    # Create elastic scattering
    elastic = openmc.data.Reaction(2)
    for T in temperatures:
        elastic.xs['{}K'.format(T)] = cross_section(1.0)
    elastic_dist = openmc.data.UncorrelatedAngleEnergy(isotropic_angle(E_min, E_max))
    product = openmc.data.Product()
    product.distribution.append(elastic_dist)
    elastic.products.append(product)
    u235_fake.reactions[2] = elastic

    # Create fission
    fission = openmc.data.Reaction(18)
    fission.center_of_mass = False
    fission.Q_value = 193.0e6
    fission_xs = (2., 4., 2.)
    for T, xs in zip(temperatures, fission_xs):
        fission.xs['{}K'.format(T)] = cross_section(xs)
    a = openmc.data.Tabulated1D([E_min, E_max], [0.988e6, 0.988e6])
    b = openmc.data.Tabulated1D([E_min, E_max], [2.249e-6, 2.249e-6])
    fission_dist = openmc.data.UncorrelatedAngleEnergy(
        isotropic_angle(E_min, E_max),
        openmc.data.WattEnergy(a, b, -E_max)
    )
    product = openmc.data.Product()
    product.distribution.append(fission_dist)
    product.yield_ = openmc.data.Polynomial((2.0,))
    fission.products.append(product)
    u235_fake.reactions[18] = fission

    # Create capture
    capture = openmc.data.Reaction(102)
    capture.q_value = 6.5e6
    capture_xs = (2., 0., 2.)
    for T, xs in zip(temperatures, capture_xs):
        capture.xs['{}K'.format(T)] = cross_section(xs)
    u235_fake.reactions[102] = capture

    # Export HDF5 file
    u235_fake.export_to_hdf5('U235_fake.h5', 'w')

    lib = openmc.data.DataLibrary()
    lib.register_file('U235_fake.h5')
    lib.export_to_xml('cross_sections_fake.xml')


def fake_thermal_scattering():
    """Create a fake thermal scattering library for U-235 at 294K and 600K
    """
    fake_tsl = openmc.data.ThermalScattering("c_U_fake", 1.9968, 4.9, [0.0253])
    fake_tsl.nuclides = ['U235']

    # Create elastic reaction
    bragg_edges = [0.00370672, 0.00494229]
    factors = [0.00375735, 0.01386287]
    coherent_xs = openmc.data.CoherentElastic(bragg_edges, factors)
    incoherent_xs_294 = openmc.data.Tabulated1D([0.00370672, 0.00370672], [0.00370672, 0.00370672])
    elastic_xs_base = openmc.data.Sum((coherent_xs, incoherent_xs_294))
    elastic_xs = {'294K': elastic_xs_base, '600K': elastic_xs_base}
    coherent_dist = openmc.data.CoherentElasticAE(coherent_xs)
    incoherent_dist_294 = openmc.data.IncoherentElasticAEDiscrete([
        [-0.6, -0.18, 0.18, 0.6], [-0.6, -0.18, 0.18, 0.6]
    ])
    incoherent_dist_600 = openmc.data.IncoherentElasticAEDiscrete([
        [-0.1, -0.2, 0.2, 0.1], [-0.1, -0.2, 0.2, 0.1]
    ])
    elastic_dist = {
        '294K': openmc.data.MixedElasticAE(coherent_dist, incoherent_dist_294),
        '600K': openmc.data.MixedElasticAE(coherent_dist, incoherent_dist_600)
        }
    fake_tsl.elastic = openmc.data.ThermalScatteringReaction(elastic_xs, elastic_dist)

    # Create inelastic reaction
    inelastic_xs = {
        '294K': openmc.data.Tabulated1D([1.0e-5, 4.9], [13.4, 3.35]),
        '600K': openmc.data.Tabulated1D([1.0e-2, 10], [1.4, 5])
        }
    breakpoints = [3]
    interpolation = [2]
    energy = [1.0e-5, 4.3e-2, 4.9]
    energy_out = [
        openmc.data.Tabular([0.0002, 0.067, 0.146, 0.366], [0.25, 0.25, 0.25, 0.25]),
        openmc.data.Tabular([0.0001, 0.009, 0.137, 0.277], [0.25, 0.25, 0.25, 0.25]),
        openmc.data.Tabular([0.0579, 4.555, 4.803, 4.874], [0.25, 0.25, 0.25, 0.25]),
    ]
    for eout in energy_out:
        eout.normalize()
        eout.c = eout.cdf()
    discrete = openmc.stats.Discrete([-0.9, -0.6, -0.3, -0.1, 0.1, 0.3, 0.6, 0.9], [1/8]*8)
    discrete.c = discrete.cdf()[1:]
    mu = [[discrete]*4]*3
    dist = openmc.data.IncoherentInelasticAE(
        breakpoints, interpolation, energy, energy_out, mu)
    inelastic_dist = {'294K': dist, '600K': dist}
    inelastic = openmc.data.ThermalScatteringReaction(inelastic_xs, inelastic_dist)
    fake_tsl.inelastic = inelastic

    return fake_tsl


def edit_fake_cross_sections():
    """Edit the test cross sections xml to include fake thermal scattering data
    """
    lib = openmc.data.DataLibrary.from_xml("cross_sections_fake.xml")
    c_U_fake = fake_thermal_scattering()
    c_U_fake.export_to_hdf5("c_U_fake.h5")
    lib.register_file("c_U_fake.h5")
    lib.export_to_xml("cross_sections_fake.xml")


@pytest.fixture(scope='module')
def model(tmp_path_factory):
    tmp_path = tmp_path_factory.mktemp("temp_interp")
    orig = Path.cwd()
    os.chdir(tmp_path)

    make_fake_cross_section()

    model = openmc.model.Model()
    mat = openmc.Material()
    mat.add_nuclide('U235', 1.0)
    model.materials.append(mat)
    model.materials.cross_sections = str(Path('cross_sections_fake.xml').resolve())

    sph = openmc.Sphere(r=100.0, boundary_type='reflective')
    cell = openmc.Cell(fill=mat, region=-sph)
    model.geometry = openmc.Geometry([cell])

    model.settings.particles = 1000
    model.settings.inactive = 0
    model.settings.batches = 10

    tally = openmc.Tally()
    tally.scores = ['absorption', 'fission', 'scatter', 'nu-fission']
    model.tallies = [tally]

    try:
        yield model
    finally:
        os.chdir(orig)


@pytest.mark.parametrize(
    ["method", "temperature", "fission_expected", "tolerance"],
    [
        ("nearest", 300.0, 0.5, 10),
        ("nearest", 600.0, 1.0, 10),
        ("nearest", 900.0, 0.5, 10),
        ("interpolation", 360.0, 0.6, 10),
        ("interpolation", 450.0, 0.75, 10),
        ("interpolation", 540.0, 0.9, 10),
        ("interpolation", 660.0, 0.9, 10),
        ("interpolation", 750.0, 0.75, 10),
        ("interpolation", 840.0, 0.6, 10),
        ("interpolation", 295.0, 0.5, 10),
        ("interpolation", 990.0, 0.5, 100),
    ]
)
def test_interpolation(model, method, temperature, fission_expected, tolerance):
    model.settings.temperature = {'method': method, 'default': temperature, "tolerance": tolerance}
    sp_filename = model.run()
    with openmc.StatePoint(sp_filename) as sp:
        t = sp.tallies[model.tallies[0].id]
        absorption_mean, fission_mean, scatter_mean, nu_fission_mean = t.mean.ravel()
        absorption_unc, fission_unc, scatter_unc, nu_fission_unc = t.std_dev.ravel()

        nu = 2.0
        assert abs(absorption_mean - 1) < 3*absorption_unc
        assert abs(fission_mean - fission_expected) < 3*fission_unc
        assert abs(scatter_mean - 1/4) < 3*scatter_unc
        assert abs(nu_fission_mean - nu*fission_expected) < 3*nu_fission_unc

        # Check that k-effective value matches expected
        k = sp.keff
        if isnan(k.s):
            assert k.n == pytest.approx(nu*fission_expected)
        else:
            assert abs(k.n - nu*fission_expected) <= 3*k.s


def test_temperature_interpolation_tolerance(model):
    """Test applying global and cell temperatures with thermal scattering libraries
    """
    edit_fake_cross_sections()
    model.materials[0].add_s_alpha_beta("c_U_fake")

    # Default k-effective, using the thermal scattering data's minimum available temperature
    model.settings.temperature = {'method': "nearest", 'default': 294, "tolerance": 50}
    sp_filename = model.run()
    with openmc.StatePoint(sp_filename) as sp:
        default_k = sp.keff.n

    # Get k-effective with temperature below the minimum but in interpolation mode
    model.settings.temperature = {'method': "interpolation", 'default': 255, "tolerance": 50}
    sp_filename = model.run()
    with openmc.StatePoint(sp_filename) as sp:
        interpolated_k = sp.keff.n

    # Get the k-effective with the temperature applied to the cell, instead of globally
    model.settings.temperature = {'method': "interpolation", 'default': 500, "tolerance": 50}
    for cell in model.geometry.get_all_cells().values():
        cell.temperature = 275
    sp_filename = model.run()
    with openmc.StatePoint(sp_filename) as sp:
        cell_k = sp.keff.n

    # All calculated k-effectives should be equal
    assert default_k == pytest.approx(interpolated_k)
    assert interpolated_k == pytest.approx(cell_k)
