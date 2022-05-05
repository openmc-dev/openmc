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
    ["method", "temperature", "fission_expected"],
    [
        ("nearest", 300.0, 0.5),
        ("nearest", 600.0, 1.0),
        ("nearest", 900.0, 0.5),
        ("interpolation", 360.0, 0.6),
        ("interpolation", 450.0, 0.75),
        ("interpolation", 540.0, 0.9),
        ("interpolation", 660.0, 0.9),
        ("interpolation", 750.0, 0.75),
        ("interpolation", 840.0, 0.6),
    ]
)
def test_interpolation(model, method, temperature, fission_expected):
    model.settings.temperature = {'method': method, 'default': temperature}
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
