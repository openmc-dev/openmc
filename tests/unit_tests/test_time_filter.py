from math import sqrt
from random import uniform

import numpy as np
import openmc
import pytest


def test_time_filter_basics():
    f = openmc.TimeFilter([1.0, 2.0, 5.0, 9.0])
    np.testing.assert_allclose(f.bins[0], [1.0, 2.0])
    np.testing.assert_allclose(f.bins[1], [2.0, 5.0])
    np.testing.assert_allclose(f.bins[-1], [5.0, 9.0])
    assert len(f.bins) == 3

    # Make sure __repr__ works
    repr(f)

    # to_xml_element()
    elem = f.to_xml_element()
    assert elem.tag == 'filter'
    assert elem.attrib['type'] == 'time'


def time(particle, distance, E):
    """Return the time it takes a particle at a given energy to travel a certain
    distance"""
    if particle == 'neutron':
        mass = 939.56542052e6  # eV/c²
    elif particle == 'photon':
        mass = 0.0

    # Calculate speed via v = c * sqrt(1 - γ^-2)
    inv_gamma = mass / (E + mass)
    velocity = 2.99792458e10 * sqrt(1 - inv_gamma * inv_gamma)  # cm/s
    return distance / velocity


@pytest.fixture(params=['neutron', 'photon'])
def model(request):
    # Select random sphere radius, source position, and source energy
    r = uniform(0., 2.)
    x = uniform(2., 10.)
    E = uniform(0., 20.0e6)

    # Create model
    model = openmc.Model()
    mat = openmc.Material()
    mat.add_nuclide('Zr90', 1.0)
    mat.set_density('g/cm3', 1.0)
    model.materials.append(mat)
    inner_sphere = openmc.Sphere(r=r)
    outer_sphere = openmc.Sphere(r=10.0, boundary_type='vacuum')
    center_cell = openmc.Cell(fill=mat, region=-inner_sphere)
    outer_void = openmc.Cell(region=+inner_sphere & -outer_sphere)
    model.geometry = openmc.Geometry([center_cell, outer_void])
    model.settings = openmc.Settings()
    model.settings.run_mode = 'fixed source'
    model.settings.particles = 1000
    model.settings.batches = 20
    particle = request.param
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point((x, 0., 0.)),
        angle=openmc.stats.Monodirectional([-1., 0., 0.]),
        energy=openmc.stats.Discrete([E], [1.0]),
        particle=particle
    )

    # Calculate time it will take neutrons to reach sphere
    t0 = time(particle, x - r, E)

    # Create tally with time filter
    tally = openmc.Tally()
    tally.filters = [openmc.TimeFilter([0.0, t0, 2*t0])]
    tally.scores = ['total']
    model.tallies.append(tally)
    return model


@pytest.fixture(params=['neutron', 'photon'])
def model_surf(request):
    # Select random distance and source energy
    x = uniform(50., 100.)
    E = uniform(0., 20.0e6)

    # Create model
    model = openmc.Model()
    mat = openmc.Material()
    mat.add_nuclide('Zr90', 1.0)
    mat.set_density('g/cm3', 1.0)
    model.materials.append(mat)
    left = openmc.XPlane(-1., boundary_type='vacuum')
    black_surface = openmc.XPlane(x, boundary_type='vacuum')
    right = openmc.XPlane(x + 1)
    void_cell = openmc.Cell(region=+left & -black_surface)
    black_cell = openmc.Cell(region=+black_surface & -right)
    model.geometry = openmc.Geometry([void_cell, black_cell])
    model.settings = openmc.Settings()
    model.settings.run_mode = 'fixed source'
    model.settings.particles = 1000
    model.settings.batches = 20
    particle = request.param
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point((0., 0., 0.)),
        angle=openmc.stats.Monodirectional([1., 0., 0.]),
        energy=openmc.stats.Discrete([E], [1.0]),
        particle=particle
    )

    # Calculate time it will take neutrons to reach purely-absorbing surface
    t0 = time(particle, x, E)

    # Create tally with surface and time filters
    tally = openmc.Tally()
    tally.filters = [
        openmc.SurfaceFilter([black_surface]),
        openmc.TimeFilter([0.0, t0*0.999, t0*1.001, 100.0])
    ]
    tally.scores = ['current']
    model.tallies.append(tally)
    return model


def test_time_filter_volume(model, run_in_tmpdir):
    sp_filename = model.run()
    with openmc.StatePoint(sp_filename) as sp:
        t = sp.tallies[model.tallies[0].id]
        values = t.mean.ravel()

        # Before t0, the reaction rate should be zero
        assert values[0] == pytest.approx(0.0)

        # After t0, the reaction rate should be positive
        assert values[1] > 0.0


def test_time_filter_surface(model_surf, run_in_tmpdir):
    sp_filename = model_surf.run()
    with openmc.StatePoint(sp_filename) as sp:
        t = sp.tallies[model_surf.tallies[0].id]
        values = t.mean.ravel()
        print(values)

        # Before t0-ε, the current should be zero
        assert values[0] == 0.0

        # Between t0-ε and t0+ε, the current should be one
        assert values[1] == 1.0

        # After t0+ε, the current should be zero
        assert values[2] == 0.0


def test_small_time_interval(run_in_tmpdir):
    # Create a model with a photon source at 1.0e8 seconds. Based on the speed
    # of the photon, the time intervals are on the order of 1e-9 seconds, which
    # are effectively 0 when compared to the starting time of the photon.
    mat = openmc.Material()
    mat.add_element('N', 1.0)
    mat.set_density('g/cm3', 0.001)
    sph = openmc.Sphere(r=5.0, boundary_type='vacuum')
    cell = openmc.Cell(fill=mat, region=-sph)
    model = openmc.Model()
    model.geometry = openmc.Geometry([cell])
    model.settings.particles = 100
    model.settings.batches = 10
    model.settings.run_mode = 'fixed source'
    model.settings.source = openmc.IndependentSource(
        time=openmc.stats.Discrete([1.0e8], [1.0]),
        particle='photon'
    )

    # Add tallies with and without a time filter that should match all particles
    time_filter = openmc.TimeFilter([0.0, 1.0e100])
    tally_with_filter = openmc.Tally()
    tally_with_filter.filters = [time_filter]
    tally_with_filter.scores = ['flux']
    tally_without_filter = openmc.Tally()
    tally_without_filter.scores = ['flux']
    model.tallies.extend([tally_with_filter, tally_without_filter])

    # Run the model and make sure the two tallies match
    sp_filename = model.run()
    with openmc.StatePoint(sp_filename) as sp:
        flux_with = sp.tallies[tally_with_filter.id].mean.ravel()[0]
        flux_without = sp.tallies[tally_without_filter.id].mean.ravel()[0]
        assert flux_with == pytest.approx(flux_without)
