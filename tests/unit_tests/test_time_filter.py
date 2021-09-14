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


@pytest.fixture
def model():
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
    model.settings.source = openmc.Source(
        space=openmc.stats.Point((x, 0., 0.)),
        angle=openmc.stats.Monodirectional([-1., 0., 0.]),
        energy=openmc.stats.Discrete([E], [1.0])
    )

    # Source particles will take a time of t = x/v to reach the inner sphere,
    # where x is the distance from the source to the sphere and and v =
    # sqrt(2E/m). Before this time, there should be no reactions,
    MASS_NEUTRON_EV = 939.56542052e6  # eV/c²
    C_LIGHT = 2.99792458e10  # cm/s
    distance = x - r
    velocity = sqrt(2*E / MASS_NEUTRON_EV) * C_LIGHT
    t0 = distance / velocity

    # Create tally with time filter
    tally = openmc.Tally()
    tally.filters = [openmc.TimeFilter([0.0, t0, 2*t0])]
    tally.scores = ['total']
    model.tallies.append(tally)
    return model


def test_time_filter_transport(model, run_in_tmpdir):
    sp_filename = model.run()
    with openmc.StatePoint(sp_filename) as sp:
        t = sp.tallies[model.tallies[0].id]
        values = t.mean.ravel()

        # Before t0, the reaction rate should be zero
        assert values[0] == 0.0

        # After t0, the reaction rate should be positive
        assert values[1] > 0.0
