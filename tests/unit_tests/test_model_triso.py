#!/usr/bin/env python

from math import pi

import numpy as np
from numpy.linalg import norm
import openmc
import openmc.model
import pytest
import scipy.spatial


_PACKING_FRACTION = 0.35
_RADIUS = 4.25e-2
domain_params = [
    {'shape': 'cube', 'length': 0.75, 'radius': 0., 'volume': 0.75**3},
    {'shape': 'cylinder', 'length': 0.5, 'radius': 0.5, 'volume': 0.5*pi*0.5**2},
    {'shape': 'sphere', 'length': 0., 'radius': 0.5, 'volume': 4/3*pi*0.5**3}
]


@pytest.fixture(scope='module', params=domain_params,
                ids=['cube', 'cylinder', 'sphere'])
def domain(request):
    return request.param


@pytest.fixture(scope='module')
def triso_universe():
    sphere = openmc.Sphere(R=_RADIUS)
    cell = openmc.Cell(region=-sphere)
    univ = openmc.Universe(cells=[cell])
    return univ


@pytest.fixture(scope='module')
def trisos(domain, triso_universe):
    trisos = openmc.model.pack_trisos(
        radius=_RADIUS,
        fill=triso_universe,
        domain_shape=domain['shape'],
        domain_length=domain['length'],
        domain_radius=domain['radius'],
        domain_center=(0., 0., 0.),
        initial_packing_fraction=0.2,
        packing_fraction=_PACKING_FRACTION
    )
    return trisos


def test_overlap(trisos):
    """Check that no TRISO particles overlap."""
    centers = [t.center for t in trisos]

    # Create KD tree for quick nearest neighbor search
    tree = scipy.spatial.cKDTree(centers)

    # Find distance to nearest neighbor for all particles
    d = tree.query(centers, k=2)[0]

    # Get the smallest distance between any two particles
    d_min = min(d[:, 1])
    assert d_min > 2*_RADIUS or d_min == pytest.approx(2*_RADIUS)


def test_contained(trisos, domain):
    """Make sure all particles are entirely contained within the domain."""
    if domain['shape'] == 'cube':
        x = max(np.hstack([abs(t.center) for t in trisos])) + _RADIUS
        assert x < 0.5*domain['length'] or x == pytest.approx(0.5*domain['length'])

    elif domain['shape'] == 'cylinder':
        r = max([norm(t.center[0:2]) for t in trisos]) + _RADIUS
        z = max([abs(t.center[2]) for t in trisos]) + _RADIUS
        assert r < domain['radius'] or r == pytest.approx(domain['radius'])
        assert z < 0.5*domain['length'] or z == pytest.approx(0.5*domain['length'])

    elif domain['shape'] == 'sphere':
        r = max([norm(t.center) for t in trisos]) + _RADIUS
        assert r < domain['radius'] or r == pytest.approx(domain['radius'])


def test_packing_fraction(trisos, domain):
    """Check that the actual PF is close to the requested PF."""
    pf = len(trisos)*4/3*pi*_RADIUS**3/domain['volume']
    assert pf == pytest.approx(_PACKING_FRACTION, rel=1e-2)


def test_n_particles(triso_universe):
    """Check that the function returns the correct number of particles"""
    trisos = openmc.model.pack_trisos(
        radius=_RADIUS, fill=triso_universe, domain_shape='cube',
        domain_length=1.0, n_particles=800
    )
    assert len(trisos) == 800
    

def test_triso_lattice(triso_universe):
    trisos = openmc.model.pack_trisos(
        radius=_RADIUS, fill=triso_universe, domain_shape='cube',
        domain_length=1.0, domain_center=(0., 0., 0.), packing_fraction=0.2
    )

    lower_left = np.array((-.5, -.5, -.5))
    upper_right = np.array((.5, .5, .5))
    shape = (3, 3, 3)
    pitch = (upper_right - lower_left)/shape
    background = openmc.Material()

    lattice = openmc.model.create_triso_lattice(
        trisos, lower_left, pitch, shape, background
    )


def test_domain_input(triso_universe):
    # Invalid domain shape
    with pytest.raises(ValueError):
        trisos = openmc.model.pack_trisos(
            radius=1, fill=triso_universe, n_particles=100,
            domain_shape='circle'
        )
    # Don't specify domain length on a cube
    with pytest.raises(ValueError):
        trisos = openmc.model.pack_trisos(
            radius=1, fill=triso_universe, n_particles=100,
            domain_shape='cube'
        )
    # Don't specify domain radius on a sphere
    with pytest.raises(ValueError):
        trisos = openmc.model.pack_trisos(
            radius=1, fill=triso_universe, n_particles=100,
            domain_shape='sphere'
        )
    

def test_packing_fraction_input(triso_universe):
    # Provide neither packing fraction nor number of particles
    with pytest.raises(ValueError):
        trisos = openmc.model.pack_trisos(
            radius=1, fill=triso_universe, domain_shape='cube',
            domain_length=10
        )
    # Provide both packing fraction and number of particles
    with pytest.raises(ValueError):
        trisos = openmc.model.pack_trisos(
            radius=1, fill=triso_universe, domain_shape='cube',
            domain_length=10, n_particles=100, packing_fraction=0.2
        )
    # Specify a packing fraction that is too high for CRP
    with pytest.raises(ValueError):
        trisos = openmc.model.pack_trisos(
            radius=1, fill=triso_universe, domain_shape='cube',
            domain_length=10, packing_fraction=1
        )
    # Specify a packing fraction that is too high for RSP
    with pytest.raises(ValueError):
        trisos = openmc.model.pack_trisos(
            radius=1, fill=triso_universe, domain_shape='cube',
            domain_length=10, packing_fraction=0.5,
            initial_packing_fraction=0.4
        )
