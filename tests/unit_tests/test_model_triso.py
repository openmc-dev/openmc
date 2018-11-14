#!/usr/bin/env python

from math import pi

import numpy as np
from numpy.linalg import norm
import openmc
import openmc.model
import pytest
import scipy.spatial


_RADIUS = 0.1
_PACKING_FRACTION = 0.35
_SHAPES = ['rectangular_prism', 'cylinder', 'sphere', 'spherical_shell']
_VOLUMES = [1.**3, 1.*pi*1.**2, 4/3*pi*1.**3, 4/3*pi*(1.**3 - 0.5**3)]


@pytest.fixture(scope='module')
def container(request):
    return request.getfixturevalue(request.param)


@pytest.fixture(scope='module')
def centers(request):
    container = request.getfixturevalue(request.param)
    return openmc.model.pack_spheres(radius=_RADIUS, region=container,
        pf=_PACKING_FRACTION, initial_pf=0.2)


@pytest.fixture(scope='module')
def rectangular_prism():
    min_x = openmc.XPlane(x0=-0.5)
    max_x = openmc.XPlane(x0=0.5)
    min_y = openmc.YPlane(y0=-0.5)
    max_y = openmc.YPlane(y0=0.5)
    min_z = openmc.ZPlane(z0=-0.5)
    max_z = openmc.ZPlane(z0=0.5)
    return +min_x & -max_x & +min_y & -max_y & +min_z & -max_z


@pytest.fixture(scope='module')
def cylinder():
    cylinder = openmc.ZCylinder(R=1.)
    min_z = openmc.ZPlane(z0=-0.5)
    max_z = openmc.ZPlane(z0=0.5)
    return +min_z & -max_z & -cylinder


@pytest.fixture(scope='module')
def sphere():
    sphere = openmc.Sphere(R=1.)
    return -sphere


@pytest.fixture(scope='module')
def spherical_shell():
    sphere = openmc.Sphere(R=1.)
    inner_sphere = openmc.Sphere(R=0.5)
    return -sphere & +inner_sphere


@pytest.fixture(scope='module')
def triso_universe():
    sphere = openmc.Sphere(R=_RADIUS)
    cell = openmc.Cell(region=-sphere)
    univ = openmc.Universe(cells=[cell])
    return univ


@pytest.mark.parametrize('centers', _SHAPES, indirect=True)
def test_overlap(centers):
    """Check that none of the spheres in the packed configuration overlap."""
    # Create KD tree for quick nearest neighbor search
    tree = scipy.spatial.cKDTree(centers)

    # Find distance to nearest neighbor for all spheres
    d = tree.query(centers, k=2)[0]

    # Get the smallest distance between any two spheres
    d_min = min(d[:, 1])
    assert d_min > 2*_RADIUS or d_min == pytest.approx(2*_RADIUS)


@pytest.mark.parametrize('centers,shape', zip(_SHAPES, _SHAPES),
                         indirect=['centers'])
def test_contained(centers, shape):
    """Make sure all spheres are entirely contained within the domain."""
    if shape == 'rectangular_prism':
        x = np.amax(abs(centers)) + _RADIUS
        assert x < 0.5 or x == pytest.approx(0.5)

    elif shape == 'cylinder':
        r = max(np.linalg.norm(centers[:,0:2], axis=1)) + _RADIUS
        z = max(abs(centers[:,2])) + _RADIUS
        assert r < 1. or r == pytest.approx(1.)
        assert z < 0.5 or z == pytest.approx(0.5)

    elif shape == 'sphere':
        r = max(np.linalg.norm(centers, axis=1)) + _RADIUS
        assert r < 1. or r == pytest.approx(1.)

    elif shape == 'spherical_shell':
        d = np.linalg.norm(centers, axis=1)
        r_max = max(d) + _RADIUS
        r_min = min(d) - _RADIUS
        assert r_max < 1. or r_max == pytest.approx(1.)
        assert r_min > 0.5 or r_min == pytest.approx(0.5)


@pytest.mark.parametrize('centers,volume', zip(_SHAPES, _VOLUMES),
                         indirect=['centers'])
def test_packing_fraction(centers, volume):
    """Check that the actual PF is close to the requested PF."""
    pf = len(centers) * 4/3 * pi *_RADIUS**3 / volume
    assert pf == pytest.approx(_PACKING_FRACTION, rel=1e-2)


@pytest.mark.parametrize('container', _SHAPES, indirect=True)
def test_num_spheres(container):
    """Check that the function returns the correct number of spheres"""
    centers = openmc.model.pack_spheres(
        radius=_RADIUS, region=container, num_spheres=50
    )
    assert len(centers) == 50
    

def test_triso_lattice(triso_universe, rectangular_prism):
    centers = openmc.model.pack_spheres(
        radius=_RADIUS, region=rectangular_prism, pf=0.2
    )
    trisos = [openmc.model.TRISO(_RADIUS, triso_universe, c) for c in centers]

    lower_left = np.array((-.5, -.5, -.5))
    upper_right = np.array((.5, .5, .5))
    shape = (3, 3, 3)
    pitch = (upper_right - lower_left)/shape
    background = openmc.Material()

    lattice = openmc.model.create_triso_lattice(
        trisos, lower_left, pitch, shape, background
    )


def test_container_input(triso_universe):
    # Invalid container shape
    with pytest.raises(ValueError):
        centers = openmc.model.pack_spheres(
            radius=_RADIUS, region=+openmc.Sphere(R=1.), num_spheres=100
        )


def test_packing_fraction_input(sphere):
    # Provide neither packing fraction nor number of spheres
    with pytest.raises(ValueError):
        centers = openmc.model.pack_spheres(
            radius=_RADIUS, region=sphere
        )

    # Specify a packing fraction that is too high for CRP
    with pytest.raises(ValueError):
        centers = openmc.model.pack_spheres(
            radius=_RADIUS, region=sphere, pf=1.
        )

    # Specify a packing fraction that is too high for RSP
    with pytest.raises(ValueError):
        centers = openmc.model.pack_spheres(
            radius=_RADIUS, region=sphere, pf=0.5, initial_pf=0.4
        )
