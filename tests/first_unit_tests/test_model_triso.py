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
_PARAMS = [
    {'shape': 'rectangular_prism', 'volume': 1**3},
    {'shape': 'x_cylinder', 'volume': 1*pi*1**2},
    {'shape': 'y_cylinder', 'volume': 1*pi*1**2},
    {'shape': 'z_cylinder', 'volume': 1*pi*1**2},
    {'shape': 'sphere', 'volume': 4/3*pi*1**3},
    {'shape': 'spherical_shell', 'volume': 4/3*pi*(1**3 - 0.5**3)}
]


@pytest.fixture(scope='module', params=_PARAMS)
def container(request):
    return request.param


@pytest.fixture(scope='module')
def centers(request, container):
    return request.getfixturevalue('centers_' + container['shape'])


@pytest.fixture(scope='module')
def centers_rectangular_prism():
    min_x = openmc.XPlane(0)
    max_x = openmc.XPlane(1)
    min_y = openmc.YPlane(0)
    max_y = openmc.YPlane(1)
    min_z = openmc.ZPlane(0)
    max_z = openmc.ZPlane(1)
    region = +min_x & -max_x & +min_y & -max_y & +min_z & -max_z
    return openmc.model.pack_spheres(radius=_RADIUS, region=region,
        pf=_PACKING_FRACTION, initial_pf=0.2)


@pytest.fixture(scope='module')
def centers_x_cylinder():
    cylinder = openmc.XCylinder(r=1, y0=1, z0=2)
    min_x = openmc.XPlane(0)
    max_x = openmc.XPlane(1)
    region = +min_x & -max_x & -cylinder
    return openmc.model.pack_spheres(radius=_RADIUS, region=region,
        pf=_PACKING_FRACTION, initial_pf=0.2)


@pytest.fixture(scope='module')
def centers_y_cylinder():
    cylinder = openmc.YCylinder(r=1, x0=1, z0=2)
    min_y = openmc.YPlane(0)
    max_y = openmc.YPlane(1)
    region = +min_y & -max_y & -cylinder
    return openmc.model.pack_spheres(radius=_RADIUS, region=region,
        pf=_PACKING_FRACTION, initial_pf=0.2)


@pytest.fixture(scope='module')
def centers_z_cylinder():
    cylinder = openmc.ZCylinder(r=1, x0=1, y0=2)
    min_z = openmc.ZPlane(0)
    max_z = openmc.ZPlane(1)
    region = +min_z & -max_z & -cylinder
    return openmc.model.pack_spheres(radius=_RADIUS, region=region,
        pf=_PACKING_FRACTION, initial_pf=0.2)


@pytest.fixture(scope='module')
def centers_sphere():
    sphere = openmc.Sphere(r=1, x0=1, y0=2, z0=3)
    region = -sphere
    return openmc.model.pack_spheres(radius=_RADIUS, region=region,
        pf=_PACKING_FRACTION, initial_pf=0.2)


@pytest.fixture(scope='module')
def centers_spherical_shell():
    sphere = openmc.Sphere(r=1, x0=1, y0=2, z0=3)
    inner_sphere = openmc.Sphere(r=0.5, x0=1, y0=2, z0=3)
    region = -sphere & +inner_sphere
    return openmc.model.pack_spheres(radius=_RADIUS, region=region,
        pf=_PACKING_FRACTION, initial_pf=0.2)


@pytest.fixture(scope='module')
def triso_universe():
    sphere = openmc.Sphere(r=_RADIUS)
    cell = openmc.Cell(region=-sphere)
    univ = openmc.Universe(cells=[cell])
    return univ


def test_overlap(centers):
    """Check that none of the spheres in the packed configuration overlap."""
    # Create KD tree for quick nearest neighbor search
    tree = scipy.spatial.cKDTree(centers)

    # Find distance to nearest neighbor for all spheres
    d = tree.query(centers, k=2)[0]

    # Get the smallest distance between any two spheres
    d_min = min(d[:, 1])
    assert d_min > 2*_RADIUS or d_min == pytest.approx(2*_RADIUS)


def test_contained_rectangular_prism(centers_rectangular_prism):
    """Make sure all spheres are entirely contained within the domain."""
    d_max = np.amax(centers_rectangular_prism) + _RADIUS
    d_min = np.amin(centers_rectangular_prism) - _RADIUS
    assert d_max < 1 or d_max == pytest.approx(1)
    assert d_min > 0 or d_min == pytest.approx(0)


def test_contained_x_cylinder(centers_x_cylinder):
    """Make sure all spheres are entirely contained within the domain."""
    d = np.linalg.norm(centers_x_cylinder[:,[1,2]] - [1, 2], axis=1)
    r_max = max(d) + _RADIUS
    x_max = max(centers_x_cylinder[:,0]) + _RADIUS
    x_min = min(centers_x_cylinder[:,0]) - _RADIUS
    assert r_max < 1 or r_max == pytest.approx(1)
    assert x_max < 1 or x_max == pytest.approx(1)
    assert x_min > 0 or x_min == pytest.approx(0)


def test_contained_y_cylinder(centers_y_cylinder):
    """Make sure all spheres are entirely contained within the domain."""
    d = np.linalg.norm(centers_y_cylinder[:,[0,2]] - [1, 2], axis=1)
    r_max = max(d) + _RADIUS
    y_max = max(centers_y_cylinder[:,1]) + _RADIUS
    y_min = min(centers_y_cylinder[:,1]) - _RADIUS
    assert r_max < 1 or r_max == pytest.approx(1)
    assert y_max < 1 or y_max == pytest.approx(1)
    assert y_min > 0 or y_min == pytest.approx(0)


def test_contained_z_cylinder(centers_z_cylinder):
    """Make sure all spheres are entirely contained within the domain."""
    d = np.linalg.norm(centers_z_cylinder[:,[0,1]] - [1, 2], axis=1)
    r_max = max(d) + _RADIUS
    z_max = max(centers_z_cylinder[:,2]) + _RADIUS
    z_min = min(centers_z_cylinder[:,2]) - _RADIUS
    assert r_max < 1 or r_max == pytest.approx(1)
    assert z_max < 1 or z_max == pytest.approx(1)
    assert z_min > 0 or z_min == pytest.approx(0)


def test_contained_sphere(centers_sphere):
    """Make sure all spheres are entirely contained within the domain."""
    d = np.linalg.norm(centers_sphere - [1, 2, 3], axis=1)
    r_max = max(d) + _RADIUS
    assert r_max < 1 or r_max == pytest.approx(1)


def test_contained_spherical_shell(centers_spherical_shell):
    """Make sure all spheres are entirely contained within the domain."""
    d = np.linalg.norm(centers_spherical_shell - [1, 2, 3], axis=1)
    r_max = max(d) + _RADIUS
    r_min = min(d) - _RADIUS
    assert r_max < 1 or r_max == pytest.approx(1)
    assert r_min > 0.5 or r_min == pytest.approx(0.5)


def test_packing_fraction(container, centers):
    """Check that the actual PF is close to the requested PF."""
    pf = len(centers) * 4/3 * pi *_RADIUS**3 / container['volume']
    assert pf == pytest.approx(_PACKING_FRACTION, rel=1e-2)


def test_num_spheres():
    """Check that the function returns the correct number of spheres"""
    centers = openmc.model.pack_spheres(
        radius=_RADIUS, region=-openmc.Sphere(r=1), num_spheres=50
    )
    assert len(centers) == 50


def test_triso_lattice(triso_universe, centers_rectangular_prism):
    trisos = [openmc.model.TRISO(_RADIUS, triso_universe, c)
              for c in centers_rectangular_prism]

    lower_left = np.array((0, 0, 0))
    upper_right = np.array((1, 1, 1))
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
            radius=_RADIUS, region=+openmc.Sphere(r=1), num_spheres=100
        )


def test_packing_fraction_input():
    # Provide neither packing fraction nor number of spheres
    with pytest.raises(ValueError):
        centers = openmc.model.pack_spheres(
            radius=_RADIUS, region=-openmc.Sphere(r=1)
        )

    # Specify a packing fraction that is too high for CRP
    with pytest.raises(ValueError):
        centers = openmc.model.pack_spheres(
            radius=_RADIUS, region=-openmc.Sphere(r=1), pf=1
        )

    # Specify a packing fraction that is too high for RSP
    with pytest.raises(ValueError):
        centers = openmc.model.pack_spheres(
            radius=_RADIUS, region=-openmc.Sphere(r=1), pf=0.5, initial_pf=0.4
        )
