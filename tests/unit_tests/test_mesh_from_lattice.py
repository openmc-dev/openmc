import numpy as np
import openmc
import pytest


@pytest.fixture(scope='module')
def pincell1(uo2, water):
    cyl = openmc.ZCylinder(r=0.35)
    fuel = openmc.Cell(fill=uo2, region=-cyl)
    moderator = openmc.Cell(fill=water, region=+cyl)

    univ = openmc.Universe(cells=[fuel, moderator])
    univ.fuel = fuel
    univ.moderator = moderator
    return univ


@pytest.fixture(scope='module')
def pincell2(uo2, water):
    cyl = openmc.ZCylinder(r=0.4)
    fuel = openmc.Cell(fill=uo2, region=-cyl)
    moderator = openmc.Cell(fill=water, region=+cyl)

    univ = openmc.Universe(cells=[fuel, moderator])
    univ.fuel = fuel
    univ.moderator = moderator
    return univ


@pytest.fixture(scope='module')
def zr():
    zr = openmc.Material()
    zr.add_element('Zr', 1.0)
    zr.set_density('g/cm3', 1.0)
    return zr


@pytest.fixture(scope='module')
def rlat2(pincell1, pincell2, uo2, water, zr):
    """2D Rectangular lattice for testing."""
    all_zr = openmc.Cell(fill=zr)
    pitch = 1.2
    n = 3
    u1, u2 = pincell1, pincell2
    lattice = openmc.RectLattice()
    lattice.lower_left = (-pitch*n/2, -pitch*n/2)
    lattice.pitch = (pitch, pitch)
    lattice.outer = openmc.Universe(cells=[all_zr])
    lattice.universes = [
        [u1, u2, u1],
        [u2, u1, u2],
        [u2, u1, u1]
    ]

    return lattice


@pytest.fixture(scope='module')
def rlat3(pincell1, pincell2, uo2, water, zr):
    """3D Rectangular lattice for testing."""

    # Create another universe for top layer
    hydrogen = openmc.Material()
    hydrogen.add_element('H', 1.0)
    hydrogen.set_density('g/cm3', 0.09)
    h_cell = openmc.Cell(fill=hydrogen)
    u3 = openmc.Universe(cells=[h_cell])

    all_zr = openmc.Cell(fill=zr)
    pitch = 1.2
    n = 3
    u1, u2 = pincell1, pincell2
    lattice = openmc.RectLattice()
    lattice.lower_left = (-pitch*n/2, -pitch*n/2, -10.0)
    lattice.pitch = (pitch, pitch, 10.0)
    lattice.outer = openmc.Universe(cells=[all_zr])
    lattice.universes = [
        [[u1, u2, u1],
         [u2, u1, u2],
         [u2, u1, u1]],
        [[u3, u1, u2],
         [u1, u3, u2],
         [u2, u1, u1]]
    ]

    return lattice


def test_mesh2d(rlat2):
    shape = np.array(rlat2.shape)
    width = shape*rlat2.pitch

    mesh1 = openmc.RegularMesh.from_rect_lattice(rlat2)
    assert np.array_equal(mesh1.dimension, (3, 3))
    assert np.array_equal(mesh1.lower_left, rlat2.lower_left)
    assert np.array_equal(mesh1.upper_right, rlat2.lower_left + width)

    mesh2 = openmc.RegularMesh.from_rect_lattice(rlat2, division=3)
    assert np.array_equal(mesh2.dimension, (9, 9))
    assert np.array_equal(mesh2.lower_left, rlat2.lower_left)
    assert np.array_equal(mesh2.upper_right, rlat2.lower_left + width)


def test_mesh3d(rlat3):
    shape = np.array(rlat3.shape)
    width = shape*rlat3.pitch

    mesh1 = openmc.RegularMesh.from_rect_lattice(rlat3)
    assert np.array_equal(mesh1.dimension, (3, 3, 2))
    assert np.array_equal(mesh1.lower_left, rlat3.lower_left)
    assert np.array_equal(mesh1.upper_right, rlat3.lower_left + width)

    mesh2 = openmc.RegularMesh.from_rect_lattice(rlat3, division=3)
    assert np.array_equal(mesh2.dimension, (9, 9, 6))
    assert np.array_equal(mesh2.lower_left, rlat3.lower_left)
    assert np.array_equal(mesh2.upper_right, rlat3.lower_left + width)
