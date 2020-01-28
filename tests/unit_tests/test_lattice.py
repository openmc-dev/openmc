from math import sqrt
import xml.etree.ElementTree as ET

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

    # Add extra attributes for comparison purpose
    lattice.cells = [u1.fuel, u1.moderator, u2.fuel, u2.moderator, all_zr]
    lattice.mats = [uo2, water, zr]
    lattice.univs = [u1, u2, lattice.outer]
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

    # Add extra attributes for comparison purpose
    lattice.cells = [u1.fuel, u1.moderator, u2.fuel, u2.moderator,
                     h_cell, all_zr]
    lattice.mats = [uo2, water, zr, hydrogen]
    lattice.univs = [u1, u2, u3, lattice.outer]
    return lattice


@pytest.fixture(scope='module')
def hlat2(pincell1, pincell2, uo2, water, zr):
    """2D Hexagonal lattice for testing."""
    all_zr = openmc.Cell(fill=zr)
    pitch = 1.2
    u1, u2 = pincell1, pincell2
    lattice = openmc.HexLattice()
    lattice.center = (0., 0.)
    lattice.pitch = (pitch,)
    lattice.outer = openmc.Universe(cells=[all_zr])
    lattice.universes = [
        [u2, u1, u1, u1, u1, u1, u1, u1, u1, u1, u1, u1],
        [u2, u1, u1, u1, u1, u1],
        [u2]
    ]

    # Add extra attributes for comparison purpose
    lattice.cells = [u1.fuel, u1.moderator, u2.fuel, u2.moderator, all_zr]
    lattice.mats = [uo2, water, zr]
    lattice.univs = [u1, u2, lattice.outer]
    return lattice


@pytest.fixture(scope='module')
def hlat3(pincell1, pincell2, uo2, water, zr):
    """3D Hexagonal lattice for testing."""

    # Create another universe for top layer
    hydrogen = openmc.Material()
    hydrogen.add_element('H', 1.0)
    hydrogen.set_density('g/cm3', 0.09)
    h_cell = openmc.Cell(fill=hydrogen)
    u3 = openmc.Universe(cells=[h_cell])

    all_zr = openmc.Cell(fill=zr)
    pitch = 1.2
    u1, u2 = pincell1, pincell2
    lattice = openmc.HexLattice()
    lattice.center = (0., 0., 0.)
    lattice.pitch = (pitch, 10.0)
    lattice.outer = openmc.Universe(cells=[all_zr])
    lattice.universes = [
        [[u2, u1, u1, u1, u1, u1, u1, u1, u1, u1, u1, u1],
         [u2, u1, u1, u1, u1, u1],
         [u2]],
        [[u1, u1, u1, u1, u1, u1, u3, u1, u1, u1, u1, u1],
         [u1, u1, u1, u3, u1, u1],
         [u3]]
    ]

    # Add extra attributes for comparison purpose
    lattice.cells = [u1.fuel, u1.moderator, u2.fuel, u2.moderator,
                     h_cell, all_zr]
    lattice.mats = [uo2, water, zr, hydrogen]
    lattice.univs = [u1, u2, u3, lattice.outer]
    return lattice


def test_get_nuclides(rlat2, rlat3, hlat2, hlat3):
    for lat in (rlat2, hlat2):
        nucs = rlat2.get_nuclides()
        assert sorted(nucs) == ['H1', 'O16', 'U235',
                                'Zr90', 'Zr91', 'Zr92', 'Zr94', 'Zr96']
    for lat in (rlat3, hlat3):
        nucs = rlat3.get_nuclides()
        assert sorted(nucs) == ['H1', 'H2', 'O16', 'U235',
                                'Zr90', 'Zr91', 'Zr92', 'Zr94', 'Zr96']


def test_get_all_cells(rlat2, rlat3, hlat2, hlat3):
    for lat in (rlat2, rlat3, hlat2, hlat3):
        cells = set(lat.get_all_cells().values())
        assert not cells ^ set(lat.cells)


def test_get_all_materials(rlat2, rlat3, hlat2, hlat3):
    for lat in (rlat2, rlat3, hlat2, hlat3):
        mats = set(lat.get_all_materials().values())
        assert not mats ^ set(lat.mats)


def test_get_all_universes(rlat2, rlat3, hlat2, hlat3):
    for lat in (rlat2, rlat3, hlat2, hlat3):
        univs = set(lat.get_all_universes().values())
        assert not univs ^ set(lat.univs)


def test_get_universe(rlat2, rlat3, hlat2, hlat3):
    u1, u2, outer = rlat2.univs
    assert rlat2.get_universe((0, 0)) == u2
    assert rlat2.get_universe((1, 0)) == u1
    assert rlat2.get_universe((0, 1)) == u2

    u1, u2, u3, outer = rlat3.univs
    assert rlat3.get_universe((0, 0, 0)) == u2
    assert rlat3.get_universe((2, 2, 0)) == u1
    assert rlat3.get_universe((0, 2, 1)) == u3
    assert rlat3.get_universe((2, 1, 1)) == u2

    u1, u2, outer = hlat2.univs
    assert hlat2.get_universe((0, 0)) == u2
    assert hlat2.get_universe((0, 2)) == u2
    assert hlat2.get_universe((1, 0)) == u1
    assert hlat2.get_universe((-2, 2)) == u1

    hlat2.orientation = 'x'
    assert hlat2.get_universe((2, 0)) == u2
    assert hlat2.get_universe((1, 0)) == u2
    assert hlat2.get_universe((1, 1)) == u1
    assert hlat2.get_universe((-1, 1)) == u1
    hlat2.orientation = 'y'

    u1, u2, u3, outer = hlat3.univs
    assert hlat3.get_universe((0, 0, 0)) == u2
    assert hlat3.get_universe((0, 0, 1)) == u3
    assert hlat3.get_universe((0, 2, 0)) == u2
    assert hlat3.get_universe((0, 2, 1)) == u1
    assert hlat3.get_universe((0, -2, 0)) == u1
    assert hlat3.get_universe((0, -2, 1)) == u3


def test_find(rlat2, rlat3, hlat2, hlat3):
    pitch = rlat2.pitch[0]
    seq = rlat2.find((0., 0., 0.))
    assert seq[-1] == rlat2.cells[0]
    seq = rlat2.find((pitch, 0., 0.))
    assert seq[-1] == rlat2.cells[2]
    seq = rlat2.find((0., -pitch, 0.))
    assert seq[-1] == rlat2.cells[0]
    seq = rlat2.find((pitch*100, 0., 0.))
    assert seq[-1] == rlat2.cells[-1]
    seq = rlat3.find((-pitch, pitch, 5.0))
    assert seq[-1] == rlat3.cells[-2]

    pitch = hlat2.pitch[0]
    seq = hlat2.find((0., 0., 0.))
    assert seq[-1] == hlat2.cells[2]
    seq = hlat2.find((0.5, 0., 0.))
    assert seq[-1] == hlat2.cells[3]
    seq = hlat2.find((sqrt(3)*pitch, 0., 0.))
    assert seq[-1] == hlat2.cells[0]
    seq = hlat2.find((0., pitch, 0.))
    assert seq[-1] == hlat2.cells[2]

    # bottom of 3D lattice
    seq = hlat3.find((0., 0., -5.))
    assert seq[-1] == hlat3.cells[2]
    seq = hlat3.find((0., pitch, -5.))
    assert seq[-1] == hlat3.cells[2]
    seq = hlat3.find((0., -pitch, -5.))
    assert seq[-1] == hlat3.cells[0]
    seq = hlat3.find((sqrt(3)*pitch, 0., -5.))
    assert seq[-1] == hlat3.cells[0]

    # top of 3D lattice
    seq = hlat3.find((0., 0., 5.))
    assert seq[-1] == hlat3.cells[-2]
    seq = hlat3.find((0., pitch, 5.))
    assert seq[-1] == hlat3.cells[0]
    seq = hlat3.find((0., -pitch, 5.))
    assert seq[-1] == hlat3.cells[-2]
    seq = hlat3.find((sqrt(3)*pitch, 0., 5.))
    assert seq[-1] == hlat3.cells[0]


def test_clone(rlat2, hlat2, hlat3):
    rlat_clone = rlat2.clone()
    assert rlat_clone.id != rlat2.id
    assert rlat_clone.lower_left == rlat2.lower_left
    assert rlat_clone.pitch == rlat2.pitch

    hlat_clone = hlat2.clone()
    assert hlat_clone.id != hlat2.id
    assert hlat_clone.center == hlat2.center
    assert hlat_clone.pitch == hlat2.pitch

    hlat_clone = hlat3.clone()
    assert hlat_clone.id != hlat3.id
    assert hlat_clone.center == hlat3.center
    assert hlat_clone.pitch == hlat3.pitch

    rlat_clone = rlat2.clone(clone_materials=False)
    assert rlat_clone.get_all_materials() == rlat2.get_all_materials()

    rlat_clone = rlat2.clone(clone_materials=False, clone_regions=False)
    for c1 in rlat_clone.cells:
        for c2 in rlat2.cells:
            if c1.fill == c2.fill:
                print(c1.fill)
                assert c1.region == c2.region


def test_repr(rlat2, rlat3, hlat2, hlat3):
    repr(rlat2)
    repr(rlat3)
    repr(hlat2)
    repr(hlat3)


def test_indices_rect(rlat2, rlat3):
    # (y, x) indices
    assert rlat2.indices == [(0, 0), (0, 1), (0, 2),
                             (1, 0), (1, 1), (1, 2),
                             (2, 0), (2, 1), (2, 2)]
    # (z, y, x) indices
    assert rlat3.indices == [
        (0, 0, 0), (0, 0, 1), (0, 0, 2),
        (0, 1, 0), (0, 1, 1), (0, 1, 2),
        (0, 2, 0), (0, 2, 1), (0, 2, 2),
        (1, 0, 0), (1, 0, 1), (1, 0, 2),
        (1, 1, 0), (1, 1, 1), (1, 1, 2),
        (1, 2, 0), (1, 2, 1), (1, 2, 2)
    ]


def test_indices_hex(hlat2, hlat3):
    # (r, i) indices
    assert hlat2.indices == (
        [(0, i) for i in range(12)] +
        [(1, i) for i in range(6)] +
        [(2, 0)]
    )

    # (z, r, i) indices
    assert hlat3.indices == (
        [(0, 0, i) for i in range(12)] +
        [(0, 1, i) for i in range(6)] +
        [(0, 2, 0)] +
        [(1, 0, i) for i in range(12)] +
        [(1, 1, i) for i in range(6)] +
        [(1, 2, 0)]
    )


def test_xml_rect(rlat2, rlat3):
    for lat in (rlat2, rlat3):
        geom = ET.Element('geometry')
        lat.create_xml_subelement(geom)
        elem = geom.find('lattice')
        assert elem.tag == 'lattice'
        assert elem.get('id') == str(lat.id)
        assert len(elem.find('pitch').text.split()) == lat.ndim
        assert len(elem.find('lower_left').text.split()) == lat.ndim
        assert len(elem.find('universes').text.split()) == len(lat.indices)


def test_xml_hex(hlat2, hlat3):
    for lat in (hlat2, hlat3):
        geom = ET.Element('geometry')
        lat.create_xml_subelement(geom)
        elem = geom.find('hex_lattice')
        assert elem.tag == 'hex_lattice'
        assert elem.get('id') == str(lat.id)
        assert len(elem.find('center').text.split()) == lat.ndim
        assert len(elem.find('pitch').text.split()) == lat.ndim - 1
        assert len(elem.find('universes').text.split()) == len(lat.indices)


def test_show_indices():
    for i in range(1, 11):
        lines = openmc.HexLattice.show_indices(i).split('\n')
        assert len(lines) == 4*i - 3
        lines_x = openmc.HexLattice.show_indices(i, 'x').split('\n')
        assert len(lines) == 4*i - 3
