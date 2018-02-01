import xml.etree. ElementTree as ET

import numpy as np
import openmc
import pytest

from tests.unit_tests import assert_unbounded


def test_contains():
    # Cell with specified region
    s = openmc.XPlane()
    c = openmc.Cell(region=+s)
    assert (1.0, 0.0, 0.0) in c
    assert (-1.0, 0.0, 0.0) not in c

    # Cell with no region
    c = openmc.Cell()
    assert (10.0, -4., 2.0) in c


def test_repr(cell_with_lattice):
    cells, mats, univ, lattice = cell_with_lattice
    repr(cells[0])  # cell with distributed materials
    repr(cells[1])  # cell with material
    repr(cells[2])  # cell with lattice

    # Empty cell
    c = openmc.Cell()
    repr(c)


def test_bounding_box():
    zcyl = openmc.ZCylinder()
    c = openmc.Cell(region=-zcyl)
    ll, ur = c.bounding_box
    assert ll == pytest.approx((-1., -1., -np.inf))
    assert ur == pytest.approx((1., 1., np.inf))

    # Cell with no region specified
    c = openmc.Cell()
    assert_unbounded(c)


def test_clone():
    m = openmc.Material()
    cyl = openmc.ZCylinder()
    c = openmc.Cell(fill=m, region=-cyl)
    c.temperature = 650.

    c2 = c.clone()
    assert c2.id != c.id
    assert c2.fill != c.fill
    assert c2.region != c.region
    assert c2.temperature == c.temperature


def test_temperature(cell_with_lattice):
    # Make sure temperature propagates through universes
    m = openmc.Material()
    s = openmc.XPlane()
    c1 = openmc.Cell(fill=m, region=+s)
    c2 = openmc.Cell(fill=m, region=-s)
    u1 = openmc.Universe(cells=[c1, c2])
    c = openmc.Cell(fill=u1)

    c.temperature = 400.0
    assert c1.temperature == 400.0
    assert c2.temperature == 400.0
    with pytest.raises(ValueError):
        c.temperature = -100.

    # distributed temperature
    cells, _, _, _ = cell_with_lattice
    c = cells[0]
    c.temperature = (300., 600., 900.)


def test_rotation():
    u = openmc.Universe()
    c = openmc.Cell(fill=u)
    c.rotation = (180.0, 0.0, 0.0)
    assert np.allclose(c.rotation_matrix, [
        [1., 0., 0.],
        [0., -1., 0.],
        [0., 0., -1.]
    ])

    c.rotation = (0.0, 90.0, 0.0)
    assert np.allclose(c.rotation_matrix, [
        [0., 0., -1.],
        [0., 1., 0.],
        [1., 0., 0.]
    ])


def test_get_nuclides(uo2):
    c = openmc.Cell(fill=uo2)
    nucs = c.get_nuclides()
    assert nucs == ['U235', 'O16']


def test_nuclide_densities(uo2):
    c = openmc.Cell(fill=uo2)
    expected_nucs = ['U235', 'O16']
    expected_density = [1.0, 2.0]
    tuples = list(c.get_nuclide_densities().values())
    for nuc, density, t in zip(expected_nucs, expected_density, tuples):
        assert nuc == t[0]
        assert density == t[1]

    # Empty cell
    c = openmc.Cell()
    assert not c.get_nuclide_densities()


def test_get_all_universes(cell_with_lattice):
    # Cell with nested universes
    c1 = openmc.Cell()
    u1 = openmc.Universe(cells=[c1])
    c2 = openmc.Cell(fill=u1)
    u2 = openmc.Universe(cells=[c2])
    c3 = openmc.Cell(fill=u2)
    univs = set(c3.get_all_universes().values())
    assert not (univs ^ {u1, u2})

    # Cell with lattice
    cells, mats, univ, lattice = cell_with_lattice
    univs = set(cells[-1].get_all_universes().values())
    assert not (univs ^ {univ})


def test_get_all_materials(cell_with_lattice):
    # Normal cell
    m = openmc.Material()
    c = openmc.Cell(fill=m)
    test_mats = set(c.get_all_materials().values())
    assert not(test_mats ^ {m})

    # Cell filled with distributed materials
    cells, mats, univ, lattice = cell_with_lattice
    c = cells[0]
    test_mats = set(c.get_all_materials().values())
    assert not (test_mats ^ set(m for m in c.fill if m is not None))

    # Cell filled with universe
    c = cells[-1]
    test_mats = set(c.get_all_materials().values())
    assert not (test_mats ^ set(mats))


def test_to_xml_element(cell_with_lattice):
    cells, mats, univ, lattice = cell_with_lattice

    c = cells[-1]
    root = ET.Element('geometry')
    elem = c.create_xml_subelement(root)
    assert elem.tag == 'cell'
    assert elem.get('id') == str(c.id)
    assert elem.get('region') is None
    surf_elem = root.find('surface')
    assert surf_elem.get('id') == str(cells[0].region.surface.id)

    c = cells[0]
    c.temperature = 900.0
    elem = c.create_xml_subelement(root)
    assert elem.get('region') == str(c.region)
    assert elem.get('temperature') == str(c.temperature)
