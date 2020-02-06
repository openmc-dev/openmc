import xml.etree.ElementTree as ET

import numpy as np
import openmc
import pytest

from tests.unit_tests import assert_unbounded


def test_basic():
    c1 = openmc.Cell()
    c2 = openmc.Cell()
    c3 = openmc.Cell()
    u = openmc.Universe(name='cool', cells=(c1, c2, c3))
    assert u.name == 'cool'

    cells = set(u.cells.values())
    assert not (cells ^ {c1, c2, c3})

    # Test __repr__
    repr(u)

    with pytest.raises(TypeError):
        u.add_cell(openmc.Material())
    with pytest.raises(TypeError):
        u.add_cells(c1)

    u.remove_cell(c3)
    cells = set(u.cells.values())
    assert not (cells ^ {c1, c2})

    u.clear_cells()
    assert not set(u.cells)


def test_bounding_box():
    cyl1 = openmc.ZCylinder(r=1.0)
    cyl2 = openmc.ZCylinder(r=2.0)
    c1 = openmc.Cell(region=-cyl1)
    c2 = openmc.Cell(region=+cyl1 & -cyl2)

    u = openmc.Universe(cells=[c1, c2])
    ll, ur = u.bounding_box
    assert ll == pytest.approx((-2., -2., -np.inf))
    assert ur == pytest.approx((2., 2., np.inf))

    u = openmc.Universe()
    assert_unbounded(u)


def test_plot(run_in_tmpdir, sphere_model):
    m = sphere_model.materials[0]
    univ = sphere_model.geometry.root_universe

    colors = {m: 'limegreen'}
    for basis in ('xy', 'yz', 'xz'):
        univ.plot(
            basis=basis,
            pixels=(10, 10),
            color_by='material',
            colors=colors,
        )


def test_get_nuclides(uo2):
    c = openmc.Cell(fill=uo2)
    univ = openmc.Universe(cells=[c])
    nucs = univ.get_nuclides()
    assert nucs == ['U235', 'O16']


def test_cells():
    cells = [openmc.Cell() for i in range(5)]
    cells2 = [openmc.Cell() for i in range(3)]
    cells[0].fill = openmc.Universe(cells=cells2)
    u = openmc.Universe(cells=cells)
    assert not (set(u.cells.values()) ^ set(cells))

    all_cells = set(u.get_all_cells().values())
    assert not (all_cells ^ set(cells + cells2))


def test_get_all_materials(cell_with_lattice):
    cells, mats, univ, lattice = cell_with_lattice
    test_mats = set(univ.get_all_materials().values())
    assert not (test_mats ^ set(mats))


def test_get_all_universes():
    c1 = openmc.Cell()
    u1 = openmc.Universe(cells=[c1])
    c2 = openmc.Cell()
    u2 = openmc.Universe(cells=[c2])
    c3 = openmc.Cell(fill=u1)
    c4 = openmc.Cell(fill=u2)
    u3 = openmc.Universe(cells=[c3, c4])

    univs = set(u3.get_all_universes().values())
    assert not (univs ^ {u1, u2})


def test_clone():

    c1 = openmc.Cell(cell_id=1)
    c1.region = -openmc.ZCylinder(r=1.0)
    c2 = openmc.Cell(cell_id=2)
    c2.fill = openmc.Material()
    c3 = openmc.Cell()
    u1 = openmc.Universe(name='cool', cells=(c1, c2, c3))

    u2 = u1.clone()
    assert u2.name == u1.name
    assert u2.cells != u1.cells
    assert u2.get_all_materials() != u1.get_all_materials()

    u2 = u1.clone(clone_materials=False)
    assert u2.get_all_materials() == u1.get_all_materials()

    u3 = u1.clone(clone_regions=False)
    assert next(iter(u3.cells.values())).region ==\
         next(iter(u1.cells.values())).region


def test_create_xml(cell_with_lattice):
    cells = [openmc.Cell() for i in range(5)]
    u = openmc.Universe(cells=cells)

    geom = ET.Element('geom')
    u.create_xml_subelement(geom)
    cell_elems = geom.findall('cell')
    assert len(cell_elems) == len(cells)
    assert all(c.get('universe') == str(u.id) for c in cell_elems)
    assert not (set(c.get('id') for c in cell_elems) ^
                set(str(c.id) for c in cells))
