import xml.etree.ElementTree as ET

import numpy as np
import openmc
import pytest


def test_volume(run_in_tmpdir, uo2):
    """Test adding volume information from a volume calculation."""
    # Create model with nested spheres
    model = openmc.model.Model()
    model.materials.append(uo2)
    inner = openmc.Sphere(r=1.)
    outer = openmc.Sphere(r=2., boundary_type='vacuum')
    c1 = openmc.Cell(fill=uo2, region=-inner)
    c2 = openmc.Cell(region=+inner & -outer)
    u = openmc.Universe(cells=[c1, c2])
    model.geometry.root_universe = u
    model.settings.particles = 100
    model.settings.batches = 10
    model.settings.run_mode = 'fixed source'
    model.settings.source = openmc.Source(space=openmc.stats.Point())

    ll, ur = model.geometry.bounding_box
    assert ll == pytest.approx((-outer.r, -outer.r, -outer.r))
    assert ur == pytest.approx((outer.r, outer.r, outer.r))
    model.settings.volume_calculations

    for domain in (c1, uo2, u):
        # Run stochastic volume calculation
        volume_calc = openmc.VolumeCalculation(
            domains=[domain], samples=1000, lower_left=ll, upper_right=ur)
        model.settings.volume_calculations = [volume_calc]
        model.export_to_xml()
        openmc.calculate_volumes()

        # Load results and add volume information
        volume_calc.load_results('volume_1.h5')
        model.geometry.add_volume_information(volume_calc)

        # get_nuclide_densities relies on volume information
        nucs = set(domain.get_nuclide_densities())
        assert not nucs ^ {'U235', 'O16'}


def test_export_xml(run_in_tmpdir, uo2):
    s1 = openmc.Sphere(r=1.)
    s2 = openmc.Sphere(r=2., boundary_type='reflective')
    c1 = openmc.Cell(fill=uo2, region=-s1)
    c2 = openmc.Cell(fill=uo2, region=+s1 & -s2)
    geom = openmc.Geometry([c1, c2])
    geom.export_to_xml()

    doc = ET.parse('geometry.xml')
    root = doc.getroot()
    assert root.tag == 'geometry'
    cells = root.findall('cell')
    assert [int(c.get('id')) for c in cells] == [c1.id, c2.id]
    surfs = root.findall('surface')
    assert [int(s.get('id')) for s in surfs] == [s1.id, s2.id]


def test_find(uo2):
    xp = openmc.XPlane()
    c1 = openmc.Cell(fill=uo2, region=+xp)
    c2 = openmc.Cell(region=-xp)
    u1 = openmc.Universe(cells=(c1, c2))

    cyl = openmc.ZCylinder()
    c3 = openmc.Cell(fill=u1, region=-cyl)
    c4 = openmc.Cell(region=+cyl)
    geom = openmc.Geometry((c3, c4))

    seq = geom.find((0.5, 0., 0.))
    assert seq[-1] == c1
    seq = geom.find((-0.5, 0., 0.))
    assert seq[-1] == c2
    seq = geom.find((-1.5, 0., 0.))
    assert seq[-1] == c4


def test_get_all_cells():
    cells = [openmc.Cell() for i in range(5)]
    cells2 = [openmc.Cell() for i in range(3)]
    cells[0].fill = openmc.Universe(cells=cells2)
    geom = openmc.Geometry(cells)

    all_cells = set(geom.get_all_cells().values())
    assert not all_cells ^ set(cells + cells2)


def test_get_all_materials():
    m1 = openmc.Material()
    m2 = openmc.Material()
    c1 = openmc.Cell(fill=m1)
    u1 = openmc.Universe(cells=[c1])

    s = openmc.Sphere()
    c2 = openmc.Cell(fill=u1, region=-s)
    c3 = openmc.Cell(fill=m2, region=+s)
    geom = openmc.Geometry([c2, c3])

    all_mats = set(geom.get_all_materials().values())
    assert not all_mats ^ {m1, m2}


def test_get_all_material_cells():
    m1 = openmc.Material()
    m2 = openmc.Material()
    c1 = openmc.Cell(fill=m1)
    u1 = openmc.Universe(cells=[c1])

    s = openmc.Sphere()
    c2 = openmc.Cell(fill=u1, region=-s)
    c3 = openmc.Cell(fill=m2, region=+s)
    geom = openmc.Geometry([c2, c3])

    all_cells = set(geom.get_all_material_cells().values())
    assert not all_cells ^ {c1, c3}


def test_get_all_material_universes():
    m1 = openmc.Material()
    m2 = openmc.Material()
    c1 = openmc.Cell(fill=m1)
    u1 = openmc.Universe(cells=[c1])

    s = openmc.Sphere()
    c2 = openmc.Cell(fill=u1, region=-s)
    c3 = openmc.Cell(fill=m2, region=+s)
    geom = openmc.Geometry([c2, c3])

    all_univs = set(geom.get_all_material_universes().values())
    assert not all_univs ^ {u1, geom.root_universe}


def test_get_all_lattices(cell_with_lattice):
    cells, mats, univ, lattice = cell_with_lattice
    geom = openmc.Geometry([cells[-1]])

    lats = list(geom.get_all_lattices().values())
    assert lats == [lattice]


def test_get_all_surfaces(uo2):
    planes = [openmc.ZPlane(z0=z) for z in np.linspace(-100., 100.)]
    slabs = []
    for region in openmc.model.subdivide(planes):
        slabs.append(openmc.Cell(fill=uo2, region=region))
    geom = openmc.Geometry(slabs)

    surfs = set(geom.get_all_surfaces().values())
    assert not surfs ^ set(planes)


def test_get_by_name():
    m1 = openmc.Material(name='zircaloy')
    m1.add_element('Zr', 1.0)
    m2 = openmc.Material(name='Zirconium')
    m2.add_element('Zr', 1.0)

    c1 = openmc.Cell(fill=m1, name='cell1')
    u1 = openmc.Universe(name='Zircaloy universe', cells=[c1])

    cyl = openmc.ZCylinder()
    c2 = openmc.Cell(fill=u1, region=-cyl, name='cell2')
    c3 = openmc.Cell(fill=m2, region=+cyl, name='Cell3')
    root = openmc.Universe(name='root Universe', cells=[c2, c3])
    geom = openmc.Geometry(root)

    mats = set(geom.get_materials_by_name('zirc'))
    assert not mats ^ {m1, m2}
    mats = set(geom.get_materials_by_name('zirc', True))
    assert not mats ^ {m1}
    mats = set(geom.get_materials_by_name('zirconium', False, True))
    assert not mats ^ {m2}
    mats = geom.get_materials_by_name('zirconium', True, True)
    assert not mats

    cells = set(geom.get_cells_by_name('cell'))
    assert not cells ^ {c1, c2, c3}
    cells = set(geom.get_cells_by_name('cell', True))
    assert not cells ^ {c1, c2}
    cells = set(geom.get_cells_by_name('cell3', False, True))
    assert not cells ^ {c3}
    cells = geom.get_cells_by_name('cell3', True, True)
    assert not cells

    cells = set(geom.get_cells_by_fill_name('Zircaloy'))
    assert not cells ^ {c1, c2}
    cells = set(geom.get_cells_by_fill_name('Zircaloy', True))
    assert not cells ^ {c2}
    cells = set(geom.get_cells_by_fill_name('Zircaloy', False, True))
    assert not cells ^ {c1}
    cells = geom.get_cells_by_fill_name('Zircaloy', True, True)
    assert not cells

    univs = set(geom.get_universes_by_name('universe'))
    assert not univs ^ {u1, root}
    univs = set(geom.get_universes_by_name('universe', True))
    assert not univs ^ {u1}
    univs = set(geom.get_universes_by_name('universe', True, True))
    assert not univs


def test_hex_prism():
    hex_prism = openmc.model.hexagonal_prism(edge_length=5.0,
                                             origin=(0.0, 0.0),
                                             orientation='y')
    # clear checks
    assert (0.0, 0.0, 0.0) in hex_prism
    assert (10.0, 10.0, 10.0) not in hex_prism
    # edge checks
    assert (0.0, 5.01, 0.0) not in hex_prism
    assert (0.0, 4.99, 0.0) in hex_prism

    rounded_hex_prism = openmc.model.hexagonal_prism(edge_length=5.0,
                                                     origin=(0.0, 0.0),
                                                     orientation='y',
                                                     corner_radius=1.0)

    # clear checks
    assert (0.0, 0.0, 0.0) in rounded_hex_prism
    assert (10.0, 10.0, 10.0) not in rounded_hex_prism
    # edge checks
    assert (0.0, 5.01, 0.0) not in rounded_hex_prism
    assert (0.0, 4.99, 0.0) not in rounded_hex_prism


def test_get_lattice_by_name(cell_with_lattice):
    cells, _, _, lattice = cell_with_lattice
    geom = openmc.Geometry([cells[-1]])

    f = geom.get_lattices_by_name
    assert f('lattice') == [lattice]
    assert f('lattice', True) == []
    assert f('Lattice', True) == [lattice]
    assert f('my lattice', False, True) == [lattice]
    assert f('my lattice', True, True) == []


def test_clone():
    c1 = openmc.Cell()
    c2 = openmc.Cell()
    root = openmc.Universe(cells=[c1, c2])
    geom = openmc.Geometry(root)

    clone = geom.clone()
    root_clone = clone.root_universe

    assert root.id != root_clone.id
    assert not (set(root.cells) & set(root_clone.cells))


def test_determine_paths(cell_with_lattice):
    cells, mats, univ, lattice = cell_with_lattice
    u = openmc.Universe(cells=[cells[-1]])
    geom = openmc.Geometry(u)

    geom.determine_paths()
    assert len(cells[0].paths) == 4
    assert len(cells[1].paths) == 4
    assert len(cells[2].paths) == 1
    assert len(mats[0].paths) == 1
    assert len(mats[-1].paths) == 4

    # Test get_instances
    for i in range(4):
        assert geom.get_instances(cells[0].paths[i]) == i
        assert geom.get_instances(mats[-1].paths[i]) == i


def test_from_xml(run_in_tmpdir, mixed_lattice_model):
    # Export model
    mixed_lattice_model.export_to_xml()

    # Import geometry
    geom = openmc.Geometry.from_xml()
    assert isinstance(geom, openmc.Geometry)
    ll, ur = geom.bounding_box
    assert ll == pytest.approx((-6.0, -6.0, -np.inf))
    assert ur == pytest.approx((6.0, 6.0, np.inf))


def test_rotation_matrix():
    """Test ability to set a rotation matrix directly"""
    y = openmc.YPlane()
    cyl1 = openmc.ZCylinder(r=1.0)
    cyl2 = openmc.ZCylinder(r=2.0, boundary_type='vacuum')

    # Create a universe and then reflect in the y-direction
    c1 = openmc.Cell(region=-cyl1 & +y)
    c2 = openmc.Cell(region=+cyl1 & +y)
    c3 = openmc.Cell(region=-y)
    univ = openmc.Universe(cells=[c1, c2, c3])
    c = openmc.Cell(fill=univ, region=-cyl2)
    c.rotation = [[1, 0, 0], [0, -1, 0], [0, 0, 1]]
    assert np.all(c.rotation_matrix == c.rotation)
    geom = openmc.Geometry([c])

    assert geom.find((0.0, 0.5, 0.0))[-1] == c3
    assert geom.find((0.0, -0.5, 0.0))[-1] == c1
    assert geom.find((0.0, -1.5, 0.0))[-1] == c2

def test_remove_redundant_surfaces():
    """Test ability to remove redundant surfaces"""

    m1 = openmc.Material()
    m1.add_nuclide('U235', 1.0, 'wo')
    m1.add_nuclide('O16', 2.0, 'wo')
    m1.set_density('g/cm3', 10.0)

    m2 = openmc.Material()
    m2.add_element('Zr', 1.0)
    m2.set_density('g/cm3', 2.0)

    m3 = openmc.Material()
    m3.add_nuclide('H1', 2.0)
    m3.add_nuclide('O16', 1.0)
    m3.set_density('g/cm3', 1.0)

    def get_cyl_cell(r1, r2, z1, z2, fill):
        """Create a finite height cylindrical cell with a fill"""

        cyl2 = openmc.ZCylinder(r=r2)
        zplane1 = openmc.ZPlane(z1)
        zplane2 = openmc.ZPlane(z2)
        if np.isclose(r1, 0):
            region = -cyl2 & +zplane1 & -zplane2
        else:
            cyl1 = openmc.ZCylinder(r=r1)
            region = +cyl1 & -cyl2 & +zplane1 & -zplane2
        return openmc.Cell(region=region, fill=fill)

    r1, r2, r3 = 1., 2., 3.
    z1, z2 = -2., 2.
    fuel = get_cyl_cell(0, r1, z1, z2, m1)
    clad = get_cyl_cell(r1, r2, z1, z2, m2)
    water = get_cyl_cell(r2, r3, z1, z2, m3)
    root = openmc.Universe(cells=[fuel, clad, water])
    geom = openmc.Geometry(root)
    
    # There should be 6 redundant surfaces in this geometry
    n_redundant_surfs = len(geom.get_redundant_surfaces().keys())
    assert n_redundant_surfs == 6
    # Remove redundant surfaces
    geom.remove_redundant_surfaces()
    # There should be 0 remaining redundant surfaces
    n_redundant_surfs = len(geom.get_redundant_surfaces().keys())
    assert n_redundant_surfs == 0
