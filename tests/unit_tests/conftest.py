import openmc
import pytest


@pytest.fixture(scope='module')
def uo2():
    m = openmc.Material(material_id=100, name='UO2')
    m.add_nuclide('U235', 1.0)
    m.add_nuclide('O16', 2.0)
    m.set_density('g/cm3', 10.0)
    m.depletable = True
    return m


@pytest.fixture(scope='module')
def water():
    m = openmc.Material(name='light water')
    m.add_nuclide('H1', 2.0)
    m.add_nuclide('O16', 1.0)
    m.set_density('g/cm3', 1.0)
    m.add_s_alpha_beta('c_H_in_H2O')
    return m


@pytest.fixture(scope='module')
def sphere_model():
    model = openmc.model.Model()
    m = openmc.Material()
    m.add_nuclide('U235', 1.0)
    m.set_density('g/cm3', 1.0)
    model.materials.append(m)

    sph = openmc.Sphere(boundary_type='vacuum')
    c = openmc.Cell(fill=m, region=-sph)
    model.geometry.root_universe = openmc.Universe(cells=[c])

    model.settings.particles = 100
    model.settings.batches = 10
    model.settings.run_mode = 'fixed source'
    model.settings.source = openmc.Source(space=openmc.stats.Point())
    return model


@pytest.fixture
def cell_with_lattice():
    m_inside = [openmc.Material(), openmc.Material(), None, openmc.Material()]
    m_outside = openmc.Material()

    cyl = openmc.ZCylinder(R=1.0)
    inside_cyl = openmc.Cell(fill=m_inside, region=-cyl)
    outside_cyl = openmc.Cell(fill=m_outside, region=+cyl)
    univ = openmc.Universe(cells=[inside_cyl, outside_cyl])

    lattice = openmc.RectLattice(name='My Lattice')
    lattice.lower_left = (-4.0, -4.0)
    lattice.pitch = (4.0, 4.0)
    lattice.universes = [[univ, univ], [univ, univ]]
    main_cell = openmc.Cell(fill=lattice)

    return ([inside_cyl, outside_cyl, main_cell],
            [m_inside[0], m_inside[1], m_inside[3], m_outside],
            univ, lattice)
