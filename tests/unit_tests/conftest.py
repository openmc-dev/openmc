import openmc
import pytest

from tests.regression_tests import config


@pytest.fixture(scope='module')
def mpi_intracomm():
    if config['mpi']:
        from mpi4py import MPI
        return MPI.COMM_WORLD
    else:
        return None


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

    cyl = openmc.ZCylinder(r=1.0)
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

@pytest.fixture
def mixed_lattice_model(uo2, water):
    cyl = openmc.ZCylinder(r=0.4)
    c1 = openmc.Cell(fill=uo2, region=-cyl)
    c1.temperature = 600.0
    c2 = openmc.Cell(fill=water, region=+cyl)
    pin = openmc.Universe(cells=[c1, c2])

    empty = openmc.Cell()
    empty_univ = openmc.Universe(cells=[empty])

    hex_lattice = openmc.HexLattice()
    hex_lattice.center = (0.0, 0.0)
    hex_lattice.pitch = (1.2, 10.0)
    outer_ring = [pin]*6
    inner_ring = [empty_univ]
    axial_level = [outer_ring, inner_ring]
    hex_lattice.universes = [axial_level]*3
    hex_lattice.outer = empty_univ

    cell_hex = openmc.Cell(fill=hex_lattice)
    u = openmc.Universe(cells=[cell_hex])
    rotated_cell_hex = openmc.Cell(fill=u)
    rotated_cell_hex.rotation = (0., 0., 30.)
    ur = openmc.Universe(cells=[rotated_cell_hex])

    d = 6.0
    rect_lattice = openmc.RectLattice()
    rect_lattice.lower_left = (-d, -d)
    rect_lattice.pitch = (d, d)
    rect_lattice.outer = empty_univ
    rect_lattice.universes = [
        [ur, empty_univ],
        [empty_univ, u]
    ]

    xmin = openmc.XPlane(-d, boundary_type='periodic')
    xmax = openmc.XPlane(d, boundary_type='periodic')
    xmin.periodic_surface = xmax
    ymin = openmc.YPlane(-d, boundary_type='periodic')
    ymax = openmc.YPlane(d, boundary_type='periodic')
    main_cell = openmc.Cell(fill=rect_lattice,
                            region=+xmin & -xmax & +ymin & -ymax)

    # Create geometry and use unique material in each fuel cell
    geometry = openmc.Geometry([main_cell])
    geometry.determine_paths()
    c1.fill = [water.clone() for i in range(c1.num_instances)]

    return openmc.model.Model(geometry)
