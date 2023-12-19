import pytest
from unittest import TestCase
import numpy as np

import openmc
import openmc.lib

from tests import cdtemp

<<<<<<< HEAD

@pytest.fixture(scope='module', autouse=True)
def double_lattice_model():
=======
@pytest.fixture(scope='class')
def double_rect_lattice_model():
>>>>>>> b5b6beafe (Adding tests for hex lattice)
    openmc.reset_auto_ids()
    model = openmc.Model()

    # Create a single material
    m = openmc.Material()
    m.add_nuclide('U235', 1.0)
    m.set_density('g/cm3', 10.0)
    model.materials.append(m)

    # Create a universe with a single infinite cell
    c = openmc.Cell(fill=m)
    u = openmc.Universe(cells=[c])

    # Create a 2x2 lattice filled with above universe
    lattice = openmc.RectLattice()
    lattice.lower_left = (0.0, 0.0)
    lattice.pitch = (1.0, 1.0)
    lattice.universes = np.full((2, 2), u)

    # Create two cells each filled with the same lattice, one from x=0..2 and
    # y=0..2 and the other from x=2..4 and y=0..2
    x0 = openmc.XPlane(0.0, boundary_type='vacuum')
    x2 = openmc.XPlane(2.0)
    x4 = openmc.XPlane(4.0, boundary_type='vacuum')
    y0 = openmc.YPlane(0.0, boundary_type='vacuum')
    y2 = openmc.YPlane(2.0, boundary_type='vacuum')
    cell_with_lattice1 = openmc.Cell(fill=lattice, region=+x0 & -x2 & +y0 & -y2)
    cell_with_lattice2 = openmc.Cell(fill=lattice, region=+x2 & -x4 & +y0 & -y2)
    cell_with_lattice2.translation = (2., 0., 0.)
    model.geometry = openmc.Geometry([cell_with_lattice1, cell_with_lattice2])

    tally = openmc.Tally()
    tally.filters = [openmc.DistribcellFilter(c)]
    tally.scores = ['flux']
    model.tallies = [tally]

    # Add box source that covers the model space well
    bbox = model.geometry.bounding_box
    bbox[0][2] = -0.5
    bbox[1][2] = 0.5
    space = openmc.stats.Box(*bbox)
    model.settings.source = openmc.IndependentSource(space=space)

    # Add necessary settings and export
    model.settings.batches = 10
    model.settings.inactive = 0
    model.settings.particles = 100

    with cdtemp():
        model.export_to_xml()
        openmc.lib.init()
        yield
        openmc.lib.finalize()

# This shows the expected cell instance numbers for each lattice position:
#      ┌─┬─┬─┬─┐
#      │2│3│6│7│
#      ├─┼─┼─┼─┤
#      │0│1│4│5│
#      └─┴─┴─┴─┘
rect_expected_results = [
    ((0.5, 0.5, 0.0), 0),
    ((1.5, 0.5, 0.0), 1),
    ((0.5, 1.5, 0.0), 2),
    ((1.5, 1.5, 0.0), 3),
    ((2.5, 0.5, 0.0), 4),
    ((3.5, 0.5, 0.0), 5),
    ((2.5, 1.5, 0.0), 6),
    ((3.5, 1.5, 0.0), 7),
]


@pytest.mark.parametrize("r,expected_cell_instance", rect_expected_results, ids=lambda p : f'{p}')
def test_cell_instance_rect_multilattice(double_rect_lattice_model, r, expected_cell_instance):
    _, cell_instance = openmc.lib.find_cell(r)
    assert cell_instance == expected_cell_instance


<<<<<<< HEAD
def test_cell_instance_multilattice_results():
    with openmc.lib.run_in_memory():
        openmc.lib.run()
        tally_results = openmc.lib.tallies[1].mean
        assert (tally_results != 0.0).all()
=======
def test_cell_instance_multilattice_results(double_rect_lattice_model):
    openmc.run()
    openmc.lib.run()
    tally_results = openmc.lib.tallies[1].mean
    assert (tally_results != 0.0).all()


@pytest.fixture(scope='function')
def double_hex_lattice_model():
    openmc.reset_auto_ids()
    radius = 0.9
    pin_lattice_pitch = 2.0
    # make the hex prism a little larger to make sure test
    # locations are definitively in the model
    hex_prism_edge = 1.2 * pin_lattice_pitch

    model = openmc.Model()

    # materials
    nat_u = openmc.Material()
    nat_u.set_density('g/cm3', 12.0)
    nat_u.add_element('U', 1.0)

    graphite = openmc.Material()
    graphite.set_density('g/cm3', 1.1995)
    graphite.add_element('C', 1.0)

    # zplanes to define lower and upper region
    z_low = openmc.ZPlane(z0=-10, boundary_type='vacuum')
    z_mid = openmc.ZPlane(z0=0)
    z_high = openmc.ZPlane(z0=10, boundary_type='vacuum')
    hex_prism = openmc.model.HexagonalPrism(edge_length=1.2*pin_lattice_pitch, boundary_type='reflective')

    # geometry
    cyl = openmc.ZCylinder(r=radius)
    univ = openmc.model.pin([cyl], [nat_u, graphite])

    # create a hexagonal lattice of compacts
    hex_lattice = openmc.HexLattice()
    hex_lattice.orientation = 'y'
    hex_lattice.pitch = (pin_lattice_pitch,)
    hex_lattice.center = (0., 0.)
    center = [univ]
    ring = [univ, univ, univ, univ, univ, univ]
    hex_lattice.universes = [ring, center]
    lower_hex_cell = openmc.Cell(fill=hex_lattice, region=-hex_prism & +z_low & -z_mid)
    upper_hex_cell = openmc.Cell(fill=hex_lattice, region=-hex_prism & +z_mid & -z_high)
    hex_cells = [lower_hex_cell, upper_hex_cell]
    model.geometry = openmc.Geometry(hex_cells)

    # moderator
    cell = next(iter(univ.get_all_cells().values()))
    tally = openmc.Tally(tally_id=1)
    filter = openmc.DistribcellFilter(cell)
    tally.filters = [filter]
    tally.scores = ['flux']
    model.tallies = [tally]

    # settings
    settings = openmc.Settings()
    # source definition. fission source given bounding box of graphite active region
    system_LL = (-pin_lattice_pitch*np.sqrt(3)/2, -pin_lattice_pitch, -5)
    system_UR = (pin_lattice_pitch*np.sqrt(3)/2, pin_lattice_pitch, 5)
    source_dist = openmc.stats.Box(system_LL, system_UR)
    source = openmc.IndependentSource(space=source_dist)
    settings.source = source
    settings.particles = 100
    settings.inactive = 2
    settings.batches = 10
    model.settings = settings

    model.export_to_xml()
    with cdtemp():
        model.export_to_xml()
        openmc.lib.init()
        yield
        openmc.lib.finalize()

# Lower cell instances
#          6
#     5         4
#          3
#     2         1
#          0
# Upper cell instances
#          13
#     12        11
#          10
#     9         8
#          7

hex_expected_results = [
    ((0.0, -2.0, -5.0), 0),
    ((1.732, -1.0, -5.0), 1),
    ((-1.732, -1.0, -5.0), 2),
    ((0.0, 0.0, -0.1), 3),
    ((1.732, 1.0, -5.0), 4),
    ((-1.732, 1.0, -5.0), 5),
    ((0.0, 2.0, -0.1), 6),
    ((0.0, -2.0, 5.0), 7),
    ((1.732, -1.0, 5.0), 8),
    ((-1.732, -1.0, 5.0), 9),
    ((0.0, 0.0, 5.0), 10),
    ((1.732, 1.0, 5.0), 11),
    ((-1.732, 1.0, 5.0), 12),
    ((0.0, 2.0, 5.0), 13),
]
@pytest.mark.parametrize("r,expected_cell_instance", hex_expected_results, ids=lambda p : f'{p}')
def test_cell_instance_hex_multilattice(double_hex_lattice_model, r, expected_cell_instance):
    print(r)
    _, cell_instance = openmc.lib.find_cell(r)
    assert cell_instance == expected_cell_instance


def test_cell_instance_hex_multilattice_results(double_hex_lattice_model):
    openmc.run()
    openmc.lib.run()
    tally_results = openmc.lib.tallies[1].mean
    assert (tally_results != 0.0).all()
>>>>>>> b5b6beafe (Adding tests for hex lattice)
