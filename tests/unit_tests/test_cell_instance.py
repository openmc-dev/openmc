import numpy as np
import pytest

import openmc
import openmc.lib

from tests import cdtemp


@pytest.fixture(scope='module', autouse=True)
def double_lattice_model():
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
expected_results = [
    ((0.5, 0.5, 0.0), 0),
    ((1.5, 0.5, 0.0), 1),
    ((0.5, 1.5, 0.0), 2),
    ((1.5, 1.5, 0.0), 3),
    ((2.5, 0.5, 0.0), 4),
    ((3.5, 0.5, 0.0), 5),
    ((2.5, 1.5, 0.0), 6),
    ((3.5, 1.5, 0.0), 7),
]
@pytest.mark.parametrize("r,expected_cell_instance", expected_results)
def test_cell_instance_multilattice(r, expected_cell_instance):
    _, cell_instance = openmc.lib.find_cell(r)
    assert cell_instance == expected_cell_instance


def test_cell_instance_multilattice_results():
    openmc.lib.run()
    tally_results = openmc.lib.tallies[1].mean
    assert (tally_results != 0.0).all()
