from math import sqrt

import pytest
import numpy as np
import openmc
import openmc.lib

from tests import cdtemp


@pytest.fixture(scope="module", autouse=True)
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
    nat_u.set_density("g/cm3", 12.0)
    nat_u.add_element("U", 1.0)

    graphite = openmc.Material()
    graphite.set_density("g/cm3", 1.1995)
    graphite.add_element("C", 1.0)

    # zplanes to define lower and upper region
    z_low = openmc.ZPlane(-10, boundary_type="vacuum")
    z_mid = openmc.ZPlane(0)
    z_high = openmc.ZPlane(10, boundary_type="vacuum")
    hex_prism = openmc.model.HexagonalPrism(
        edge_length=hex_prism_edge, boundary_type="reflective"
    )

    # geometry
    cyl = openmc.ZCylinder(r=radius)
    univ = openmc.model.pin([cyl], [nat_u, graphite])

    # create a hexagonal lattice of compacts
    hex_lattice = openmc.HexLattice()
    hex_lattice.orientation = "y"
    hex_lattice.pitch = (pin_lattice_pitch,)
    hex_lattice.center = (0.0, 0.0)
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
    tally.scores = ["flux"]
    model.tallies = [tally]

    # settings
    # source definition. fission source given bounding box of graphite active region
    system_LL = (-pin_lattice_pitch * sqrt(3) / 2, -pin_lattice_pitch, -5)
    system_UR = (pin_lattice_pitch * sqrt(3) / 2, pin_lattice_pitch, 5)
    source_dist = openmc.stats.Box(system_LL, system_UR)
    model.settings.source = openmc.IndependentSource(space=source_dist)
    model.settings.particles = 100
    model.settings.inactive = 2
    model.settings.batches = 10

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


@pytest.mark.parametrize("r,expected_cell_instance", hex_expected_results, ids=str)
def test_cell_instance_hex_multilattice(r, expected_cell_instance):
    _, cell_instance = openmc.lib.find_cell(r)
    assert cell_instance == expected_cell_instance


def test_cell_instance_multilattice_results():
    openmc.lib.run()
    tally_results = openmc.lib.tallies[1].mean
    assert (tally_results != 0.0).all()
