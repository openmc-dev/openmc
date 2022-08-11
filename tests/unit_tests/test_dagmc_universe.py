import openmc
import pytest


def test_bounding_box():

    u = openmc.DAGMCUniverse("tests/regression_tests/dagmc/universes/dagmc.h5m")

    ll, ur = u.bounding_box
    assert ll == pytest.approx((-25., -25., -25))
    assert ur == pytest.approx((25., 25., 25))


def test_bounding_region():

    u = openmc.DAGMCUniverse("tests/regression_tests/dagmc/universes/dagmc.h5m")

    region = u.bounding_region()  # should default to bounded_type='box'
    assert isinstance(region, openmc.Region)
    assert len(region) == 6
    assert region[0].surface.type == 'x-plane'
    assert region[1].surface.type == 'x-plane'
    assert region[2].surface.type == 'y-plane'
    assert region[3].surface.type == 'y-plane'
    assert region[4].surface.type == 'z-plane'
    assert region[5].surface.type == 'z-plane'
    assert region[0].surface.boundary_type == 'vacuum'
    assert region[1].surface.boundary_type == 'vacuum'
    assert region[2].surface.boundary_type == 'vacuum'
    assert region[3].surface.boundary_type == 'vacuum'
    assert region[4].surface.boundary_type == 'vacuum'
    assert region[5].surface.boundary_type == 'vacuum'

    region = u.bounding_region(bounded_type='sphere', boundary_type='reflective')
    assert isinstance(region, openmc.Region)
    assert isinstance(region, openmc.Halfspace)
    assert region.surface.type == 'sphere'
    assert region.surface.boundary_type == 'reflective'

def test_bounded_universe():

    u = openmc.DAGMCUniverse("tests/regression_tests/dagmc/universes/dagmc.h5m")

    # bounded with defaults
    bu = u.bounded_universe()

    cells = list(bu.get_all_cells().items())
    assert len(cells) == 1
    assert cells[0][0] == 10000  # default bounding_cell_id is 10000
    assert cells[0][1].id == 10000  # default bounding_cell_id is 10000
    surfaces = list(cells[0][1].region.get_surfaces().items())
    assert len(surfaces) == 6
    assert surfaces[0][1].id == 10000

    # bounded with non defaults
    bu = u.bounded_universe(
        bounding_cell_id=42,
        bounded_type='sphere',
        starting_id=43
    )

    cells = list(bu.get_all_cells().items())
    assert len(cells) == 1
    assert cells[0][0] == 42  # default bounding_cell_id is 10000
    assert cells[0][1].id == 42  # default bounding_cell_id is 10000
    surfaces = list(cells[0][1].region.get_surfaces().items())
    assert surfaces[0][1].type == 'sphere'
    assert surfaces[0][1].id == 43
