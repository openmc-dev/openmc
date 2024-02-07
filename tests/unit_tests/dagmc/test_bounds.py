import openmc
import pytest
from pathlib import Path


def test_bounding_box(request):
    """Checks that the DAGMCUniverse.bounding_box returns the correct values"""

    u = openmc.DAGMCUniverse(Path(request.fspath).parent / "dagmc.h5m")

    ll, ur = u.bounding_box
    assert ll == pytest.approx((-25.0, -25.0, -25))
    assert ur == pytest.approx((25.0, 25.0, 25))


def test_bounding_region(request):
    """Checks that the DAGMCUniverse.bounding_region() returns a region with
    correct surfaces and boundary types"""

    u = openmc.DAGMCUniverse(Path(request.fspath).parent / "dagmc.h5m")

    region = u.bounding_region()  # should default to bounded_type='box'
    assert isinstance(region, openmc.Region)
    assert len(region) == 6
    assert region[0].surface.type == "x-plane"
    assert region[0].surface.x0 == -25.
    assert region[1].surface.type == "x-plane"
    assert region[1].surface.x0 == 25.
    assert region[2].surface.type == "y-plane"
    assert region[2].surface.y0 == -25.
    assert region[3].surface.type == "y-plane"
    assert region[3].surface.y0 == 25.
    assert region[4].surface.type == "z-plane"
    assert region[4].surface.z0 == -25.
    assert region[5].surface.type == "z-plane"
    assert region[5].surface.z0 == 25.
    assert region[0].surface.boundary_type == "vacuum"
    assert region[1].surface.boundary_type == "vacuum"
    assert region[2].surface.boundary_type == "vacuum"
    assert region[3].surface.boundary_type == "vacuum"
    assert region[4].surface.boundary_type == "vacuum"
    assert region[5].surface.boundary_type == "vacuum"
    region = u.bounding_region(padding_distance=5)
    assert region[0].surface.x0 == -30.
    assert region[1].surface.x0 == 30.
    assert region[2].surface.y0 == -30.
    assert region[3].surface.y0 == 30.
    assert region[4].surface.z0 == -30.
    assert region[5].surface.z0 == 30.

    region = u.bounding_region(bounded_type="sphere", boundary_type="reflective")
    assert isinstance(region, openmc.Region)
    assert isinstance(region, openmc.Halfspace)
    assert region.surface.type == "sphere"
    assert region.surface.boundary_type == "reflective"
    larger_region = u.bounding_region(bounded_type="sphere", padding_distance=10)
    assert larger_region.surface.r > region.surface.r


def test_bounded_universe(request):
    """Checks that the DAGMCUniverse.bounded_universe() returns a
    openmc.Universe with correct surface ids and cell ids"""

    u = openmc.DAGMCUniverse(Path(request.fspath).parent / "dagmc.h5m")

    # bounded with defaults
    bu = u.bounded_universe()

    cells = list(bu.get_all_cells().items())
    assert isinstance(bu, openmc.Universe)
    assert len(cells) == 1
    assert cells[0][0] == 10000  # default bounding_cell_id is 10000
    assert cells[0][1].id == 10000  # default bounding_cell_id is 10000
    surfaces = list(cells[0][1].region.get_surfaces().items())
    assert len(surfaces) == 6
    assert surfaces[0][1].id == 10000

    # bounded with non defaults
    bu = u.bounded_universe(bounding_cell_id=42, bounded_type="sphere", starting_id=43)

    cells = list(bu.get_all_cells().items())
    assert isinstance(bu, openmc.Universe)
    assert len(cells) == 1
    assert cells[0][0] == 42  # default bounding_cell_id is 10000
    assert cells[0][1].id == 42  # default bounding_cell_id is 10000
    surfaces = list(cells[0][1].region.get_surfaces().items())
    assert surfaces[0][1].type == "sphere"
    assert surfaces[0][1].id == 43


def test_material_names(request):
    """Checks that the DAGMCUniverse.material_names() returns a list of the
    name present in the dagmc.h5m file in the expected order"""

    u = openmc.DAGMCUniverse(Path(request.fspath).parent / "dagmc.h5m")

    assert u.material_names == ['41', 'Graveyard', 'no-void fuel']
