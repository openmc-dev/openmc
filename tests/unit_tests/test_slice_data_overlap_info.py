"""Tests for GUI-facing slice data overlap diagnostics.

These tests exercise the overlap data that the plotter can show in the GUI
after a user clicks a plotted overlap pixel. The GUI itself does not need to be
opened here; we test the Python/C API that supplies the GUI message.
"""

import numpy as np
import pytest

import openmc


@pytest.fixture
def overlapping_slice_model():
    """Model with one known overlapping region and one known clean region.

    Cells 1 and 2 overlap in the left half of a sphere. Cell 3 occupies the
    right half without overlap. All cells live in root universe 10.
    """
    openmc.reset_auto_ids()

    material = openmc.Material(name="water")
    material.add_nuclide("H1", 2.0)
    material.add_nuclide("O16", 1.0)
    material.set_density("g/cm3", 1.0)

    sphere = openmc.Sphere(r=1.0, boundary_type="vacuum")
    x_plane = openmc.XPlane(x0=0.0)

    overlap_region = -sphere & -x_plane
    clean_region = -sphere & +x_plane

    cell_1 = openmc.Cell(cell_id=1, region=overlap_region, fill=material)
    cell_2 = openmc.Cell(cell_id=2, region=overlap_region, fill=material)
    cell_3 = openmc.Cell(cell_id=3, region=clean_region, fill=material)
    root = openmc.Universe(universe_id=10, cells=[cell_1, cell_2, cell_3])

    model = openmc.Model()
    model.geometry = openmc.Geometry(root)
    model.materials = openmc.Materials([material])
    return model


@pytest.fixture
def triple_overlapping_slice_model():
    """Model with three cells overlapping in one region.

    Cells 1, 2, and 3 all occupy the left half of the sphere. Cell 4 occupies
    the right half without overlap. All cells live in root universe 10.
    """
    openmc.reset_auto_ids()

    material = openmc.Material(name="water")
    material.add_nuclide("H1", 2.0)
    material.add_nuclide("O16", 1.0)
    material.set_density("g/cm3", 1.0)

    sphere = openmc.Sphere(r=1.0, boundary_type="vacuum")
    x_plane = openmc.XPlane(x0=0.0)

    overlap_region = -sphere & -x_plane
    clean_region = -sphere & +x_plane

    cell_1 = openmc.Cell(cell_id=1, region=overlap_region, fill=material)
    cell_2 = openmc.Cell(cell_id=2, region=overlap_region, fill=material)
    cell_3 = openmc.Cell(cell_id=3, region=overlap_region, fill=material)
    cell_4 = openmc.Cell(cell_id=4, region=clean_region, fill=material)
    root = openmc.Universe(universe_id=10, cells=[cell_1, cell_2, cell_3, cell_4])

    model = openmc.Model()
    model.geometry = openmc.Geometry(root)
    model.materials = openmc.Materials([material])
    return model


def _run_overlap_slice(model):
    """Run a small slice plot and return the geometry data while initialized."""
    import openmc.lib

    geom_data, _ = openmc.lib.slice_data(
        origin=(0.0, 0.0, 0.0),
        width=(2.0, 2.0),
        pixels=(20, 20),
        basis="xy",
        show_overlaps=True,
        include_properties=False,
    )
    return geom_data


def _first_pixel_with_cell_id(geom_data, cell_id):
    """Return the first (x, y) pixel whose plotted cell field equals cell_id."""
    matches = np.argwhere(geom_data[:, :, 0] == cell_id)
    assert matches.size > 0, f"No pixel found for cell ID {cell_id}"
    y, x = matches[0]
    return int(x), int(y)


def _first_overlap_pixel(geom_data):
    """Return the first (x, y) pixel marked as an overlap."""
    matches = np.argwhere(geom_data[:, :, 0] == -3)
    assert matches.size > 0, "No overlap pixel was found"
    y, x = matches[0]
    return int(x), int(y)


def _overlap_records(x, y):
    """Return overlap data for one pixel as (cell1, cell2, universe) tuples."""
    import openmc.lib

    cell1, cell2, universe = openmc.lib.slice_data_overlap_info(x, y)
    return list(zip(cell1.tolist(), cell2.tolist(), universe.tolist()))


def test_slice_data_overlap_info_reports_cells_and_universe(
    run_in_tmpdir, overlapping_slice_model
):
    """An overlap pixel should report the exact cells and universe to the GUI."""
    import openmc.lib

    with openmc.lib.TemporarySession(
        overlapping_slice_model, output=False, args=["-c"]
    ):
        geom_data = _run_overlap_slice(overlapping_slice_model)
        x, y = _first_overlap_pixel(geom_data)

        assert openmc.lib.slice_data_overlap_count(x, y) == 1

        cell1, cell2, universe = openmc.lib.slice_data_overlap_info(x, y)

        assert cell1.tolist() == [1]
        assert cell2.tolist() == [2]
        assert universe.tolist() == [10]


def test_slice_data_overlap_info_reports_triple_overlap_pairs(
    run_in_tmpdir, triple_overlapping_slice_model
):
    """A triple-overlap pixel should report both overlapping cell pairs."""
    import openmc.lib

    with openmc.lib.TemporarySession(
        triple_overlapping_slice_model, output=False, args=["-c"]
    ):
        geom_data = _run_overlap_slice(triple_overlapping_slice_model)
        x, y = _first_overlap_pixel(geom_data)

        assert openmc.lib.slice_data_overlap_count(x, y) == 2

        records = _overlap_records(x, y)

        assert records == [
            (1, 2, 10),
            (1, 3, 10),
        ]


def test_slice_plot_overlap_data_is_empty_for_non_overlap_pixel(
    run_in_tmpdir, overlapping_slice_model
):
    """A normal plotted pixel should not produce GUI overlap details."""
    import openmc.lib

    with openmc.lib.TemporarySession(
        overlapping_slice_model, output=False, args=["-c"]
    ):
        geom_data = _run_overlap_slice(overlapping_slice_model)
        x, y = _first_pixel_with_cell_id(geom_data, 3)

        assert openmc.lib.slice_data_overlap_count(x, y) == 0

        cell1, cell2, universe = openmc.lib.slice_data_overlap_info(x, y)

        assert cell1.size == 0
        assert cell2.size == 0
        assert universe.size == 0