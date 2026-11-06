import numpy as np
import openmc
from openmc.examples import pwr_pin_cell


def test_slice_data_basic(run_in_tmpdir):
    """Test basic slice_data functionality."""
    model = pwr_pin_cell()
    geom_data, prop_data = model.slice_data(
        origin=(0, 0, 0),
        width=(1.0, 1.0),
        pixels=(100, 100),
        basis='xy'
    )

    # Without filter, should have 3 fields
    assert geom_data.shape == (100, 100, 3)
    assert geom_data.dtype == np.int32
    assert prop_data.shape == (100, 100, 2)
    assert prop_data.dtype == np.float64

    # Check we have valid geometry
    assert np.any(geom_data[:, :, 0] >= 0)  # Valid cell IDs
    assert np.any(prop_data[:, :, 0] > 0)   # Valid temperatures


def test_slice_data_no_properties(run_in_tmpdir):
    """Test slice_data without property data."""
    model = pwr_pin_cell()
    geom_data, prop_data = model.slice_data(
        origin=(0, 0, 0),
        width=(1.0, 1.0),
        pixels=(50, 50),
        include_properties=False
    )

    # Without filter, should have 3 fields
    assert geom_data.shape == (50, 50, 3)
    assert prop_data is None


def test_slice_data_with_filter(run_in_tmpdir):
    """Test slice_data with a cell filter."""
    model = pwr_pin_cell()
    cell_ids = [c.id for c in model.geometry.get_all_cells().values()]
    cell_filter = openmc.CellFilter(cell_ids)

    geom_data, _ = model.slice_data(
        origin=(0, 0, 0),
        width=(1.0, 1.0),
        pixels=(50, 50),
        filter=cell_filter,
        include_properties=False
    )

    # With filter, should have 4 fields
    assert geom_data.shape == (50, 50, 4)

    # Filter bin index should be populated where cells exist
    filter_bins = geom_data[:, :, 3]
    valid_cells = geom_data[:, :, 0] >= 0
    assert np.any(filter_bins[valid_cells] >= 0)


def test_slice_data_overlaps(run_in_tmpdir):
    """Test slice_data with overlap detection."""
    model = pwr_pin_cell()
    geom_data, _ = model.slice_data(
        origin=(0, 0, 0),
        width=(1.0, 1.0),
        pixels=(50, 50),
        show_overlaps=True,
        include_properties=False
    )

    # Without filter, should have 3 fields
    assert geom_data.shape == (50, 50, 3)
    # Check for overlap markers (-3) if any exist
    # Note: This test may pass without finding overlaps if geometry is correct


def test_slice_data_overlaps_with_filter(run_in_tmpdir):
    """Test that overlaps don't overwrite filter bin data."""
    model = pwr_pin_cell()
    cell_ids = [c.id for c in model.geometry.get_all_cells().values()]
    cell_filter = openmc.CellFilter(cell_ids)

    geom_data, _ = model.slice_data(
        origin=(0, 0, 0),
        width=(1.0, 1.0),
        pixels=(50, 50),
        filter=cell_filter,
        show_overlaps=True,
        include_properties=False
    )

    assert geom_data.shape == (50, 50, 4)

    # If any overlaps exist, verify filter bin is still valid (not -3)
    overlap_pixels = geom_data[:, :, 0] == -3
    if np.any(overlap_pixels):
        # Filter bins at overlap locations should NOT be -3
        filter_bins_at_overlaps = geom_data[overlap_pixels, 3]
        assert not np.all(filter_bins_at_overlaps == -3), \
            "Filter bins should be preserved even where overlaps are detected"


def test_slice_data_different_bases(run_in_tmpdir):
    """Test slice_data with different basis planes."""
    model = pwr_pin_cell()

    for basis in ['xy', 'xz', 'yz']:
        geom_data, prop_data = model.slice_data(
            origin=(0, 0, 0),
            width=(1.0, 1.0),
            pixels=(25, 25),
            basis=basis
        )

        assert geom_data.shape == (25, 25, 3)
        assert prop_data.shape == (25, 25, 2)


def test_slice_data_oriented_spans(run_in_tmpdir):
    """Test slice_data with oriented span vectors."""
    model = pwr_pin_cell()

    geom_data, prop_data = model.slice_data(
        origin=(0, 0, 0),
        u_span=(1.0, 0.0, 0.0),
        v_span=(0.0, 0.0, 1.0),
        pixels=(25, 25)
    )

    assert geom_data.shape == (25, 25, 3)
    assert prop_data.shape == (25, 25, 2)


def test_slice_data_level(run_in_tmpdir):
    """Test slice_data with specific universe level."""
    model = pwr_pin_cell()
    geom_data, _ = model.slice_data(
        origin=(0, 0, 0),
        width=(1.0, 1.0),
        pixels=(50, 50),
        level=0,  # Root universe only
        include_properties=False
    )

    assert geom_data.shape == (50, 50, 3)


def test_id_map_reverted(run_in_tmpdir):
    """Test that id_map returns 3D array without filter support."""
    model = pwr_pin_cell()
    id_data = model.id_map(
        origin=(0, 0, 0),
        width=(1.0, 1.0),
        pixels=(50, 50),
        basis='xy'
    )

    # Should have 3 fields (cell_id, cell_instance, material_id)
    assert id_data.shape == (50, 50, 3)
    assert id_data.dtype == np.int32

    # Check valid data
    assert np.any(id_data[:, :, 0] >= 0)  # Valid cell IDs
