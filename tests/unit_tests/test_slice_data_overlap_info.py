import pytest
import numpy as np
import openmc
import openmc.lib

# Sentinel value matching _OVERLAP in plotmodel.py and OVERLAP in plot.cpp
_OVERLAP = -3


@pytest.fixture(scope='module')
def overlap_model():
    openmc.reset_auto_ids()

    # Three cylinders: cyl1 and cyl2 overlap near x=0, cyl2 and cyl3 overlap
    # near x=4. This gives us two spatially distinct overlap regions in one model.
    mat1 = openmc.Material(components={'H1': 1.0})
    mat2 = openmc.Material(components={'H1': 1.0})
    mat3 = openmc.Material(components={'H1': 1.0})

    # cyl1 and cyl2 overlap on the left, cyl2 and cyl3 overlap on the right
    cyl1 = openmc.ZCylinder(x0=-2.0, r=2.5)
    cyl2 = openmc.ZCylinder(x0=0.0, r=2.5)
    cyl3 = openmc.ZCylinder(x0=2.0, r=2.5)
    boundary = openmc.Sphere(r=20.0, boundary_type='vacuum')
    cell1 = openmc.Cell(region=-cyl1, fill=mat1)
    cell2 = openmc.Cell(region=-cyl2, fill=mat2)
    cell3 = openmc.Cell(region=-cyl3, fill=mat3)
    cell_outside = openmc.Cell(region=+cyl1 & +cyl2 & +cyl3 & -boundary)
    geometry = openmc.Geometry([cell1, cell2, cell3, cell_outside])

    settings = openmc.Settings()
    settings.run_mode = 'fixed source'
    settings.particles = 100
    settings.batches = 1
    model = openmc.Model(geometry=geometry, settings=settings)

    with openmc.lib.TemporarySession(model, args=['-s', '1']):
        yield


def run_slice(origin=(0.0, 0.0, 0.0), width=(10.0, 6.0), show_overlaps=True):
    # Helper that runs a slice over a region covering both overlap zones
    geom_data, _ = openmc.lib.slice_data(
        origin=origin,
        width=width,
        basis='xy',
        pixels=(100, 60),
        show_overlaps=show_overlaps,
        include_properties=False,
    )
    return geom_data


def test_overlaps_enabled(overlap_model):
    # Run a single slice with overlap detection enabled and check all
    # expected properties in one pass.
    geom_data = run_slice()
    overlap_info = openmc.lib.slice_data_overlap_info()
    n = overlap_info.shape[0]
    cell_ids = geom_data[:, :, 0]

    # cell_ids should contain values more negative than _OVERLAP; RasterData
    # encodes each unique overlap as OVERLAP - overlap_idx - 1 into slot 2.
    assert np.any(cell_ids < _OVERLAP)

    # overlap_keys should have 2 entries for the two distinct overlapping
    # cylinder pairs in this model.
    assert n == 2, f"Expected exactly 2 overlap entries, got {n}"

    # Each entry is a (universe_id, cell1_id, cell2_id) triple; verify values.
    for i in range(n):
        universe_id = int(overlap_info[i, 0])
        cell1_id = int(overlap_info[i, 1])
        cell2_id = int(overlap_info[i, 2])
        assert universe_id == 1
        assert cell1_id in {1, 2, 3}
        assert cell2_id in {1, 2, 3}
        assert cell1_id != cell2_id


def test_overlaps_disabled(overlap_model):
    # With show_overlaps=False, set_overlap is never called and overlap_keys
    # is never written to, so the image and map should both be clean.
    geom_data = run_slice(show_overlaps=False)
    overlap_info = openmc.lib.slice_data_overlap_info()
    n = overlap_info.shape[0]

    assert not np.any(geom_data[:, :, 2] < _OVERLAP)
    assert n == 0
