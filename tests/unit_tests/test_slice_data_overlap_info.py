import pytest
import numpy as np
import openmc
import openmc.lib

# Sentinel value matching _OVERLAP in plotmodel.py and OVERLAP in plot.cpp
_OVERLAP = -3


@pytest.fixture
def overlap_model(tmp_path):
    openmc.reset_auto_ids()
    # Three cylinders: cyl1 and cyl2 overlap near x=0, cyl2 and cyl3 overlap
    # near x=4. This gives us two spatially distinct overlap regions in one model.
    mat_water = openmc.Material()
    mat_water.add_nuclide('H1', 2.0)
    mat_water.add_nuclide('O16', 1.0)
    mat_water.set_density('g/cm3', 1.0)

    mat_iron = openmc.Material()
    mat_iron.add_nuclide('Fe56', 1.0)
    mat_iron.set_density('g/cm3', 7.87)

    mat_lead = openmc.Material()
    mat_lead.add_nuclide('Pb208', 1.0)
    mat_lead.set_density('g/cm3', 11.34)

    # cyl1 and cyl2 overlap on the left, cyl2 and cyl3 overlap on the right
    cyl1 = openmc.ZCylinder(x0=-2.0, y0=0.0, r=2.5)
    cyl2 = openmc.ZCylinder(x0=0.0, y0=0.0, r=2.5)
    cyl3 = openmc.ZCylinder(x0=2.0, y0=0.0, r=2.5)
    boundary = openmc.Sphere(r=20.0, boundary_type='vacuum')

    cell1 = openmc.Cell(region=-cyl1, fill=mat_water)
    cell2 = openmc.Cell(region=-cyl2, fill=mat_iron)
    cell3 = openmc.Cell(region=-cyl3, fill=mat_lead)
    cell_outside = openmc.Cell(region=+cyl1 & +cyl2 & +cyl3 & -boundary)

    geometry = openmc.Geometry([cell1, cell2, cell3, cell_outside])

    settings = openmc.Settings()
    settings.run_mode = 'fixed source'
    settings.particles = 100
    settings.batches = 1
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Point((0, 0, 0))
    )

    model = openmc.Model(geometry=geometry, settings=settings)
    model.export_to_xml(tmp_path)
    return tmp_path


@pytest.fixture
def lib_session(overlap_model, monkeypatch):
    # Set OMP_NUM_THREADS=1 BEFORE init so OpenMP reads it at thread pool
    # creation time.
    monkeypatch.setenv('OMP_NUM_THREADS', '1')
    # Change to the model directory so openmc.lib.init finds the XML files
    monkeypatch.chdir(overlap_model)
    openmc.lib.init()
    yield  # tests run here
    openmc.lib.finalize()  # always runs, even if the test fails


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


def test_overlap_pixels_present(lib_session):
    # mat_ids (geom_data[:,:,2]) should contain values more negative than -3
    # when overlaps are detected. RasterData::set_overlap in plot.cpp
    # encodes each unique overlap as OVERLAP - overlap_idx - 1 into slot 2.
    geom_data = run_slice()
    assert np.any(geom_data[:, :, 2] < _OVERLAP)


def test_no_overlap_pixels_when_disabled(lib_session):
    # With show_overlaps=False, set_overlap is never called so no encoded
    # indices should appear in mat_ids
    geom_data = run_slice(show_overlaps=False)
    assert not np.any(geom_data[:, :, 2] < _OVERLAP)


def test_overlap_map_populated(lib_session):
    # After a slice with overlaps, overlap_keys in plot.cpp should have
    # at least one entry. slice_data_overlap_info reads from that vector
    run_slice()
    _, n = openmc.lib.slice_data_overlap_info()
    assert n > 0


def test_two_distinct_overlap_regions(lib_session):
    # With three cylinders arranged so that cyl1/cyl2 overlap on the left and
    # cyl2/cyl3 overlap on the right, we expect at least 2 entries in the
    # overlap_map — one for each distinct overlapping pair of cells.
    run_slice()
    _, n = openmc.lib.slice_data_overlap_info()
    assert n >= 2, f"Expected at least 2 overlap entries, got {n}"


def test_overlap_map_empty_when_disabled(lib_session):
    # overlap_keys is cleared at the start of every openmc_slice_data call
    # and never written to when show_overlaps=False, so n should be 0
    run_slice(show_overlaps=False)
    _, n = openmc.lib.slice_data_overlap_info()
    assert n == 0


def test_overlap_info_correct(lib_session):
    # Each entry in overlap_info is a (universe_id, cell1_id, cell2_id) triple
    # packed as [i*3, i*3+1, i*3+2]. Verify values make sense for our model:
    # only universe 1 exists, and only cells 1, 2, 3 can overlap.
    run_slice()
    overlap_info, n = openmc.lib.slice_data_overlap_info()
    assert n > 0
    for i in range(n):
        universe_id = int(overlap_info[i * 3])
        cell1_id    = int(overlap_info[i * 3 + 1])
        cell2_id    = int(overlap_info[i * 3 + 2])
        assert universe_id == 1
        assert cell1_id in {1, 2, 3}
        assert cell2_id in {1, 2, 3}
        assert cell1_id != cell2_id


def test_overlap_pairs_are_correct(lib_session):
    # The two expected overlapping pairs are {1,2} (cyl1/cyl2) and {2,3}
    # (cyl2/cyl3). Verify both are present in the overlap_map.
    run_slice()
    overlap_info, n = openmc.lib.slice_data_overlap_info()
    assert n >= 2

    pairs = set()
    for i in range(n):
        c1 = int(overlap_info[i * 3 + 1])
        c2 = int(overlap_info[i * 3 + 2])
        pairs.add(frozenset({c1, c2}))

    assert frozenset({1, 2}) in pairs, "Expected overlap between cells 1 and 2"
    assert frozenset({2, 3}) in pairs, "Expected overlap between cells 2 and 3"


def test_overlap_index_decodes_correctly(lib_session):
    # Core end-to-end test mirroring what getIDinfo does in plotgui.py:
    # read raw mat_id from geom_data, decode the overlap index via
    # _OVERLAP - id - 1, and verify it maps to a valid overlap_map entry.
    geom_data = run_slice()
    overlap_info, n = openmc.lib.slice_data_overlap_info()
    mat_ids = geom_data[:, :, 2]

    # Check every unique encoded overlap value found in the image
    for raw in np.unique(mat_ids[mat_ids < _OVERLAP]):
        idx = int(_OVERLAP - int(raw) - 1)
        # Decoded index must be within the range returned by overlap_info
        assert 0 <= idx < n, f"Decoded index {idx} out of range [0, {n})"
        universe_id = int(overlap_info[idx * 3])
        cell1_id    = int(overlap_info[idx * 3 + 1])
        cell2_id    = int(overlap_info[idx * 3 + 2])
        assert universe_id == 1
        assert cell1_id != cell2_id
        assert cell1_id in {1, 2, 3}
        assert cell2_id in {1, 2, 3}
