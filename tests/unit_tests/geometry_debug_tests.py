import numpy as np
import pytest
import warnings
from collections import deque
from openmc.model import Model  # adjust import as needed

_NOT_FOUND = -2
_CELL_ID = 1


def make_grid(h, w, fill=_CELL_ID):
    """Helper: make a fully-defined grid."""
    return np.full((h, w), fill, dtype=np.int32)


# --- _classify_undefined_regions tests (no openmc.lib needed) ---

def test_no_undefined_pixels():
    # All defined: all three outputs should be empty
    cell_ids = make_grid(5, 5)
    undefined, outside, internal = Model._classify_undefined_regions(cell_ids)
    assert not undefined.any()
    assert not outside.any()
    assert not internal.any()


def test_none_input():
    undefined, outside, internal = Model._classify_undefined_regions(None)
    assert undefined is None
    assert outside is None
    assert internal is None


def test_border_undefined_is_outside():
    # Undefined pixel on the border should be classified as outside
    cell_ids = make_grid(5, 5)
    cell_ids[0, 0] = _NOT_FOUND
    undefined, outside, internal = Model._classify_undefined_regions(cell_ids)
    assert outside[0, 0]
    assert not internal[0, 0]


def test_interior_undefined_is_internal():
    # Undefined pixel fully surrounded by defined pixels is internal
    cell_ids = make_grid(5, 5)
    cell_ids[2, 2] = _NOT_FOUND
    undefined, outside, internal = Model._classify_undefined_regions(cell_ids)
    assert internal[2, 2]
    assert not outside[2, 2]


def test_connected_border_undefined_floods_inward():
    # A chain of undefined pixels from the border should all be outside
    cell_ids = make_grid(5, 5)
    cell_ids[0, 0] = _NOT_FOUND
    cell_ids[1, 0] = _NOT_FOUND
    cell_ids[2, 0] = _NOT_FOUND
    _, outside, internal = Model._classify_undefined_regions(cell_ids)
    assert outside[0, 0] and outside[1, 0] and outside[2, 0]
    assert not internal.any()


def test_internal_not_connected_to_border():
    # A ring of defined pixels isolates internal undefined pixels
    cell_ids = make_grid(7, 7)
    # Create a hollow interior undefined region
    cell_ids[2:5, 2:5] = _NOT_FOUND
    cell_ids[3, 3] = _CELL_ID  # plug the center back
    _, outside, internal = Model._classify_undefined_regions(cell_ids)
    # The ring pixels (2:5, 2:5) minus center are internal
    assert internal[2, 2]
    assert not outside[2, 2]


def test_warn_when_all_undefined_are_internal():
    # If undefined pixels exist but none touch the border, warn
    cell_ids = make_grid(5, 5)
    cell_ids[2, 2] = _NOT_FOUND
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        Model._classify_undefined_regions(cell_ids)
        assert len(w) == 1
        assert "internal" in str(w[0].message).lower()


def test_undefined_and_outside_are_disjoint():
    cell_ids = make_grid(5, 5)
    cell_ids[0, 2] = _NOT_FOUND  # border
    cell_ids[2, 2] = _NOT_FOUND  # interior
    _, outside, internal = Model._classify_undefined_regions(cell_ids)
    # outside and internal must not overlap
    assert not (outside & internal).any()


def test_all_undefined_grid():
    # Every pixel undefined: all should be outside (all border-connected)
    cell_ids = np.full((4, 4), _NOT_FOUND, dtype=np.int32)
    undefined, outside, internal = Model._classify_undefined_regions(cell_ids)
    assert undefined.all()
    assert outside.all()
    assert not internal.any()


# --- geometry_debug integration tests (require openmc.lib) ---

@pytest.fixture
def overlap_model(tmp_path):
    import openmc
    mat_water = openmc.Material(material_id=1, name='water')
    mat_water.add_nuclide('H1', 2.0)
    mat_water.add_nuclide('O16', 1.0)
    mat_water.set_density('g/cm3', 1.0)

    mat_iron = openmc.Material(material_id=2, name='iron')
    mat_iron.add_nuclide('Fe56', 1.0)
    mat_iron.set_density('g/cm3', 7.87)

    cyl1 = openmc.ZCylinder(x0=-1.0, y0=0.0, r=2.5)
    cyl2 = openmc.ZCylinder(x0= 1.0, y0=0.0, r=2.5)
    boundary = openmc.Sphere(r=15.0, boundary_type='vacuum')

    cell1 = openmc.Cell(cell_id=1, region=-cyl1, fill=mat_water)
    cell2 = openmc.Cell(cell_id=2, region=-cyl2, fill=mat_iron)
    cell_outside = openmc.Cell(cell_id=3, region=~(-cyl1 | -cyl2) & -boundary)

    universe = openmc.Universe(universe_id=1, cells=[cell1, cell2, cell_outside])
    geometry = openmc.Geometry(universe)

    settings = openmc.Settings()
    settings.run_mode = 'fixed source'
    settings.particles = 100
    settings.batches = 1
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Point((0, 0, 0))
    )

    model = openmc.Model(
        geometry=geometry,
        materials=openmc.Materials([mat_water, mat_iron]),
        settings=settings,
    )
    model.export_to_xml(tmp_path)
    return model, tmp_path


def test_geometry_debug_finds_overlaps(overlap_model, monkeypatch, tmp_path):
    model, path = overlap_model
    monkeypatch.setenv('OMP_NUM_THREADS', '1')
    monkeypatch.chdir(path)

    result = model.geometry_debug(
        lower_left=(-5.0, -5.0, -1.0),
        upper_right=(5.0, 5.0, 1.0),
        n_samples=(20, 20, 1),
    )

    assert result['n_overlap_samples'] > 0
    assert len(result['overlap_points']) > 0
    # Each point must be within the bounding box
    for pt in result['overlap_points']:
        assert -5.0 <= pt[0] <= 5.0
        assert -5.0 <= pt[1] <= 5.0
        assert -1.0 <= pt[2] <= 1.0


def test_geometry_debug_result_keys(overlap_model, monkeypatch, tmp_path):
    model, path = overlap_model
    monkeypatch.setenv('OMP_NUM_THREADS', '1')
    monkeypatch.chdir(path)

    result = model.geometry_debug(
        lower_left=(-5.0, -5.0, -1.0),
        upper_right=(5.0, 5.0, 1.0),
        n_samples=(10, 10, 1),
    )

    expected_keys = {
        'n_overlap_samples', 'n_undefined_samples',
        'overlap_points', 'undefined_points',
        'n_more_overlap_points', 'n_more_undefined_points',
    }
    assert expected_keys == set(result.keys())


def test_geometry_debug_n_samples_scalar_vs_tuple(overlap_model, monkeypatch, tmp_path):
    # Scalar and equivalent tuple should give the same result
    model, path = overlap_model
    monkeypatch.setenv('OMP_NUM_THREADS', '1')
    monkeypatch.chdir(path)

    r1 = model.geometry_debug((-5, -5, -1), (5, 5, 1), n_samples=5)
    r2 = model.geometry_debug((-5, -5, -1), (5, 5, 1), n_samples=(5, 5, 5))
    assert r1['n_overlap_samples'] == r2['n_overlap_samples']


def test_geometry_debug_invalid_n_samples(overlap_model, tmp_path):
    model, _ = overlap_model
    with pytest.raises(ValueError, match="length-3"):
        model.geometry_debug((-5, -5, -1), (5, 5, 1), n_samples=(5, 5))


def test_geometry_debug_max_examples_capped(overlap_model, monkeypatch, tmp_path):
    # With high resolution, more than 10 overlap samples but only 10 points returned
    model, path = overlap_model
    monkeypatch.setenv('OMP_NUM_THREADS', '1')
    monkeypatch.chdir(path)

    result = model.geometry_debug(
        lower_left=(-5.0, -5.0, -1.0),
        upper_right=(5.0, 5.0, 1.0),
        n_samples=(50, 50, 1),
    )

    assert len(result['overlap_points']) <= 10
    assert len(result['undefined_points']) <= 10
    assert result['n_more_overlap_points'] == max(
        0, result['n_overlap_samples'] - len(result['overlap_points'])
    )