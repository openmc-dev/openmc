import os
import numpy as np
import pytest
import h5py
from openmc.plots import lost_particle_tracks_to_vtk

def make_tracks_h5(path, tracks):
    """Write a synthetic tracks.h5 with known coordinates.
    
    tracks : dict mapping dataset name to list of (x, y, z) tuples
    """
    dtype = np.dtype([
        ('r', [('x', '<f8'), ('y', '<f8'), ('z', '<f8')]),
        ('u', [('x', '<f8'), ('y', '<f8'), ('z', '<f8')]),
        ('E', '<f8'), ('wgt', '<f8'), ('time', '<f8'),
    ])
    with h5py.File(path, 'w') as f:
        for name, points in tracks.items():
            arr = np.zeros(len(points), dtype=dtype)
            for i, (x, y, z) in enumerate(points):
                arr[i]['r']['x'] = x
                arr[i]['r']['y'] = y
                arr[i]['r']['z'] = z
            f.create_dataset(name, data=arr)


def read_vtp_points(path):
    """Read point coordinates back from a .vtp file using VTK."""
    import vtk
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(path)
    reader.Update()
    polydata = reader.GetOutput()
    points = polydata.GetPoints()
    return [points.GetPoint(i) for i in range(points.GetNumberOfPoints())]


def test_point_coordinates_match(tmp_path):
    # Known coordinates to write and then verify round-trip
    expected_points = [(0.0, 0.0, 0.0), (1.0, 2.0, 3.0), (4.0, 5.0, 6.0)]
    tracks_file = str(tmp_path / "tracks.h5")
    output_dir  = str(tmp_path / "vtk_files")

    make_tracks_h5(tracks_file, {"track_1_1_1": expected_points})
    lost_particle_tracks_to_vtk(tracks_file=tracks_file, output_dir=output_dir)

    actual_points = read_vtp_points(os.path.join(output_dir, "track_1_1_1.vtp"))

    assert len(actual_points) == len(expected_points)
    for actual, expected in zip(actual_points, expected_points):
        assert actual == pytest.approx(expected, abs=1e-10)


def test_point_order_preserved(tmp_path):
    # Order of points along the track must be preserved exactly
    points = [(float(i), float(i*2), float(i*3)) for i in range(5)]
    tracks_file = str(tmp_path / "tracks.h5")
    output_dir  = str(tmp_path / "vtk_files")

    make_tracks_h5(tracks_file, {"track_1_1_1": points})
    lost_particle_tracks_to_vtk(tracks_file=tracks_file, output_dir=output_dir)

    actual = read_vtp_points(os.path.join(output_dir, "track_1_1_1.vtp"))
    for actual_pt, expected_pt in zip(actual, points):
        assert actual_pt == pytest.approx(expected_pt, abs=1e-10)


def test_multiple_tracks_coordinates(tmp_path):
    # Each track's points must appear in its own file without cross-contamination
    track_data = {
        "track_1_1_1": [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)],
        "track_1_1_2": [(9.0, 9.0, 9.0), (8.0, 8.0, 8.0)],
    }
    tracks_file = str(tmp_path / "tracks.h5")
    output_dir  = str(tmp_path / "vtk_files")

    make_tracks_h5(tracks_file, track_data)
    lost_particle_tracks_to_vtk(tracks_file=tracks_file, output_dir=output_dir)

    for name, expected_points in track_data.items():
        actual = read_vtp_points(os.path.join(output_dir, f"{name}.vtp"))
        assert len(actual) == len(expected_points)
        for actual_pt, expected_pt in zip(actual, expected_points):
            assert actual_pt == pytest.approx(expected_pt, abs=1e-10)