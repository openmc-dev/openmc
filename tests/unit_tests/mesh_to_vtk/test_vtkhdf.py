"""Tests for writing structured mesh data to VTKHDF files"""

from math import pi
from pathlib import Path

import h5py
import numpy as np
import pytest

import openmc


def _read_vtk(filename):
    """Read a legacy VTK or VTKHDF file, skipping the test without VTK"""
    vtk = pytest.importorskip("vtk")
    if Path(filename).suffix == ".vtkhdf":
        reader = vtk.vtkHDFReader()
    else:
        reader = vtk.vtkGenericDataObjectReader()
    reader.SetFileName(str(filename))
    reader.Update()
    return reader.GetOutput()


def _cell_data(grid, name):
    """Return a cell data array of a VTK grid as a NumPy array"""
    array = grid.GetCellData().GetArray(name)
    assert array is not None
    return np.array([array.GetTuple1(i)
                     for i in range(array.GetNumberOfTuples())])


def _rectilinear_mesh():
    mesh = openmc.RectilinearMesh()
    mesh.x_grid = [0., 1., 3.]
    mesh.y_grid = [0., 2., 3., 7.]
    mesh.z_grid = [-1., 0., 4., 5., 9.]
    return mesh


# Each mesh has a different number of elements along each axis so that a
# transposed write cannot match by accident. Except for the regular mesh, the
# element volumes also differ.
MESHES = {
    "regular": lambda: openmc.RegularMesh.from_domain(
        openmc.BoundingBox([0.0, 0.0, 0.0], [2.0, 3.0, 4.0]),
        dimension=[2, 3, 4],
    ),
    "rectilinear": _rectilinear_mesh,
    "cylindrical": lambda: openmc.CylindricalMesh(
        r_grid=[0.0, 1.0, 2.0],
        phi_grid=np.linspace(0.0, 2 * pi, 4),
        z_grid=[0.0, 1.0, 2.0, 3.0, 4.0],
    ),
    "spherical": lambda: openmc.SphericalMesh(
        r_grid=[0.0, 1.0, 2.0],
        theta_grid=np.linspace(0.0, pi, 4),
        phi_grid=np.linspace(0.0, 2 * pi, 6),
    ),
}
NON_REGULAR = ["rectilinear", "cylindrical", "spherical"]


@pytest.mark.parametrize("suffix", [".vtkhdf", ".vtk"])
@pytest.mark.parametrize("kind", list(MESHES))
def test_write_vtk_round_trip(run_in_tmpdir, kind, suffix):
    """VTK reads back every element value in mesh index order.

    VTK rejects VTKHDF files that declare an unsupported dataset type, so the
    files are read with VTK rather than only inspected with h5py. The legacy
    format is included so that both writers are held to the same ordering.
    """
    pytest.importorskip("vtk")
    mesh = MESHES[kind]()
    # distinct values per element so any reordering is detected
    data = np.arange(float(mesh.n_elements)).reshape(mesh.dimension)
    filename = f"round_trip_{kind}{suffix}"
    mesh.write_data_to_vtk(filename, {"values": data},
                           volume_normalization=False)

    grid = _read_vtk(filename)
    assert grid.GetNumberOfCells() == mesh.n_elements
    assert grid.GetNumberOfPoints() == np.prod([d + 1 for d in mesh.dimension])
    np.testing.assert_allclose(_cell_data(grid, "values"), data.T.ravel())


@pytest.mark.parametrize("kind", list(MESHES))
def test_write_vtkhdf_element_geometry(run_in_tmpdir, kind):
    """Each element read back by VTK must span the expected mesh vertices"""
    mesh = MESHES[kind]()
    mesh.write_data_to_vtk(datasets=None, filename=f"test_geom_{kind}.vtkhdf")

    grid = _read_vtk(f"test_geom_{kind}.vtkhdf")
    n_i, n_j, _ = mesh.dimension
    vertices = np.asarray(mesh.vertices)

    for cell_id in range(grid.GetNumberOfCells()):
        i = cell_id % n_i
        j = (cell_id // n_i) % n_j
        k = cell_id // (n_i * n_j)
        cell = grid.GetCell(cell_id)
        found = np.array([cell.GetPoints().GetPoint(n) for n in range(8)])
        expected = vertices[i:i+2, j:j+2, k:k+2].reshape(-1, 3)
        np.testing.assert_allclose(
            np.array(sorted(map(tuple, np.round(found, 8)))),
            np.array(sorted(map(tuple, np.round(expected, 8)))),
            atol=1e-8,
        )


def test_write_vtkhdf_regular_mesh(run_in_tmpdir):
    """A regular mesh is uniform, so it is written as VTKHDF ImageData."""
    mesh = openmc.RegularMesh()
    mesh.lower_left = [-5., -5., -5.]
    mesh.upper_right = [5., 5., 5.]
    mesh.dimension = [2, 3, 4]

    rng = np.random.default_rng(42)
    ref_data = rng.random(mesh.dimension)
    filename = "test_regular_mesh.vtkhdf"
    mesh.write_data_to_vtk(datasets={"data": ref_data}, filename=filename)

    assert Path(filename).exists()

    with h5py.File(filename, "r") as f:
        assert "VTKHDF" in f
        root = f["VTKHDF"]
        assert root.attrs["Type"] == b"ImageData"
        assert tuple(root.attrs["Version"]) == (2, 1)
        assert "CellData" in root
        assert "data" in root["CellData"]

        # ImageData is described by an extent, origin and spacing rather than
        # explicit point coordinates
        assert "Points" not in root
        np.testing.assert_array_equal(root.attrs["WholeExtent"],
                                      [0, 2, 0, 3, 0, 4])
        np.testing.assert_allclose(root.attrs["Origin"], [-5., -5., -5.])
        np.testing.assert_allclose(root.attrs["Spacing"], [5., 10/3, 2.5])

        # element data is stored with the last index varying fastest
        assert root["CellData"]["data"].shape == (4, 3, 2)


@pytest.mark.parametrize("kind", NON_REGULAR)
def test_write_vtkhdf_unstructured_grid(run_in_tmpdir, kind):
    """Other structured meshes are written as hexahedra on the mesh vertices."""
    mesh = MESHES[kind]()
    filename = f"unstructured_{kind}.vtkhdf"
    mesh.write_data_to_vtk(filename, {"values": np.ones(mesh.dimension)})

    n_cells = mesh.n_elements
    with h5py.File(filename, "r") as f:
        root = f["VTKHDF"]
        assert root.attrs["Type"] == b"UnstructuredGrid"
        assert tuple(root.attrs["Version"]) == (2, 1)
        np.testing.assert_allclose(
            root["Points"][()], np.asarray(mesh.vertices).reshape(-1, 3)
        )
        assert root["NumberOfPoints"][0] == np.prod(
            [d + 1 for d in mesh.dimension])
        assert root["NumberOfCells"][0] == n_cells
        assert root["NumberOfConnectivityIds"][0] == 8 * n_cells
        np.testing.assert_array_equal(root["Offsets"][()],
                                      np.arange(0, 8 * n_cells + 1, 8))
        # all elements are linear hexahedra
        np.testing.assert_array_equal(root["Types"][()], 12)
        assert root["CellData"]["values"].shape == (n_cells,)


def test_write_vtkhdf_volume_normalization(run_in_tmpdir):
    """Data are divided by element volumes by default and kept otherwise."""
    mesh = openmc.RegularMesh()
    mesh.lower_left = [0., 0., 0.]
    mesh.upper_right = [10., 10., 10.]
    mesh.dimension = [2, 2, 2]

    # non-uniform data so that ordering mistakes change the result
    data = np.arange(8, dtype=float).reshape(mesh.dimension)
    mesh.write_data_to_vtk("default.vtkhdf", {"flux": data})
    mesh.write_data_to_vtk("unnormalized.vtkhdf", {"flux": data},
                           volume_normalization=False)

    with h5py.File("default.vtkhdf", "r") as f:
        normalized = f["VTKHDF"]["CellData"]["flux"][()]
    with h5py.File("unnormalized.vtkhdf", "r") as f:
        unnormalized = f["VTKHDF"]["CellData"]["flux"][()]

    expected = data.T.ravel().reshape(2, 2, 2)
    # every element has a volume of 5 x 5 x 5 = 125
    np.testing.assert_allclose(normalized, expected / 125.0)
    np.testing.assert_allclose(unnormalized, expected)


@pytest.mark.parametrize("flat", [False, True], ids=["3d", "flat"])
@pytest.mark.parametrize("kind", NON_REGULAR)
def test_write_vtkhdf_nonuniform_volume_normalization(run_in_tmpdir, kind,
                                                      flat):
    """Each value must be divided by the volume of its own element.

    The element volumes differ, so pairing values and volumes in the wrong
    order changes the written data.
    """
    mesh = MESHES[kind]()
    assert np.unique(np.round(mesh.volumes, 12)).size > 1

    rng = np.random.default_rng(7)
    data = rng.random(mesh.dimension)
    # flat tally data is ordered with the first mesh index varying fastest
    dataset = data.T.ravel() if flat else data
    filename = f"nonuniform_{kind}.vtkhdf"
    mesh.write_data_to_vtk(filename, {"values": dataset})

    with h5py.File(filename, "r") as f:
        written = f["VTKHDF"]["CellData"]["values"][()]
    np.testing.assert_allclose(written, (data / mesh.volumes).T.ravel())


def _omitted_axis_mesh(kind):
    """Return a mesh whose last axis has a single element"""
    if kind == "regular":
        mesh = openmc.RegularMesh()
        mesh.lower_left = [0., 0., 0.]
        mesh.upper_right = [4., 3., 2.]
        mesh.dimension = [4, 3, 1]
    else:
        mesh = openmc.RectilinearMesh()
        mesh.x_grid = [0., 1., 3., 4., 7.]
        mesh.y_grid = [0., 2., 3., 7.]
        mesh.z_grid = [0., 2.]
    return mesh


@pytest.mark.parametrize("volume_normalization", [False, True])
@pytest.mark.parametrize("suffix", [".vtk", ".vtkhdf"])
@pytest.mark.parametrize("kind", ["regular", "rectilinear"])
def test_write_vtk_omitted_trailing_axis(run_in_tmpdir, kind, suffix,
                                         volume_normalization):
    """A dataset may omit trailing mesh axes of size one in either format.

    The values must be written in the same element order as a dataset with
    the full mesh shape.
    """
    pytest.importorskip("vtk")
    mesh = _omitted_axis_mesh(kind)
    data = np.arange(12, dtype=float).reshape(4, 3)
    filename = f"omitted_axis_{kind}{suffix}"
    mesh.write_data_to_vtk(filename, {"data": data},
                           volume_normalization=volume_normalization)

    expected = data.T.ravel()
    if volume_normalization:
        expected = expected / mesh.volumes.T.ravel()
    np.testing.assert_allclose(_cell_data(_read_vtk(filename), "data"),
                               expected)


@pytest.mark.parametrize(
    "upper_right,dimension,extent", [
        ([10.], [5], [0, 5, 0, 1, 0, 1]),
        ([10., 10.], [5, 3], [0, 5, 0, 3, 0, 1]),
    ],
    ids=["1d", "2d"],
)
def test_write_vtkhdf_low_dimensional_regular_mesh(run_in_tmpdir, upper_right,
                                                   dimension, extent):
    """1D and 2D regular meshes are padded to three dimensions."""
    mesh = openmc.RegularMesh()
    mesh.lower_left = [0.] * len(dimension)
    mesh.upper_right = upper_right
    mesh.dimension = dimension

    data = np.arange(float(mesh.n_elements)).reshape(dimension)
    filename = "low_dimensional.vtkhdf"
    mesh.write_data_to_vtk(filename, {"data": data})

    with h5py.File(filename, "r") as f:
        root = f["VTKHDF"]
        # a VTKHDF extent always has six entries; flat axes have one element
        np.testing.assert_array_equal(root.attrs["WholeExtent"], extent)
        saved = root["CellData"]["data"][()]

    # default volume normalization divides by the element volume
    expected = (data / np.prod(mesh.width)).T.ravel()
    np.testing.assert_allclose(
        saved, expected.reshape(extent[5], extent[3], extent[1]))

    assert _read_vtk(filename).GetNumberOfCells() == mesh.n_elements


def test_write_vtkhdf_multiple_datasets(run_in_tmpdir):
    """Test writing multiple datasets to VTKHDF format."""
    mesh = openmc.RegularMesh()
    mesh.lower_left = [0., 0., 0.]
    mesh.upper_right = [1., 1., 1.]
    mesh.dimension = [2, 2, 2]

    rng = np.random.default_rng(42)
    data1 = rng.random(mesh.dimension)
    data2 = rng.random(mesh.dimension)
    data3 = rng.random(mesh.dimension)

    filename = "test_multiple_datasets.vtkhdf"
    mesh.write_data_to_vtk(
        datasets={"flux": data1, "power": data2, "heating": data3},
        filename=filename,
        volume_normalization=False
    )

    assert Path(filename).exists()

    shape = (2, 2, 2)
    with h5py.File(filename, "r") as f:
        root = f["VTKHDF"]
        assert "flux" in root["CellData"]
        assert "power" in root["CellData"]
        assert "heating" in root["CellData"]

        np.testing.assert_allclose(
            root["CellData"]["flux"][()], data1.T.ravel().reshape(shape))
        np.testing.assert_allclose(
            root["CellData"]["power"][()], data2.T.ravel().reshape(shape))
        np.testing.assert_allclose(
            root["CellData"]["heating"][()], data3.T.ravel().reshape(shape))


def test_write_vtkhdf_invalid_data_shape(run_in_tmpdir):
    """Test that VTKHDF raises error for mismatched data shape."""
    mesh = openmc.RegularMesh()
    mesh.lower_left = [0., 0., 0.]
    mesh.upper_right = [1., 1., 1.]
    mesh.dimension = [2, 2, 2]

    # Create data with wrong shape but correct total size
    wrong_data = np.ones((4, 2, 1))
    filename = "test_invalid_shape.vtkhdf"

    with pytest.raises(ValueError, match="Cannot apply multidimensional dataset"):
        mesh.write_data_to_vtk(
            datasets={"bad_data": wrong_data},
            filename=filename
        )


@pytest.mark.parametrize("kind", list(MESHES))
def test_write_vtkhdf_curvilinear_warning(run_in_tmpdir, recwarn, kind):
    """VTKHDF output warns that curvilinear elements are not written.

    The option only changes the legacy output of cylindrical and spherical
    meshes, so only those meshes warn.
    """
    mesh = MESHES[kind]()
    filename = f"curvilinear_{kind}.vtkhdf"
    mesh.write_data_to_vtk(filename, curvilinear=True)

    curvilinear_warnings = [w for w in recwarn
                            if "Curvilinear elements" in str(w.message)]
    expect_warning = kind in ("cylindrical", "spherical")
    assert len(curvilinear_warnings) == int(expect_warning)
    assert Path(filename).exists()
