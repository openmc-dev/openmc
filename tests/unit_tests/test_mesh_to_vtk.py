import numpy as np
from pathlib import Path
import pytest

vtk = pytest.importorskip("vtk")
from vtk.util import numpy_support as nps

import openmc


regular_mesh = openmc.RegularMesh()
regular_mesh.lower_left = (0, 0, 0)
regular_mesh.upper_right = (1, 1, 1)
regular_mesh.dimension = [30, 20, 10]

rectilinear_mesh = openmc.RectilinearMesh()
rectilinear_mesh.x_grid = np.linspace(1, 2, num=30)
rectilinear_mesh.y_grid = np.linspace(1, 2, num=30)
rectilinear_mesh.z_grid = np.linspace(1, 2, num=30)

cylinder_mesh = openmc.CylindricalMesh()
cylinder_mesh.r_grid = np.linspace(1, 2, num=30)
cylinder_mesh.phi_grid = np.linspace(0, np.pi, num=50)
cylinder_mesh.z_grid = np.linspace(0, 1, num=30)

spherical_mesh = openmc.SphericalMesh()
spherical_mesh.r_grid = np.linspace(1, 2, num=30)
spherical_mesh.phi_grid = np.linspace(0, np.pi, num=50)
spherical_mesh.theta_grid = np.linspace(0, np.pi / 2, num=30)


@pytest.mark.parametrize("mesh", [cylinder_mesh, regular_mesh, rectilinear_mesh, spherical_mesh])
def test_write_data_to_vtk(mesh, tmpdir):
    # BUILD
    filename = Path(tmpdir) / "out.vtk"

    # use mesh element volumes as data to check volume-normalization ordering
    # kji (i changing fastest) orering is expected for input data
    # by using the volumes transposed as the data here, we can ensure the
    # normalization is happening correctly
    data = mesh.volumes.T

    # RUN
    mesh.write_data_to_vtk(filename=filename, datasets={"label1": data, "label2": data})

    # TEST
    assert filename.is_file()

    # read file
    reader = vtk.vtkStructuredGridReader()
    reader.SetFileName(str(filename))
    reader.Update()

    # check name of datasets
    vtk_grid = reader.GetOutput()
    array1 = vtk_grid.GetCellData().GetArray(0)
    array2 = vtk_grid.GetCellData().GetArray(1)

    assert array1.GetName() == "label1"
    assert array2.GetName() == "label2"

    # check size of datasets
    data1 = nps.vtk_to_numpy(array1)
    data2 = nps.vtk_to_numpy(array2)
    assert data1.size == data.size
    assert data2.size == data.size

    assert all(data1 == data2)
    assert all(data1 == 1.0)


@pytest.mark.parametrize("mesh", [cylinder_mesh, regular_mesh, rectilinear_mesh, spherical_mesh])
def test_write_data_to_vtk_size_mismatch(mesh):
    """Checks that an error is raised when the size of the dataset
    doesn't match the mesh number of cells

    Parameters
    ----------
    mesh : openmc.StructuredMesh
        The mesh to test
    """
    right_size = mesh.num_mesh_cells
    data = np.random.random(right_size + 1)

    # Error message has \ in to escape characters that are otherwise recognized
    # by regex. These are needed to make the test string match the error message
    # string when using the match argument as that uses regular expression
    expected_error_msg = (
        fr"The size of the dataset 'label' \({len(data)}\) should be equal to "
        fr"the number of mesh cells \({mesh.num_mesh_cells}\)"
    )
    with pytest.raises(ValueError, match=expected_error_msg):
        mesh.write_data_to_vtk(filename="out.vtk", datasets={"label": data})

def test_write_data_to_vtk_round_trip(run_in_tmpdir):
    cmesh = openmc.CylindricalMesh()
    cmesh.r_grid = (0.0, 1.0, 2.0)
    cmesh.z_grid = (0.0, 2.0, 4.0, 5.0)
    cmesh.phi_grid = (0.0, 3.0, 6.0)

    smesh = openmc.SphericalMesh()
    smesh.r_grid = (0.0, 1.0, 2.0)
    smesh.theta_grid = (0.0, 2.0, 4.0, 5.0)
    smesh.phi_grid = (0.0, 3.0, 6.0)

    rmesh = openmc.RegularMesh()
    rmesh.lower_left = (0.0, 0.0, 0.0)
    rmesh.upper_right = (1.0, 3.0, 5.0)
    rmesh.dimension = (2, 1, 6)

    for mesh in [smesh, cmesh, rmesh]:

        filename = "mesh.vtk"
        data = np.array([1.0] * 12)  # there are 12 voxels in each mesh
        mesh.write_data_to_vtk(
            filename=filename,
            datasets={"normalized": data},
            volume_normalization=True
        )

        reader = vtk.vtkStructuredGridReader()
        reader.SetFileName(filename)
        reader.ReadAllFieldsOn()
        reader.Update()

        cell_data = reader.GetOutput().GetCellData()
        uniform_array = cell_data.GetArray("normalized")
        num_tuples = uniform_array.GetNumberOfTuples()
        vtk_values = [uniform_array.GetValue(i) for i in range(num_tuples)]

        # checks that the vtk cell values are equal to the data / mesh volumes
        assert np.allclose(vtk_values, data / mesh.volumes.T.flatten())

        mesh.write_data_to_vtk(
            filename=filename,
            datasets={"not_normalized": data},
            volume_normalization=False,
        )

        reader = vtk.vtkStructuredGridReader()
        reader.SetFileName(filename)
        reader.ReadAllFieldsOn()
        reader.Update()

        cell_data = reader.GetOutput().GetCellData()
        uniform_array = cell_data.GetArray("not_normalized")
        num_tuples = uniform_array.GetNumberOfTuples()
        vtk_values = [uniform_array.GetValue(i) for i in range(num_tuples)]

        # checks that the vtk cell values are equal to the data
        assert np.array_equal(vtk_values, data)
