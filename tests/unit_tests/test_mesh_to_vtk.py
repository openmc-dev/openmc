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

    data = np.random.random(mesh.num_mesh_cells)

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
    assert nps.vtk_to_numpy(array1).size == data.size
    assert nps.vtk_to_numpy(array2).size == data.size


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
        f"The size of the dataset 'label' \({len(data)}\) should be equal to "
        f"the number of mesh cells \({mesh.num_mesh_cells}\)"
    )
    with pytest.raises(ValueError, match=expected_error_msg):
        mesh.write_data_to_vtk(filename="out.vtk", datasets={"label": data})
