import numpy as np
from os.path import exists
import pytest

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
    filename = tmpdir / "out.vtk"

    data = np.random.random(mesh.dimension[0]*mesh.dimension[1]*mesh.dimension[2])

    mesh.write_data_to_vtk(filename=filename, datasets={"label": data})
    assert exists(filename)