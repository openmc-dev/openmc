import numpy as np
import openmc.lib
import pytest


@pytest.mark.parametrize('rotation', [
    np.eye(3),
    np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]]),
])
def test_mesh_filter_rotation_matrix(rotation):
    # Allocating a mesh/filter through the C API does not require transport.
    mesh = openmc.lib.RegularMesh()
    mesh.dimension = (2, 3, 4)
    mesh.set_parameters(lower_left=(0., 0., 0.), upper_right=(2., 3., 4.))
    mesh_filter = openmc.lib.MeshFilter(mesh)
    try:
        np.testing.assert_array_equal(mesh_filter.rotation, np.zeros(3))
        mesh_filter.rotation = (0., 0., 90.)
        np.testing.assert_array_equal(mesh_filter.rotation, (0., 0., 90.))

        mesh_filter.rotation = rotation
        np.testing.assert_array_equal(mesh_filter.rotation, rotation)
    finally:
        openmc.lib.finalize()
