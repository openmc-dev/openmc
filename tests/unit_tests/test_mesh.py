import openmc
import pytest
import numpy as np

@pytest.mark.parametrize("val_left,val_right", [(0, 0), (-1., -1.), (2.0, 2)])
def test_raises_error_when_flat(val_left, val_right):
    """Checks that an error is raised when a mesh is flat"""
    mesh = openmc.RegularMesh()

    # Same X
    with pytest.raises(ValueError):
        mesh.lower_left = [val_left, -25, -25]
        mesh.upper_right = [val_right, 25, 25]

    with pytest.raises(ValueError):
        mesh.upper_right = [val_right, 25, 25]
        mesh.lower_left = [val_left, -25, -25]

    # Same Y
    with pytest.raises(ValueError):
        mesh.lower_left = [-25, val_left, -25]
        mesh.upper_right = [25, val_right, 25]

    with pytest.raises(ValueError):
        mesh.upper_right = [25, val_right, 25]
        mesh.lower_left = [-25, val_left, -25]

    # Same Z
    with pytest.raises(ValueError):
        mesh.lower_left = [-25, -25, val_left]
        mesh.upper_right = [25, 25, val_right]

    with pytest.raises(ValueError):
        mesh.upper_right = [25, 25, val_right]
        mesh.lower_left = [-25, -25, val_left]


def test_get_index_where():

    mesh = openmc.RegularMesh()
    mesh.dimension = [20, 30, 50]
    mesh.lower_left = (-2, -2, -2)
    mesh.upper_right = (2, 2, 2)

    #  checks center of mesh
    assert mesh.get_index_where(value=0, basis='xy') == 24
    assert mesh.get_index_where(value=0, basis='xz') == 14
    assert mesh.get_index_where(value=0, basis='yz') == 10

    #  checks edge of mesh
    mesh.upper_right = (4, 6, 8)
    assert mesh.get_index_where(value=8, basis='xy') == 49
    assert mesh.get_index_where(value=6, basis='xz') == 29
    assert mesh.get_index_where(value=4, basis='yz') == 19

    # checks outside of the mesh
    with pytest.raises(ValueError):
        assert mesh.get_index_where(value=9, basis='xy')
    with pytest.raises(ValueError):
        assert mesh.get_index_where(value=7, basis='xz')
    with pytest.raises(ValueError):
        assert mesh.get_index_where(value=5, basis='yz')


def test_mesh_bounding_box():
    mesh = openmc.RegularMesh()
    mesh.lower_left = [-2, -3 ,-5]
    mesh.upper_right = [2, 3, 5]
    bb = mesh.bounding_box
    assert isinstance(bb, openmc.BoundingBox)
    np.testing.assert_array_equal(bb.lower_left, np.array([-2, -3 ,-5]))
    np.testing.assert_array_equal(bb.upper_right, np.array([2, 3, 5]))
