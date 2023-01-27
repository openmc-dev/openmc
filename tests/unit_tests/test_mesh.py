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


def test_corner_none_returns_false():
    """Checks mesh is not considered flat when one
    corner is None
    """
    mesh = openmc.RegularMesh()
    mesh.lower_left = [-25, -25, -25]
    assert not mesh.is_flat()

    mesh = openmc.RegularMesh()
    mesh.upper_right = [-25, -25, -25]
    assert not mesh.is_flat()

    # test with np array
    mesh = openmc.RegularMesh()
    mesh.lower_left = np.array([-25, -25, -25])
    assert not mesh.is_flat()

    mesh = openmc.RegularMesh()
    mesh.upper_right = np.array([-25, -25, -25])
    assert not mesh.is_flat()