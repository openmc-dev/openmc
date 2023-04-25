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


def test_get_data_slice():

    mesh = openmc.RegularMesh()
    mesh.dimension = [2, 3, 5]

    data = np.linspace(1, 30, 30)

    xy_slice = mesh.get_data_slice(dataset=data, basis='xy', slice_index=0)
    assert xy_slice.shape == (3, 2)
    np.testing.assert_array_equal(
        xy_slice,
        np.array(
            [
                [5., 6.],
                [3., 4.],
                [1., 2.]
            ]
        )
    )

    yz_slice = mesh.get_data_slice(dataset=data, basis='yz', slice_index=0)
    assert yz_slice.shape == (5, 3)
    np.testing.assert_array_equal(
        yz_slice,
        np.array(
            [
                [25., 27., 29.],
                [19., 21., 23.],
                [13., 15., 17.],
                [7.,  9., 11.],
                [1.,  3.,  5.]
            ]
        )
    )

    xz_slice = mesh.get_data_slice(dataset=data, basis='xz', slice_index=0)
    assert xz_slice.shape == (5, 2)
    np.testing.assert_array_equal(
        xz_slice,
        np.array(
            [
                [25., 26.],
                [19., 20.],
                [13., 14.],
                [7.,  8.],
                [1.,  2.]
            ]
        )
    )
