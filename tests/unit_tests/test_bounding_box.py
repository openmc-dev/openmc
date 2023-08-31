import numpy as np
import openmc
import pytest


test_bb_1 = openmc.BoundingBox((-10.0, -20.0, -30.0), (1.0, 2.0, 3.0))
test_bb_2 = openmc.BoundingBox((1.0, 2.0, 3.0), (11.0, 22.0, 33.0))
test_bb_3 = openmc.BoundingBox((-10.0, -20.0, -30.0), (-1.0, -2.0, -3.0))


@pytest.mark.parametrize(
    "bb, expected",
    [
        (test_bb_1, 7986),  # 11 * 22 * 33
        (test_bb_2, 6000),  # 10 * 20 * 30
        (test_bb_3, 4374),  # 9 * 18 * 27
    ],
)
def test_bounding_box_volume(bb, expected):
    assert bb.volume == expected


@pytest.mark.parametrize(
    "bb, expected",
    [
        (test_bb_1, np.array([-10.0, -20.0, -30.0])),
        (test_bb_2, np.array([1.0, 2.0, 3.0])),
        (test_bb_3, np.array([-10.0, -20.0, -30.0])),
    ],
)
def test_bounding_lower_left(bb, expected):
    assert np.array_equiv(expected, bb.lower_left)


@pytest.mark.parametrize(
    "bb, expected",
    [
        (test_bb_1, np.array([1.0, 2.0, 3.0])),
        (test_bb_2, np.array([11.0, 22.0, 33.0])),
        (test_bb_3, np.array([-1.0, -2.0, -3.0])),
    ],
)
def test_bounding_upper_right(bb, expected):
    assert np.array_equiv(expected, bb.upper_right)


@pytest.mark.parametrize(
    "bb, expected",
    [
        (test_bb_1, np.array([-4.5, -9.0, -13.5])),
        (test_bb_2, np.array([6.0, 12.0, 18.0])),
        (test_bb_3, np.array([-5.5, -11.0, -16.5])),
    ],
)
def test_bounding_box_center(bb, expected):
    assert np.array_equiv(expected, bb.center)


def test_bounding_box_input_checking():
    # checks that only passing lower_left is not accepted
    with pytest.raises(TypeError):
        openmc.BoundingBox((-10, -20, -3))
    # checks that a tuple with three entry is not accepted
    with pytest.raises(TypeError):
        openmc.BoundingBox((-1, -2, -3), (-1, -2, -3), (-1, -2, -3))
    # checks that a numpy array with two entries is not accepted
    with pytest.raises(ValueError):
        openmc.BoundingBox(np.array([-10, -30]), np.array([1, 2, 3]))
    # checks that a numpy array with two entries is not accepted
    with pytest.raises(ValueError):
        openmc.BoundingBox(np.array([-10, -20, -30]), np.array([1, 3]))
    # checks that a numpy array with four entries is not accepted
    with pytest.raises(ValueError):
        openmc.BoundingBox(np.array([-10, -20, -3, -4]), np.array([1, 2, 3]))
    # checks that a numpy array with four entries is not accepted
    with pytest.raises(ValueError):
        openmc.BoundingBox(np.array([-10, -20, -4]), np.array([1, 2, 3, 4]))


def test_bounding_box_extents():
    assert test_bb_1.extent['xy'] == (-10., 1., -20., 2.)
    assert test_bb_1.extent['xz'] == (-10., 1., -30., 3.)
    assert test_bb_1.extent['yz'] == (-20., 2., -30., 3.)
