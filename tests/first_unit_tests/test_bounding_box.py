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
    assert test_bb_1.extent["xy"] == (-10.0, 1.0, -20.0, 2.0)
    assert test_bb_1.extent["xz"] == (-10.0, 1.0, -30.0, 3.0)
    assert test_bb_1.extent["yz"] == (-20.0, 2.0, -30.0, 3.0)


def test_bounding_box_methods():
    test_bb = openmc.BoundingBox.infinite()

    # check assignment operator
    test_bb[0] = [-10, -11, -12]
    test_bb[1] = [13, 14, 15]

    assert all(test_bb[0] == [-10, -11, -12])
    assert all(test_bb[1] == [13, 14, 15])

    # check length and iteration
    assert len(test_bb) == 2
    ll, ur = test_bb
    assert all(ll == [-10, -11, -12])
    assert all(ur == [13, 14, 15])

    # test expand/reduce methods
    other_bb = openmc.BoundingBox([-5, -5, -50], [5, 50, 5])

    reduced_bb = test_bb & other_bb

    # inplace was False by default. BoundingBox.reduce should return a new object
    assert test_bb is not reduced_bb

    # the original bounding box should be unchanged
    assert all(test_bb[0] == [-10, -11, -12])
    assert all(test_bb[1] == [13, 14, 15])

    assert all(reduced_bb[0] == [-5, -5, -12])
    assert all(reduced_bb[1] == [5, 14, 5])

    test_bb &= other_bb

    assert all(test_bb[0] == [-5, -5, -12])
    assert all(test_bb[1] == [5, 14, 5])

    other_bb = openmc.BoundingBox([-50, -50, -1], [50, 1, 50])

    expanded_bb = test_bb | other_bb

    # inplace was False by default. BoundingBox.expand should return a new object
    assert test_bb is not expanded_bb

    # the original bounding box should be unchanged
    assert all(test_bb[0] == [-5, -5, -12])
    assert all(test_bb[1] == [5, 14, 5])

    assert all(expanded_bb[0] == [-50, -50, -12])
    assert all(expanded_bb[1] == [50, 14, 50])

    test_bb |= other_bb

    assert all(test_bb[0] == [-50, -50, -12])
    assert all(test_bb[1] == [50, 14, 50])

    extended_bbox = test_bb.expand(0.1)

    assert extended_bbox is not test_bb

    # the original bounding box should not be changed with inplace as False
    assert all(test_bb[0] == [-50, -50, -12])
    assert all(test_bb[1] == [50, 14, 50])

    assert all(extended_bbox[0] == [-50.1, -50.1, -12.1])
    assert all(extended_bbox[1] == [50.1, 14.1, 50.1])

    extended_bbox = test_bb.expand(0.1, True)

    # inplace was set to True. BoundingBox.reduce should return the same object
    assert extended_bbox is test_bb

    assert all(test_bb[0] == [-50.1, -50.1, -12.1])
    assert all(test_bb[1] == [50.1, 14.1, 50.1])


@pytest.mark.parametrize(
    "bb, other, expected",
    [
        (test_bb_1, (0, 0, 0), True),
        (test_bb_2, (3, 3, 3), False),
        # completely disjoint
        (test_bb_1, test_bb_2, False),
        # contained but touching border
        (test_bb_1, test_bb_3, False),
        # Fully contained
        (test_bb_1, openmc.BoundingBox((-9, -19, -29), (0, 0, 0)), True),
        # intersecting boxes
        (test_bb_1, openmc.BoundingBox((-9, -19, -29), (1, 2, 5)), False),
    ],
)
def test_bounding_box_contains(bb, other, expected):
    assert (other in bb) == expected


@pytest.mark.parametrize(
    "invalid, ex",
    [
        ((1, 0), ValueError),
        ((1, 2, 3, 4), ValueError),
        ("foo", TypeError),
    ],
)
def test_bounding_box_contains_checking(invalid, ex):
    with pytest.raises(ex):
        invalid in test_bb_1
