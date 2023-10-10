from __future__ import annotations
from typing import Iterable

import numpy as np

from .checkvalue import check_length


class BoundingBox:
    """Axis-aligned bounding box.

    .. versionadded:: 0.13.4

    Parameters
    ----------
    lower_left : iterable of float
        The x, y, z coordinates of the lower left corner of the bounding box in [cm]
    upper_right : iterable of float
        The x, y, z coordinates of the upper right corner of the bounding box in [cm]

    Attributes
    ----------
    center : numpy.ndarray
        x, y, z coordinates of the center of the bounding box in [cm]
    lower_left : numpy.ndarray
        The x, y, z coordinates of the lower left corner of the bounding box in [cm]
    upper_right : numpy.ndarray
        The x, y, z coordinates of the upper right corner of the bounding box in [cm]
    volume : float
        The volume of the bounding box in [cm^3]
    extent : dict
        A dictionary of basis as keys and the extent (left, right, bottom, top)
        as values. Intended use in Matplotlib plots when setting extent
    width : iterable of float
        The width of the x, y and z axis in [cm]
    """

    def __init__(self, lower_left: Iterable[float], upper_right: Iterable[float]):
        check_length("lower_left", lower_left, 3, 3)
        check_length("upper_right", upper_right, 3, 3)
        self._bounds = np.asarray([lower_left, upper_right], dtype=float)

    def __repr__(self) -> str:
        return "BoundingBox(lower_left={}, upper_right={})".format(
            tuple(self.lower_left), tuple(self.upper_right))

    def __getitem__(self, key) -> np.ndarray:
        return self._bounds[key]

    def __len__(self):
        return 2

    def __setitem__(self, key, val):
        self._bounds[key] = val

    @property
    def center(self) -> np.ndarray:
        return (self[0] + self[1]) / 2

    @property
    def lower_left(self) -> np.ndarray:
        return self[0]

    @lower_left.setter
    def lower_left(self, llc):
        check_length('lower_left', llc, 3, 3)
        self[0] = llc

    @property
    def upper_right(self) -> np.ndarray:
        return self[1]

    @upper_right.setter
    def upper_right(self, urc):
        check_length('upper_right', urc, 3, 3)
        self[1] = urc

    @property
    def volume(self) -> float:
        return np.abs(np.prod(self[1] - self[0]))

    @property
    def extent(self):
        return {
            "xy": (
                self.lower_left[0],
                self.upper_right[0],
                self.lower_left[1],
                self.upper_right[1],
            ),
            "xz": (
                self.lower_left[0],
                self.upper_right[0],
                self.lower_left[2],
                self.upper_right[2],
            ),
            "yz": (
                self.lower_left[1],
                self.upper_right[1],
                self.lower_left[2],
                self.upper_right[2],
            ),
        }

    @property
    def width(self):
        return self.upper_right - self.lower_left

    def extend(self, padding_distance: float, inplace: bool = False) -> BoundingBox:
        """Returns an extended bounding box

        Parameters
        ----------
        padding_distance : float
            The distance to enlarge the bounding box by
        inplace : bool
            Whether or not to return a new BoundingBox instance or to modify the
            current BoundingBox object.

        Returns
        -------
        An extended bounding box
        """
        if inplace:
            self[0] -= padding_distance
            self[1] += padding_distance
            return self
        else:
            return BoundingBox(self[0] - padding_distance, self[1] + padding_distance)

    def expand(self, other_box: BoundingBox, inplace: bool = False) -> BoundingBox:
        """Expand the box to contain another box

        Parameters
        ----------
        other_box : BoundingBox
            The box used to resize this box
        inplace : bool
            Whether or not to return a new BoundingBox instance or to modify the
            current BoundingBox object.

        Returns
        -------
        An expanded bounding box
        """
        llc = np.minimum(self.lower_left, other_box.lower_left)
        urc = np.maximum(self.upper_right, other_box.upper_right)
        if inplace:
            self.lower_left = llc
            self.upper_right = urc
            return self
        else:
            return BoundingBox(llc, urc)

    def reduce(self, other_box: BoundingBox, inplace: bool = False) -> BoundingBox:
        """Reduce the box to match the dimensions of the input bounding box

        Parameters
        ----------
        other_box : BoundingBox
            The box used to resize this box
        inplace : bool
            Whether or not to return a new BoundingBox instance or to modify the
            current BoundingBox object.

        Returns
        -------
        A reduced bounding box
        """
        llc = np.maximum(self.lower_left, other_box.lower_left)
        urc = np.minimum(self.upper_right, other_box.upper_right)
        if inplace:
            self.lower_left = llc
            self.upper_right = urc
            return self
        else:
            return BoundingBox(llc, urc)

    @classmethod
    def infinite(cls):
        """Create an infinite box. Useful as a starting point for determining
           geometry bounds.

        Returns
        -------
        An infinitely large bounding box.
        """
        infs = np.full((3,), np.inf)
        return cls(-infs, infs)
