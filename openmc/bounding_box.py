from __future__ import annotations
from collections.abc import Iterable

import numpy as np

from .checkvalue import check_length


class BoundingBox:
    """Axis-aligned bounding box.

    .. versionadded:: 0.14.0

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
            tuple(float(x) for x in self.lower_left),
            tuple(float(x) for x in self.upper_right),
        )

    def __getitem__(self, key) -> np.ndarray:
        return self._bounds[key]

    def __len__(self):
        return 2

    def __setitem__(self, key, val):
        self._bounds[key] = val

    def __iand__(self, other: BoundingBox) -> BoundingBox:
        """Updates the box be the intersection of itself and another box

        Parameters
        ----------
        other : BoundingBox
            The box used to resize this box

        Returns
        -------
        An updated bounding box
        """
        self.lower_left = np.maximum(self.lower_left, other.lower_left)
        self.upper_right = np.minimum(self.upper_right, other.upper_right)
        return self

    def __and__(self, other: BoundingBox) -> BoundingBox:
        new = BoundingBox(*self)
        new &= other
        return new

    def __ior__(self, other: BoundingBox) -> BoundingBox:
        """Updates the box be the union of itself and another box

        Parameters
        ----------
        other : BoundingBox
            The box used to resize this box

        Returns
        -------
        An updated bounding box
        """
        self.lower_left = np.minimum(self.lower_left, other.lower_left)
        self.upper_right = np.maximum(self.upper_right, other.upper_right)
        return self

    def __or__(self, other: BoundingBox) -> BoundingBox:
        new = BoundingBox(*self)
        new |= other
        return new

    def __contains__(self, other):
        """Check whether or not a point or another bounding box is in the bounding box.

        For another bounding box to be in the parent it must lie fully inside of it.
        """
        # test for a single point
        if isinstance(other, (tuple, list, np.ndarray)):
            point = other
            check_length("Point", point, 3, 3)
            return all(point > self.lower_left) and all(point < self.upper_right)
        elif isinstance(other, BoundingBox):
            return all([p in self for p in [other.lower_left, other.upper_right]])
        else:
            raise TypeError(
                f"Unable to determine if {other} is in the bounding box."
                f" Expected a tuple or a bounding box, but {type(other)} given"
            )

    @property
    def center(self) -> np.ndarray:
        return (self[0] + self[1]) / 2

    @property
    def lower_left(self) -> np.ndarray:
        return self[0]

    @lower_left.setter
    def lower_left(self, llc):
        check_length("lower_left", llc, 3, 3)
        self[0] = llc

    @property
    def upper_right(self) -> np.ndarray:
        return self[1]

    @upper_right.setter
    def upper_right(self, urc):
        check_length("upper_right", urc, 3, 3)
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

    def expand(self, padding_distance: float, inplace: bool = False) -> BoundingBox:
        """Returns an expanded bounding box

        Parameters
        ----------
        padding_distance : float
            The distance to enlarge the bounding box by
        inplace : bool
            Whether or not to return a new BoundingBox instance or to modify the
            current BoundingBox object.

        Returns
        -------
        An expanded bounding box
        """
        if inplace:
            self[0] -= padding_distance
            self[1] += padding_distance
            return self
        else:
            return BoundingBox(self[0] - padding_distance, self[1] + padding_distance)

    @classmethod
    def infinite(cls) -> BoundingBox:
        """Create an infinite box. Useful as a starting point for determining
           geometry bounds.

        Returns
        -------
        An infinitely large bounding box.
        """
        infs = np.full((3,), np.inf)
        return cls(-infs, infs)
