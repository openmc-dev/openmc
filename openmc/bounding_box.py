from __future__ import annotations
from typing import Iterable

import numpy as np

from .checkvalue import check_length


class BoundingBox(tuple):
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

    def __new__(cls, lower_left: Iterable[float], upper_right: Iterable[float]):
        check_length("lower_left", lower_left, 3, 3)
        check_length("upper_right", upper_right, 3, 3)
        lower_left = np.array(lower_left, dtype=float)
        upper_right = np.array(upper_right, dtype=float)
        return tuple.__new__(cls, (lower_left, upper_right))

    def __repr__(self) -> str:
        return "BoundingBox(lower_left={}, upper_right={})".format(
            tuple(self.lower_left), tuple(self.upper_right))

    @property
    def center(self) -> np.ndarray:
        return (self[0] + self[1]) / 2

    @property
    def lower_left(self) -> np.ndarray:
        return self[0]

    @property
    def upper_right(self) -> np.ndarray:
        return self[1]

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

    def extend(self, padding_distance: float) -> BoundingBox:
        """Returns an extended bounding box

        Parameters
        ----------
        padding_distance : float
            The distance to enlarge the bounding box by

        Returns
        -------
        An enlarged bounding box

        """
        return BoundingBox(np.array(
            [
                self[0][0] - padding_distance,
                self[0][1] - padding_distance,
                self[0][2] - padding_distance
            ]
        ), np.array(
            [
                self[1][0] + padding_distance,
                self[1][1] + padding_distance,
                self[1][2] + padding_distance
            ]
        ))
