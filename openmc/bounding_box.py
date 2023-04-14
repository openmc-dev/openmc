import numpy as np
from .checkvalue import check_type, check_length


class BoundingBox(tuple):
    """Axis-aligned bounding box.

    Parameters
    ----------
    corners : 2-tuple of numpy.array
        Lower-left and upper-right coordinates of an axis-aligned bounding box
        of the domain.

    Attributes
    ----------
    center: numpy.array
        x, y, z coordinates of the center of the bounding box in cm.
    lower_left: numpy.array
        The x, y, z coordinates of the lower left corner of the bounding box
    upper_right
        The x, y, z coordinates of the upper right corner of the bounding box
    volume: float
        The volume of the bounding box in cm3
    """

    def __init__(self, corners):

        check_type("corners", corners, tuple)
        check_type("corners", corners[0], np.ndarray)
        check_type("corners", corners[1], np.ndarray)
        check_length("corners", corners, 2, 2)
        check_length("corners", corners[0], 3, 3)
        check_length("corners", corners[1], 3, 3)
        self.corners = corners

    @property
    def center(self):
        """The center x, y, z coordinates of the bounding box"""
        return (self.corners[0] + self.corners[1]) / 2

    @property
    def lower_left(self):
        """The x, y, z coordinates of the lower left corner of the bounding box"""
        return self.corners[0]

    @property
    def upper_right(self):
        """The x, y, z coordinates of the upper right corner of the bounding box"""
        return self.corners[1]

    @property
    def volume(self):
        """The volume of the bounding box"""
        return np.abs(np.prod(self.corners[1] - self.corners[0]))
