from collections import Iterable
from numbers import Real, Integral
from xml.etree import ElementTree as ET
import sys

import numpy as np

import openmc.checkvalue as cv


if sys.version_info[0] >= 3:
    basestring = str

# "Static" variable for auto-generated and Mesh IDs
AUTO_MESH_ID = 10000


def reset_auto_mesh_id():
    """Reset counter for auto-generated mesh IDs."""
    global AUTO_MESH_ID
    AUTO_MESH_ID = 10000


class Mesh(object):
    """A structured Cartesian mesh in two or three dimensions

    Parameters
    ----------
    mesh_id : int
        Unique identifier for the mesh
    name : str
        Name of the mesh

    Attributes
    ----------
    id : int
        Unique identifier for the mesh
    name : str
        Name of the mesh
    type : str
        Type of the mesh
    dimension : Iterable of int
        The number of mesh cells in each direction.
    lower_left : Iterable of float
        The lower-left corner of the structured mesh. If only two coordinate are
        given, it is assumed that the mesh is an x-y mesh.
    upper_right : Iterable of float
        The upper-right corner of the structrued mesh. If only two coordinate
        are given, it is assumed that the mesh is an x-y mesh.
    width : Iterable of float
        The width of mesh cells in each direction.

    """

    def __init__(self, mesh_id=None, name=''):
        # Initialize Mesh class attributes
        self.id = mesh_id
        self.name = name
        self._type = 'regular'
        self._dimension = None
        self._lower_left = None
        self._upper_right = None
        self._width = None

    def __eq__(self, mesh2):
        # Check type
        if self._type != mesh2.type:
            return False

        # Check dimension
        elif self._dimension != mesh2.dimension:
            return False

        # Check width
        elif self._width != mesh2.width:
            return False

        # Check lower left / upper right
        elif self._lower_left != mesh2.lower_left and \
             self._upper_right != mesh2.upper_right:
            return False

        else:
            return True

    @property
    def id(self):
        return self._id

    @property
    def name(self):
        return self._name

    @property
    def type(self):
        return self._type

    @property
    def dimension(self):
        return self._dimension

    @property
    def lower_left(self):
        return self._lower_left

    @property
    def upper_right(self):
        return self._upper_right

    @property
    def width(self):
        return self._width

    @property
    def num_mesh_cells(self):
        return np.prod(self._dimension)

    @id.setter
    def id(self, mesh_id):
        if mesh_id is None:
            global AUTO_MESH_ID
            self._id = AUTO_MESH_ID
            AUTO_MESH_ID += 1
        else:
            cv.check_type('mesh ID', mesh_id, Integral)
            cv.check_greater_than('mesh ID', mesh_id, 0, equality=True)
            self._id = mesh_id

    @name.setter
    def name(self, name):
        if name is not None:
            cv.check_type('name for mesh ID="{0}"'.format(self._id),
                          name, basestring)
            self._name = name
        else:
            self._name = ''

    @type.setter
    def type(self, meshtype):
        cv.check_type('type for mesh ID="{0}"'.format(self._id),
                      meshtype, basestring)
        cv.check_value('type for mesh ID="{0}"'.format(self._id),
                       meshtype, ['regular'])
        self._type = meshtype

    @dimension.setter
    def dimension(self, dimension):
        cv.check_type('mesh dimension', dimension, Iterable, Integral)
        cv.check_length('mesh dimension', dimension, 2, 3)
        self._dimension = dimension

    @lower_left.setter
    def lower_left(self, lower_left):
        cv.check_type('mesh lower_left', lower_left, Iterable, Real)
        cv.check_length('mesh lower_left', lower_left, 2, 3)
        self._lower_left = lower_left

    @upper_right.setter
    def upper_right(self, upper_right):
        cv.check_type('mesh upper_right', upper_right, Iterable, Real)
        cv.check_length('mesh upper_right', upper_right, 2, 3)
        self._upper_right = upper_right

    @width.setter
    def width(self, width):
        cv.check_type('mesh width', width, Iterable, Real)
        cv.check_length('mesh width', width, 2, 3)
        self._width = width

    def __hash__(self):
        return hash(repr(self))

    def __repr__(self):
        string = 'Mesh\n'
        string += '{0: <16}{1}{2}\n'.format('\tID', '=\t', self._id)
        string += '{0: <16}{1}{2}\n'.format('\tName', '=\t', self._name)
        string += '{0: <16}{1}{2}\n'.format('\tType', '=\t', self._type)
        string += '{0: <16}{1}{2}\n'.format('\tBasis', '=\t', self._dimension)
        string += '{0: <16}{1}{2}\n'.format('\tWidth', '=\t', self._lower_left)
        string += '{0: <16}{1}{2}\n'.format('\tOrigin', '=\t', self._upper_right)
        string += '{0: <16}{1}{2}\n'.format('\tPixels', '=\t', self._width)
        return string

    def get_mesh_xml(self):
        """Return XML representation of the mesh

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing mesh data

        """

        element = ET.Element("mesh")
        element.set("id", str(self._id))
        element.set("type", self._type)

        subelement = ET.SubElement(element, "dimension")
        subelement.text = ' '.join(map(str, self._dimension))

        subelement = ET.SubElement(element, "lower_left")
        subelement.text = ' '.join(map(str, self._lower_left))

        if self._upper_right is not None:
            subelement = ET.SubElement(element, "upper_right")
            subelement.text = ' '.join(map(str, self._upper_right))

        if self._width is not None:
            subelement = ET.SubElement(element, "width")
            subelement.text = ' '.join(map(str, self._width))

        return element
