import copy
from xml.etree import ElementTree as ET

from openmc.checkvalue import *


# "Static" variable for auto-generated and Mesh IDs
AUTO_MESH_ID = 10000


def reset_auto_mesh_id():
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
    dimension : tuple or list or ndarray
        The number of mesh cells in each direction.
    lower_left : tuple or list or ndarray
        The lower-left corner of the structured mesh. If only two coordinate are
        given, it is assumed that the mesh is an x-y mesh.
    upper_right : tuple or list or ndarray
        The upper-right corner of the structrued mesh. If only two coordinate
        are given, it is assumed that the mesh is an x-y mesh.
    width : tuple or list or ndarray
        The width of mesh cells in each direction.

    """

    def __init__(self, mesh_id=None, name=''):
        # Initialize Mesh class attributes
        self.id = mesh_id
        self.name = name
        self._type = 'rectangular'
        self._dimension = None
        self._lower_left = None
        self._upper_right = None
        self._width = None

    def __eq__(self, mesh2):
        # Check type
        if self._type != mesh2._type:
            return False

        # Check dimension
        elif self._dimension != mesh2._dimension:
            return False

        # Check width
        elif self._width != mesh2._width:
            return False

        # Check lower left / upper right
        elif self._lower_left != mesh2._lower_left and \
             self._upper_right != mesh2._upper_right:
            return False

        else:
            return True

    def __deepcopy__(self, memo):
        existing = memo.get(id(self))

        # If this is the first time we have tried to copy this object, create a copy
        if existing is None:
            clone = type(self).__new__(type(self))
            clone._id = self._id
            clone._name = self._name
            clone._type = self._type
            clone._dimension = copy.deepcopy(self._dimension, memo)
            clone._lower_left = copy.deepcopy(self._lower_left, memo)
            clone._upper_right = copy.deepcopy(self._upper_right, memo)
            clone._width = copy.deepcopy(self._width, memo)

            memo[id(self)] = clone

            return clone

        # If this object has been copied before, return the first copy made
        else:
            return existing

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

        # Check that the ID is an integer and wasn't already used
        elif not is_integer(mesh_id):
            msg = 'Unable to set a non-integer Mesh ID {0}'.format(mesh_id)
            raise ValueError(msg)

        elif mesh_id < 0:
            msg = 'Unable to set Mesh ID to {0} since it must be a ' \
                  'non-negative integer'.format(mesh_id)
            raise ValueError(msg)

        else:
            self._id = mesh_id

    @name.setter
    def name(self, name):
        if not is_string(name):
            msg = 'Unable to set name for Mesh ID={0} with a non-string ' \
                  'value {1}'.format(self._id, name)
            raise ValueError(msg)
        else:
            self._name = name

    @type.setter
    def type(self, type):
        if not is_string(type):
            msg = 'Unable to set Mesh ID={0} for type {1} which is not ' \
                  'a string'.format(self._id, type)
            raise ValueError(msg)
        elif type not in ['rectangular', 'hexagonal']:
            msg = 'Unable to set Mesh ID={0} for type {1} which since ' \
                  'only rectangular and hexagonal meshes are ' \
                  'supported '.format(self._id, type)
            raise ValueError(msg)

        self._type = type

    @dimension.setter
    def dimension(self, dimension):
        if not isinstance(dimension, (tuple, list, np.ndarray)):
            msg = 'Unable to set Mesh ID={0} with dimension {1} which is ' \
                  'not a Python list, tuple or NumPy ' \
                  'array'.format(self._id, dimension)
            raise ValueError(msg)
        elif len(dimension) != 2 and len(dimension) != 3:
            msg = 'Unable to set Mesh ID={0} with dimension {1} since it ' \
                  'must include 2 or 3 dimensions'.format(self._id, dimension)
            raise ValueError(msg)

        for dim in dimension:
            if not is_integer(dim):
                msg = 'Unable to set Mesh ID={0} with dimension {1} which ' \
                      'is a non-integer'.format(self._id, dim)
                raise ValueError(msg)

        self._dimension = dimension

    @lower_left.setter
    def lower_left(self, lower_left):
        if not isinstance(lower_left, (tuple, list, np.ndarray)):
            msg = 'Unable to set Mesh ID={0} with lower_left {1} which is ' \
                  'not a Python list, tuple or NumPy ' \
                  'array'.format(self._id, lower_left)
            raise ValueError(msg)

        elif len(lower_left) != 2 and len(lower_left) != 3:
            msg = 'Unable to set Mesh ID={0} with lower_left {1} since it ' \
                   'must include 2 or 3 dimensions'.format(self._id, lower_left)
            raise ValueError(msg)

        for coord in lower_left:
            if not is_integer(coord) and not is_float(coord):
                msg = 'Unable to set Mesh ID={0} with lower_left {1} which ' \
                      'is neither neither an integer nor a floating point ' \
                      'value'.format(self._id, coord)
                raise ValueError(msg)

        self._lower_left = lower_left

    @upper_right.setter
    def upper_right(self, upper_right):
        if not isinstance(upper_right, (tuple, list, np.ndarray)):
            msg = 'Unable to set Mesh ID={0} with upper_right {1} which ' \
                  'is not a Python list, tuple or NumPy ' \
                  'array'.format(self._id, upper_right)
            raise ValueError(msg)

        if len(upper_right) != 2 and len(upper_right) != 3:
            msg = 'Unable to set Mesh ID={0} with upper_right {1} since it ' \
                  'must include 2 or 3 dimensions'.format(self._id, upper_right)
            raise ValueError(msg)

        for coord in upper_right:
            if not is_integer(coord) and not is_float(coord):
                msg = 'Unable to set Mesh ID={0} with upper_right {1} which ' \
                      'is neither an integer nor a floating point ' \
                      'value'.format(self._id, coord)
                raise ValueError(msg)

        self._upper_right = upper_right

    @width.setter
    def width(self, width):
        if width is not None:
            if not isinstance(width, (tuple, list, np.ndarray)):
                msg = 'Unable to set Mesh ID={0} with width {1} which ' \
                      'is not a Python list, tuple or NumPy ' \
                      'array'.format(self._id, width)
                raise ValueError(msg)

        if len(width) != 2 and len(width) != 3:
            msg = 'Unable to set Mesh ID={0} with width {1} since it must ' \
                  'include 2 or 3 dimensions'.format(self._id, width)
            raise ValueError(msg)

        for dim in width:
            if not is_integer(dim) and not is_float(dim):
                msg = 'Unable to set Mesh ID={0} with width {1} which is ' \
                      'neither an integer nor a floating point ' \
                      'value'.format(self._id, width)
                raise ValueError(msg)

        self._width = width

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
