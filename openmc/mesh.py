from collections import Iterable
from numbers import Real, Integral
from xml.etree import ElementTree as ET
import sys

from six import string_types
import numpy as np

import openmc.checkvalue as cv
import openmc
from openmc.mixin import EqualityMixin, IDManagerMixin


class Mesh(EqualityMixin, IDManagerMixin):
    """A structured Cartesian mesh in one, two, or three dimensions

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

    next_id = 1
    used_ids = set()

    def __init__(self, mesh_id=None, name=''):
        # Initialize Mesh class attributes
        self.id = mesh_id
        self.name = name
        self._type = 'regular'
        self._dimension = None
        self._lower_left = None
        self._upper_right = None
        self._width = None

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

    @name.setter
    def name(self, name):
        if name is not None:
            cv.check_type('name for mesh ID="{0}"'.format(self._id),
                          name, string_types)
            self._name = name
        else:
            self._name = ''

    @type.setter
    def type(self, meshtype):
        cv.check_type('type for mesh ID="{0}"'.format(self._id),
                      meshtype, string_types)
        cv.check_value('type for mesh ID="{0}"'.format(self._id),
                       meshtype, ['regular'])
        self._type = meshtype

    @dimension.setter
    def dimension(self, dimension):
        cv.check_type('mesh dimension', dimension, Iterable, Integral)
        cv.check_length('mesh dimension', dimension, 1, 3)
        self._dimension = dimension

    @lower_left.setter
    def lower_left(self, lower_left):
        cv.check_type('mesh lower_left', lower_left, Iterable, Real)
        cv.check_length('mesh lower_left', lower_left, 1, 3)
        self._lower_left = lower_left

    @upper_right.setter
    def upper_right(self, upper_right):
        cv.check_type('mesh upper_right', upper_right, Iterable, Real)
        cv.check_length('mesh upper_right', upper_right, 1, 3)
        self._upper_right = upper_right

    @width.setter
    def width(self, width):
        cv.check_type('mesh width', width, Iterable, Real)
        cv.check_length('mesh width', width, 1, 3)
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

    @classmethod
    def from_hdf5(cls, group):
        """Create mesh from HDF5 group

        Parameters
        ----------
        group : h5py.Group
            Group in HDF5 file

        Returns
        -------
        openmc.Mesh
            Mesh instance

        """
        mesh_id = int(group.name.split('/')[-1].lstrip('mesh '))

        # Read and assign mesh properties
        mesh = cls(mesh_id)
        mesh.type = group['type'].value.decode()
        mesh.dimension = group['dimension'].value
        mesh.lower_left = group['lower_left'].value
        mesh.upper_right = group['upper_right'].value
        mesh.width = group['width'].value

        return mesh

    def cell_generator(self):
        """Generator function to traverse through every [i,j,k] index of the
        mesh

        For example the following code:

        .. code-block:: python

            for mesh_index in mymesh.cell_generator():
                print(mesh_index)

        will produce the following output for a 3-D 2x2x2 mesh in mymesh::

            [1, 1, 1]
            [2, 1, 1]
            [1, 2, 1]
            [2, 2, 1]
            ...


        """

        if len(self.dimension) == 1:
            for x in range(self.dimension[0]):
                    yield [x + 1, 1, 1]
        elif len(self.dimension) == 2:
            for y in range(self.dimension[1]):
                for x in range(self.dimension[0]):
                    yield [x + 1, y + 1, 1]
        else:
            for z in range(self.dimension[2]):
                for y in range(self.dimension[1]):
                    for x in range(self.dimension[0]):
                        yield [x + 1, y + 1, z + 1]

    def to_xml_element(self):
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

    def build_cells(self, bc=['reflective'] * 6):
        """Generates a lattice of universes with the same dimensionality
        as the mesh object.  The individual cells/universes produced
        will not have material definitions applied and so downstream code
        will have to apply that information.

        Parameters
        ----------
        bc : iterable of {'reflective', 'periodic', 'transmission', or 'vacuum'}
            Boundary conditions for each of the four faces of a rectangle
            (if aplying to a 2D mesh) or six faces of a parallelepiped
            (if applying to a 3D mesh) provided in the following order:
            [x min, x max, y min, y max, z min, z max].  2-D cells do not
            contain the z min and z max entries.

        Returns
        -------
        root_cell : openmc.Cell
            The cell containing the lattice representing the mesh geometry;
            this cell is a single parallelepiped with boundaries matching
            the outermost mesh boundary with the boundary conditions from bc
            applied.
        cells : iterable of openmc.Cell
            The list of cells within each lattice position mimicking the mesh
            geometry.

        """

        cv.check_length('bc', bc, length_min=4, length_max=6)
        for entry in bc:
            cv.check_value('bc', entry, ['transmission', 'vacuum',
                                         'reflective', 'periodic'])

        # Build the cell which will contain the lattice
        xplanes = [openmc.XPlane(x0=self.lower_left[0],
                                 boundary_type=bc[0]),
                   openmc.XPlane(x0=self.upper_right[0],
                                 boundary_type=bc[1])]
        if len(self.dimension) == 1:
            yplanes = [openmc.YPlane(y0=-1e10, boundary_type='reflective'),
                       openmc.YPlane(y0=1e10, boundary_type='reflective')]
        else:
            yplanes = [openmc.YPlane(y0=self.lower_left[1],
                                     boundary_type=bc[2]),
                       openmc.YPlane(y0=self.upper_right[1],
                                     boundary_type=bc[3])]

        if len(self.dimension) <= 2:
            # Would prefer to have the z ranges be the max supported float, but
            # these values are apparently different between python and Fortran.
            # Choosing a safe and sane default.
            # Values of +/-1e10 are used here as there seems to be an
            # inconsistency between what numpy uses as the max float and what
            # Fortran expects for a real(8), so this avoids code complication
            # and achieves the same goal.
            zplanes = [openmc.ZPlane(z0=-1e10, boundary_type='reflective'),
                       openmc.ZPlane(z0=1e10, boundary_type='reflective')]
        else:
            zplanes = [openmc.ZPlane(z0=self.lower_left[2],
                                     boundary_type=bc[4]),
                       openmc.ZPlane(z0=self.upper_right[2],
                                     boundary_type=bc[5])]
        root_cell = openmc.Cell()
        root_cell.region = ((+xplanes[0] & -xplanes[1]) &
                            (+yplanes[0] & -yplanes[1]) &
                            (+zplanes[0] & -zplanes[1]))

        # Build the universes which will be used for each of the [i,j,k]
        # locations within the mesh.
        # We will concurrently build cells to assign to these universes
        cells = []
        universes = []
        for [i, j, k] in self.cell_generator():
            cells.append(openmc.Cell())
            universes.append(openmc.Universe())
            universes[-1].add_cell(cells[-1])

        lattice = openmc.RectLattice()
        lattice.lower_left = self.lower_left

        # Assign the universe and rotate to match the indexing expected for
        # the lattice
        lattice.universes = np.rot90(np.reshape(universes, self.dimension))

        if self.width is not None:
            lattice.pitch = self.width
        else:
            dx = ((self.upper_right[0] - self.lower_left[0]) /
                  self.dimension[0])

            if len(self.dimension) == 1:
                lattice.pitch = [dx]
            elif len(self.dimension) == 2:
                dy = ((self.upper_right[1] - self.lower_left[1]) /
                      self.dimension[1])
                lattice.pitch = [dx, dy]
            else:
                dy = ((self.upper_right[1] - self.lower_left[1]) /
                      self.dimension[1])
                dz = ((self.upper_right[2] - self.lower_left[2]) /
                      self.dimension[2])
                lattice.pitch = [dx, dy, dz]

        # Fill Cell with the Lattice
        root_cell.fill = lattice

        return root_cell, cells
