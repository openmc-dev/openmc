from abc import ABCMeta
from collections.abc import Iterable
from numbers import Real, Integral
from xml.etree import ElementTree as ET
import sys
import warnings

import numpy as np

import openmc.checkvalue as cv
import openmc
from openmc._xml import get_text
from openmc.mixin import EqualityMixin, IDManagerMixin


class MeshBase(IDManagerMixin, metaclass=ABCMeta):
    """A mesh that partitions geometry for tallying purposes.

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

    """

    next_id = 1
    used_ids = set()

    def __init__(self, mesh_id=None, name=''):
        # Initialize Mesh class attributes
        self.id = mesh_id
        self.name = name

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        if name is not None:
            cv.check_type('name for mesh ID="{0}"'.format(self._id),
                          name, str)
            self._name = name
        else:
            self._name = ''

    def __repr__(self):
        string = type(self).__name__ + '\n'
        string += '{0: <16}{1}{2}\n'.format('\tID', '=\t', self._id)
        string += '{0: <16}{1}{2}\n'.format('\tName', '=\t', self._name)
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
        openmc.MeshBase
            Instance of a MeshBase subclass

        """

        mesh_type = group['type'][()].decode()
        if mesh_type == 'regular':
            return RegularMesh.from_hdf5(group)
        elif mesh_type == 'rectilinear':
            return RectilinearMesh.from_hdf5(group)
        elif mesh_type == 'unstructured':
            return UnstructuredMesh.from_hdf5(group)
        else:
            raise ValueError('Unrecognized mesh type: "' + mesh_type + '"')


class RegularMesh(MeshBase):
    """A regular Cartesian mesh in one, two, or three dimensions

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
    dimension : Iterable of int
        The number of mesh cells in each direction.
    n_dimension : int
        Number of mesh dimensions.
    lower_left : Iterable of float
        The lower-left corner of the structured mesh. If only two coordinate
        are given, it is assumed that the mesh is an x-y mesh.
    upper_right : Iterable of float
        The upper-right corner of the structrued mesh. If only two coordinate
        are given, it is assumed that the mesh is an x-y mesh.
    width : Iterable of float
        The width of mesh cells in each direction.
    indices : Iterable of tuple
        An iterable of mesh indices for each mesh element, e.g. [(1, 1, 1),
        (2, 1, 1), ...]

    """

    def __init__(self, mesh_id=None, name=''):
        super().__init__(mesh_id, name)

        self._dimension = None
        self._lower_left = None
        self._upper_right = None
        self._width = None

    @property
    def dimension(self):
        return self._dimension

    @property
    def n_dimension(self):
        if self._dimension is not None:
            return len(self._dimension)
        else:
            return None

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

    @property
    def indices(self):
        ndim = len(self._dimension)
        if ndim == 3:
            nx, ny, nz = self.dimension
            return ((x, y, z)
                    for z in range(1, nz + 1)
                    for y in range(1, ny + 1)
                    for x in range(1, nx + 1))
        elif ndim == 2:
            nx, ny = self.dimension
            return ((x, y)
                    for y in range(1, ny + 1)
                    for x in range(1, nx + 1))
        else:
            nx, = self.dimension
            return ((x,) for x in range(1, nx + 1))

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

    def __repr__(self):
        string = super().__repr__()
        string += '{0: <16}{1}{2}\n'.format('\tDimensions', '=\t', self.n_dimension)
        string += '{0: <16}{1}{2}\n'.format('\tMesh Cells', '=\t', self._dimension)
        string += '{0: <16}{1}{2}\n'.format('\tWidth', '=\t', self._lower_left)
        string += '{0: <16}{1}{2}\n'.format('\tOrigin', '=\t', self._upper_right)
        string += '{0: <16}{1}{2}\n'.format('\tPixels', '=\t', self._width)
        return string

    @classmethod
    def from_hdf5(cls, group):
        mesh_id = int(group.name.split('/')[-1].lstrip('mesh '))

        # Read and assign mesh properties
        mesh = cls(mesh_id)
        mesh.dimension = group['dimension'][()]
        mesh.lower_left = group['lower_left'][()]
        mesh.upper_right = group['upper_right'][()]
        mesh.width = group['width'][()]

        return mesh

    @classmethod
    def from_rect_lattice(cls, lattice, division=1, mesh_id=None, name=''):
        """Create mesh from an existing rectangular lattice

        Parameters
        ----------
        lattice : openmc.RectLattice
            Rectangular lattice used as a template for this mesh
        division : int
            Number of mesh cells per lattice cell.
            If not specified, there will be 1 mesh cell per lattice cell.
        mesh_id : int
            Unique identifier for the mesh
        name : str
            Name of the mesh

        Returns
        -------
        openmc.RegularMesh
            RegularMesh instance

        """
        cv.check_type('rectangular lattice', lattice, openmc.RectLattice)

        shape = np.array(lattice.shape)
        width = lattice.pitch*shape

        mesh = cls(mesh_id, name)
        mesh.lower_left = lattice.lower_left
        mesh.upper_right = lattice.lower_left + width
        mesh.dimension = shape*division

        return mesh

    def to_xml_element(self):
        """Return XML representation of the mesh

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing mesh data

        """

        element = ET.Element("mesh")
        element.set("id", str(self._id))

        if self._dimension is not None:
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

    @classmethod
    def from_xml_element(cls, elem):
        """Generate mesh from an XML element

        Parameters
        ----------
        elem : xml.etree.ElementTree.Element
            XML element

        Returns
        -------
        openmc.Mesh
            Mesh generated from XML element

        """
        mesh_id = int(get_text(elem, 'id'))
        mesh = cls(mesh_id)

        mesh_type = get_text(elem, 'type')
        if mesh_type is not None:
            mesh.type = mesh_type

        dimension = get_text(elem, 'dimension')
        if dimension is not None:
            mesh.dimension = [int(x) for x in dimension.split()]

        lower_left = get_text(elem, 'lower_left')
        if lower_left is not None:
            mesh.lower_left = [float(x) for x in lower_left.split()]

        upper_right = get_text(elem, 'upper_right')
        if upper_right is not None:
            mesh.upper_right = [float(x) for x in upper_right.split()]

        width = get_text(elem, 'width')
        if width is not None:
            mesh.width = [float(x) for x in width.split()]

        return mesh

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

        n_dim = len(self.dimension)

        # Build the cell which will contain the lattice
        xplanes = [openmc.XPlane(self.lower_left[0], bc[0]),
                   openmc.XPlane(self.upper_right[0], bc[1])]
        if n_dim == 1:
            yplanes = [openmc.YPlane(-1e10, 'reflective'),
                       openmc.YPlane(1e10, 'reflective')]
        else:
            yplanes = [openmc.YPlane(self.lower_left[1], bc[2]),
                       openmc.YPlane(self.upper_right[1], bc[3])]

        if n_dim <= 2:
            # Would prefer to have the z ranges be the max supported float, but
            # these values are apparently different between python and Fortran.
            # Choosing a safe and sane default.
            # Values of +/-1e10 are used here as there seems to be an
            # inconsistency between what numpy uses as the max float and what
            # Fortran expects for a real(8), so this avoids code complication
            # and achieves the same goal.
            zplanes = [openmc.ZPlane(-1e10, 'reflective'),
                       openmc.ZPlane(1e10, 'reflective')]
        else:
            zplanes = [openmc.ZPlane(self.lower_left[2], bc[4]),
                       openmc.ZPlane(self.upper_right[2], bc[5])]
        root_cell = openmc.Cell()
        root_cell.region = ((+xplanes[0] & -xplanes[1]) &
                            (+yplanes[0] & -yplanes[1]) &
                            (+zplanes[0] & -zplanes[1]))

        # Build the universes which will be used for each of the (i,j,k)
        # locations within the mesh.
        # We will concurrently build cells to assign to these universes
        cells = []
        universes = []
        for index in self.indices:
            cells.append(openmc.Cell())
            universes.append(openmc.Universe())
            universes[-1].add_cell(cells[-1])

        lattice = openmc.RectLattice()
        lattice.lower_left = self.lower_left

        # Assign the universe and rotate to match the indexing expected for
        # the lattice
        if n_dim == 1:
            universe_array = np.array([universes])
        elif n_dim == 2:
            universe_array = np.empty(self.dimension[::-1],
                                      dtype=openmc.Universe)
            i = 0
            for y in range(self.dimension[1] - 1, -1, -1):
                for x in range(self.dimension[0]):
                    universe_array[y][x] = universes[i]
                    i += 1
        else:
            universe_array = np.empty(self.dimension[::-1],
                                      dtype=openmc.Universe)
            i = 0
            for z in range(self.dimension[2]):
                for y in range(self.dimension[1] - 1, -1, -1):
                    for x in range(self.dimension[0]):
                        universe_array[z][y][x] = universes[i]
                        i += 1
        lattice.universes = universe_array

        if self.width is not None:
            lattice.pitch = self.width
        else:
            dx = ((self.upper_right[0] - self.lower_left[0]) /
                  self.dimension[0])

            if n_dim == 1:
                lattice.pitch = [dx]
            elif n_dim == 2:
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


def Mesh(*args, **kwargs):
    warnings.warn("Mesh has been renamed RegularMesh. Future versions of "
                  "OpenMC will not accept the name Mesh.")
    return RegularMesh(*args, **kwargs)


class RectilinearMesh(MeshBase):
    """A 3D rectilinear Cartesian mesh

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
    n_dimension : int
        Number of mesh dimensions (always 3 for a RectilinearMesh).
    x_grid : Iterable of float
        Mesh boundary points along the x-axis.
    y_grid : Iterable of float
        Mesh boundary points along the y-axis.
    z_grid : Iterable of float
        Mesh boundary points along the z-axis.
    indices : Iterable of tuple
        An iterable of mesh indices for each mesh element, e.g. [(1, 1, 1),
        (2, 1, 1), ...]

    """

    def __init__(self, mesh_id=None, name=''):
        super().__init__(mesh_id, name)

        self._x_grid = None
        self._y_grid = None
        self._z_grid = None

    @property
    def n_dimension(self):
        return 3

    @property
    def x_grid(self):
        return self._x_grid

    @property
    def y_grid(self):
        return self._y_grid

    @property
    def z_grid(self):
        return self._z_grid

    @property
    def indices(self):
        nx = len(self.x_grid) - 1
        ny = len(self.y_grid) - 1
        nz = len(self.z_grid) - 1
        return ((x, y, z)
                for z in range(1, nz + 1)
                for y in range(1, ny + 1)
                for x in range(1, nx + 1))

    @x_grid.setter
    def x_grid(self, grid):
        cv.check_type('mesh x_grid', grid, Iterable, Real)
        self._x_grid = grid

    @y_grid.setter
    def y_grid(self, grid):
        cv.check_type('mesh y_grid', grid, Iterable, Real)
        self._y_grid = grid

    @z_grid.setter
    def z_grid(self, grid):
        cv.check_type('mesh z_grid', grid, Iterable, Real)
        self._z_grid = grid

    def __repr__(self):
        fmt = '{0: <16}{1}{2}\n'
        string = super().__repr__()
        string += fmt.format('\tDimensions', '=\t', self.n_dimension)
        x_grid_str = str(self._x_grid) if not self._x_grid else len(self._x_grid)
        string += fmt.format('\tN X pnts:', '=\t', x_grid_str)
        if self._x_grid:
            string += fmt.format('\tX Min:', '=\t', self._x_grid[0])
            string += fmt.format('\tX Max:', '=\t', self._x_grid[-1])
        y_grid_str = str(self._y_grid) if not self._y_grid else len(self._y_grid)
        string += fmt.format('\tN Y pnts:', '=\t', y_grid_str)
        if self._y_grid:
            string += fmt.format('\tY Min:', '=\t', self._y_grid[0])
            string += fmt.format('\tY Max:', '=\t', self._y_grid[-1])
        z_grid_str = str(self._z_grid) if not self._z_grid else len(self._z_grid)
        string += fmt.format('\tN Z pnts:', '=\t', z_grid_str)
        if self._z_grid:
            string += fmt.format('\tZ Min:', '=\t', self._z_grid[0])
            string += fmt.format('\tZ Max:', '=\t', self._z_grid[-1])
        return string

    @classmethod
    def from_hdf5(cls, group):
        mesh_id = int(group.name.split('/')[-1].lstrip('mesh '))

        # Read and assign mesh properties
        mesh = cls(mesh_id)
        mesh.x_grid = group['x_grid'][()]
        mesh.y_grid = group['y_grid'][()]
        mesh.z_grid = group['z_grid'][()]

        return mesh

    def to_xml_element(self):
        """Return XML representation of the mesh

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing mesh data

        """

        element = ET.Element("mesh")
        element.set("id", str(self._id))
        element.set("type", "rectilinear")

        subelement = ET.SubElement(element, "x_grid")
        subelement.text = ' '.join(map(str, self.x_grid))

        subelement = ET.SubElement(element, "y_grid")
        subelement.text = ' '.join(map(str, self.y_grid))

        subelement = ET.SubElement(element, "z_grid")
        subelement.text = ' '.join(map(str, self.z_grid))

        return element


class UnstructuredMesh(MeshBase):
    """A 3D unstructured mesh

    Parameters
    ----------
    filename : str
        Location of the unstructured mesh file
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
    filename : str
        Name of the file containing the unstructured mesh
    volumes : Iterable of float
        Volumes of the unstructured mesh elements
    total_volume : float
        Volume of the unstructured mesh in total
    centroids : Iterable of tuple
        An iterable of element centroid coordinates, e.g. [(0.0, 0.0, 0.0),
        (1.0, 1.0, 1.0), ...]
    """

    def __init__(self, filename, mesh_id=None, name=''):
        super().__init__(mesh_id, name)
        self.filename = filename
        self._volumes = []
        self._centroids = []

    @property
    def filename(self):
        return self._filename

    @filename.setter
    def filename(self, filename):
        cv.check_type('Unstructured Mesh filename', filename, str)
        self._filename = filename

    @property
    def volumes(self):
        return self._volumes

    @volumes.setter
    def volumes(self, volumes):
        cv.check_type("Unstructured mesh volumes", volumes, Iterable, Real)
        self._volumes = volumes

    @property
    def total_volume(self):
        return np.sum(self.volumes)

    @property
    def centroids(self):
        return self._centroids

    @centroids.setter
    def centroids(self, centroids):
        cv.check_type("Unstructured mesh centroids", centroids,
                      Iterable, Real)
        self._centroids = centroids

    def __repr__(self):
        string = super().__repr__()
        return string + '{: <16}=\t{}\n'.format('\tFilename', self.filename)

    @classmethod
    def from_hdf5(cls, group):
        mesh_id = int(group.name.split('/')[-1].lstrip('mesh '))
        filename = group['filename'][()].decode()

        mesh = cls(filename, mesh_id=mesh_id)
        vol_data = group['volumes'][()]
        centroids = group['centroids'][()]
        mesh.volumes = np.reshape(vol_data, (vol_data.shape[0],))
        mesh.centroids = np.reshape(centroids, (vol_data.shape[0], 3))

        return mesh

    def to_xml_element(self):
        """Return XML representation of the mesh

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing mesh data

        """

        element = ET.Element("mesh")
        element.set("id", str(self._id))
        element.set("type", "unstructured")

        subelement = ET.SubElement(element, "filename")
        subelement.text = self.filename

        return element

    @classmethod
    def from_xml_element(cls, elem):
        """Generate unstructured mesh object from XML element

        Parameters
        ----------
        elem : xml.etree.ElementTree.Element
            XML element

        Returns
        -------
        openmc.UnstructuredMesh
            UnstructuredMesh generated from an XML element
        """
        mesh_id = int(get_text(elem, 'id'))
        filename = get_text(elem, 'filename')

        mesh = cls(filename, mesh_id)

        return mesh
