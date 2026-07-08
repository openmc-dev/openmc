from collections.abc import Iterable
from numbers import Integral
from functools import wraps

import h5py
import lxml.etree as ET

import openmc
import openmc.checkvalue as cv
from ._xml import get_text
from .bounding_box import BoundingBox
from .checkvalue import PathLike
from .utility_funcs import input_path


def require_statepoint_data(func):
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        if not self._has_statepoint_data:
            raise AttributeError(f'The "{func.__name__}" property requires '
                                 'information about this mesh to be loaded '
                                 'from a statepoint file.')
        return func(self, *args, **kwargs)
    return wrapper


class XDGMesh(openmc.MeshBase):
    """A 3D unstructured mesh supported by the XDG library. Various mesh
    backends are supported in XDG, including MOAB and libMesh. For more
    information on supported mesh types and libraries, see
    https://xdg-org.github.io/xdg/.

    This class is meant to serve as a reference to an unstructured mesh to be
    applied in an OpenMC model in one of the following ways:

    1. As a volumetric mesh for an tally MeshFilter (see
       :class:`openmc.MeshFilter`)
    2. As a surface mesh for an XDGUniverse (see :class:`openmc.XDGUniverse`)
    3. As a volumetric mesh for an XDGUniverse (see :class:`openmc.XDGUniverse`)


    Parameters
    ----------
    filename : path-like
        Location of the unstructured mesh file. Supported files for 'moab'
        library are .h5, .vtk and .e (exodus) depending on the MOAB build
        configuration. Supported files for 'libmesh' library are exodus mesh
        files .exo.
    library : {'moab', 'libmesh'}
        Mesh library used for the unstructured mesh tally
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
    library : {'moab', 'libmesh'}
        Mesh library used for the unstructured mesh tally
    """
    def __init__(self, filename: PathLike, library: str | None = None, mesh_id: int | None = None,
                 name: str = ''):
        super().__init__(mesh_id, name)
        self.filename = filename
        self._has_statepoint_data = False

        # set library based on file extension if not provided by the user
        if library is None:
            if filename.suffix in ('.h5m', '.vtk', '.e'):
                library = 'moab'
            elif filename.suffix in ('.exo',):
                library = 'libmesh'

        if library is None:
            raise ValueError(f"Unable to determine mesh library from file extension {filename.suffix}. "
                             "Please specify the mesh library type explicitly.")
        self.library = library

    @property
    def filename(self):
        return self._filename

    @filename.setter
    def filename(self, filename):
        cv.check_type('Unstructured Mesh filename', filename, PathLike)
        self._filename = input_path(filename)

    @property
    def library(self):
        return self._library

    @library.setter
    def library(self, lib: str | None):
        if lib is not None:
            cv.check_value('Unstructured mesh library', lib, ('moab', 'libmesh'))
        self._library = lib

    @property
    def dimension(self):
        return (self.n_elements,)

    @property
    def n_dimension(self):
        return 3

    @property
    def lower_left(self):
        return self.vertices.min(axis=0)

    @property
    def upper_right(self):
        self.vertices.max(axis=0)

    @property
    @require_statepoint_data
    def n_elements(self):
        return self._n_elements

    @n_elements.setter
    def n_elements(self, n):
        cv.check_type("XDGMesh n_elements", n, Integral)
        cv.check_greater_than("XDGMesh n_elements", n, 0, equality=False)
        self._n_elements = n

    @property
    def indices(self):
        return [(i,) for i in range(self.n_elements)]

    def __repr__(self):
        string = super().__repr__()
        string += '{: <16}=\t{}\n'.format('\tFilename', self.filename)
        string += '{: <16}=\t{}\n'.format('\tMesh Library', self.library)
        return string

    @property
    @require_statepoint_data
    def vertices(self):
        return self._vertices

    @property
    @require_statepoint_data
    def connectivity(self):
        return self._connectivity

    @property
    @require_statepoint_data
    def element_types(self):
        return self._element_types

    @property
    @require_statepoint_data
    def centroids(self):
        return np.array([self.centroid(i) for i in range(self.n_elements)])

    @require_statepoint_data
    def centroid(self, bin: int):
        """Return the vertex averaged centroid of an element

        Parameters
        ----------
        bin : int
            Bin ID for the returned centroid

        Returns
        -------
        numpy.ndarray
            x, y, z values of the element centroid

        """
        conn = self.connectivity[bin]
        # remove invalid connectivity values
        conn = conn[conn >= 0]
        coords = self.vertices[conn]
        return coords.mean(axis=0)

    @property
    def has_statepoint_data(self):
        return self._has_statepoint_data

    @has_statepoint_data.setter
    def has_statepoint_data(self, val):
        cv.check_type("XDGMesh has_statepoint_data", val, bool)
        self._has_statepoint_data = val

    @classmethod
    def from_hdf5(cls, group: h5py.Group, mesh_id: int, name: str):
        filename = group['filename'][()].decode()
        library = group['library'][()].decode()

        mesh = cls(filename=filename,
                  library=library,
                  mesh_id=mesh_id,
                  name=name)

        openmc.UnstructuredMesh._read_hdf5_mesh_data(group, mesh)
        mesh.has_statepoint_data = True

        return mesh

    def to_xml_element(self):
        """Return XML representation of the mesh

        Returns
        -------
        element : lxml.etree._Element
            XML element containing mesh data

        """
        element = super().to_xml_element()
        element.set("type", "xdg")
        element.set("library", self._library)
        subelement = ET.SubElement(element, "filename")
        subelement.text = str(self.filename)

        return element

    @classmethod
    def from_xml_element(cls, elem: ET.Element):
        """Generate unstructured mesh object from XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element

        Returns
        -------
        openmc.UnstructuredMesh
            UnstructuredMesh generated from an XML element
        """
        mesh_id = int(get_text(elem, 'id'))
        filename = get_text(elem, 'filename')
        library = get_text(elem, 'library')
        return cls(filename, library, mesh_id)
