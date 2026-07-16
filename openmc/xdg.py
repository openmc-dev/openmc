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


class XDGMesh(openmc.UnstructuredMesh):
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
        super().__init__(mesh_id=mesh_id,
                         name=name,
                         filename=filename,
                         library=library)

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
