from collections.abc import Iterable


import h5py
import lxml.etree as ET

import openmc
import openmc.checkvalue as cv
from ._xml import get_text
from .bounding_box import BoundingBox
from .checkvalue import PathLike
from .utility_funcs import input_path


class XDGMesh(openmc.MeshBase):
    """A 3D unstructured mesh

    Parameters
    ----------
    filename : path-like
        Location of the unstructured mesh file. Supported files for 'moab'
        library are .h5 and .vtk. Supported files for 'libmesh' library are
        exodus mesh files .exo.
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
    def __init__(self, filename: PathLike, library: str, mesh_id: int | None = None,
                 name: str = ''):
        super().__init__(mesh_id, name)
        self.filename = filename
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
    def library(self, lib: str):
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
        raise NotImplementedError("XDGMesh.lower_left is not implemented")

    @property
    def upper_right(self):
        raise NotImplementedError("XDGMesh.upper_right is not implemented")

    @property
    def n_elements(self):
        raise NotImplementedError("XDGMesh.n_elements is not implemented")

    @property
    def indices(self):
        raise NotImplementedError("XDGMesh.indices is not implemented")


    def __repr__(self):
        string = super().__repr__()
        string += '{: <16}=\t{}\n'.format('\tFilename', self.filename)
        string += '{: <16}=\t{}\n'.format('\tMesh Library', self.library)
        return string

    @classmethod
    def from_hdf5(cls, group: h5py.Group, mesh_id: int, name: str):
        raise NotImplementedError("XDGMesh.from_hdf5 is not implemented")
        # TODO: add length_multiplier

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


class XDGUniverse(openmc.UniverseBase):
    """A reference to a XDG file to be used in the model.

    Parameters
    ----------
    mesh : openmc.XDGMesh
        Mesh to use for the XDG universe
    type : str
        Type of XDG file to use. Options are 'volume_mesh' or 'surface_mesh'.
    universe_id : int, optional
        Unique identifier of the universe. If not specified, an identifier will
        automatically be assigned.
    name : str, optional
        Name of the universe. If not specified, the name is the empty string.
    auto_geom_ids : bool
        Set IDs automatically on initialization (True) or report overlaps in ID
        space between CSG and DAGMC (False)

    Attributes
    ----------
    id : int
        Unique identifier of the universe
    name : str
        Name of the universe
    mesh : openmc.XDGMesh
        Mesh to use for the XDG universe
    auto_geom_ids : bool
        Set IDs automatically on initialization (True) or report overlaps in ID
        space between CSG and XDG geometry (False)
    """

    def __init__(self,
                 mesh: XDGMesh,
                 universe_id=None,
                 name='',
                 auto_geom_ids=False,
                 auto_mat_ids=False):
        super().__init__(universe_id, name)
        # Initialize class attributes
        self.mesh = mesh
        self.auto_geom_ids = auto_geom_ids
        self._type = 'surface_mesh'
        self._background_material = None

    def __repr__(self):
        string = super().__repr__()
        string += '{: <16}=\t{}\n'.format('\tGeom', 'XDG')
        string += '{: <16}=\t{}\n'.format('\tAuto Geom IDs', self.auto_geom_ids)
        string += '{: <16}=\t{}\n'.format('\tType', self._type)
        string += '{: <16}=\t{}\n'.format('\tMesh', self._mesh)
        return string

    @property
    def bounding_box(self):
        return BoundingBox.infinite()

    @property
    def filename(self):
        return self._filename

    @filename.setter
    def filename(self, val: cv.PathLike):
        cv.check_type('XDG filename', val, cv.PathLike)
        self._filename = input_path(val)

    @property
    def material_overrides(self):
        raise NotImplementedError("Material overrides are not implemented for XDG")

    @property
    def mesh(self):
        return self._mesh

    @mesh.setter
    def mesh(self, val):
        cv.check_type('XDG mesh', val, XDGMesh)
        self._mesh = val

    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, val):
        cv.check_value('XDG type', val, ('surface_mesh', 'volume_mesh'))
        self._type = val

    def replace_material_assignment(self, material_name: str, material: openmc.Material):
        """Replace a material assignment within the DAGMC universe.

        Replace the material assignment of all cells filled with a material in
        the DAGMC universe. The universe must be synchronized in an initialized
        Model (see :meth:`~openmc.DAGMCUniverse.sync_dagmc_cells`) before
        calling this method.

        .. versionadded:: 0.15.1

        Parameters
        ----------
        material_name : str
            Material name to replace
        material : openmc.Material
            Material to replace the material_name with

        """
        raise NotImplementedError("Material overrides are not implemented for XDG")

    def add_material_override(self, key, overrides=None):
        """Add a material override to the universe.

        .. versionadded:: 0.15

        Parameters
        ----------
        key : openmc.DAGMCCell or int
            Cell object or ID of the Cell to override
        value : openmc.Material or Iterable of openmc.Material
            Material(s) to be applied to the Cell passed as the key

        """
        raise NotImplementedError("Material overrides are not implemented for XDG")

    @property
    def auto_geom_ids(self):
        return self._auto_geom_ids

    @auto_geom_ids.setter
    def auto_geom_ids(self, val):
        cv.check_type('DAGMC automatic geometry ids', val, bool)
        self._auto_geom_ids = val

    @property
    def material_names(self):
        raise NotImplementedError("Material names is not implemented for XDG")

    @property
    def background_material(self):
        return self._background_material

    @background_material.setter
    def background_material(self, val):
        cv.check_type('XDG background material', val, openmc.Material)
        self._background_material = val

    def _n_geom_elements(self, geom_type):
        """
        Helper function for retrieving the number geometric entities in a DAGMC
        file

        Parameters
        ----------
        geom_type : str
            The type of geometric entity to count. One of {'Volume', 'Surface'}. Returns
            the runtime number of voumes in the DAGMC model (includes implicit complement).

        Returns
        -------
        int
            Number of geometry elements of the specified type
        """
        raise NotImplementedError("Number of cells is not implemented for XDG")

    @property
    def n_cells(self):
        raise NotImplementedError("Number of cells is not implemented for XDG")

    @property
    def n_surfaces(self):
        raise NotImplementedError("Number of surfaces is not implemented for XDG")

    def create_xml_subelement(self, xml_element, memo=None):
        if memo is None:
            memo = set()

        if self in memo:
            return

        memo.add(self)

        # Set xml element values
        xdg_element = ET.Element('xdg_universe')
        xdg_element.set('id', str(self.id))

        if self.auto_geom_ids:
            xdg_element.set('auto_geom_ids', 'true')

        xdg_element.set('mesh', str(self.mesh.id))

        xdg_element.set('type', self._type)

        if self.background_material is not None:
            xdg_element.set('background_material', str(self.background_material.id))

        # add mesh element
        xml_element.append(xdg_element)

        # if this mesh has already been added to the XML element
        if self.mesh in memo:
            return
        memo.add(self.mesh)

        xml_element.append(self.mesh.to_xml_element())

    def bounding_region(
        self,
        bounded_type: str = 'box',
        boundary_type: str = 'vacuum',
        starting_id: int = 10000,
        padding_distance: float = 0.
    ):
        """Creates a either a spherical or box shaped bounding region around
        the DAGMC geometry.

        .. versionadded:: 0.13.1

        Parameters
        ----------
        bounded_type : str
            The type of bounding surface(s) to use when constructing the region.
            Options include a single spherical surface (sphere) or a rectangle
            made from six planes (box).
        boundary_type : str
            Boundary condition that defines the behavior for particles hitting
            the surface. Defaults to vacuum boundary condition. Passed into the
            surface construction.
        starting_id : int
            Starting ID of the surface(s) used in the region. For bounded_type
            'box', the next 5 IDs will also be used. Defaults to 10000 to reduce
            the chance of an overlap of surface IDs with the DAGMC geometry.
        padding_distance : float
            Distance between the bounding region surfaces and the minimal
            bounding box. Allows for the region to be larger than the DAGMC
            geometry.

        Returns
        -------
        openmc.Region
            Region instance
        """
        raise NotImplementedError("Bounded region is not implemented for XDG")

    def bounded_universe(self, bounding_cell_id=10000, **kwargs):
        """Returns an openmc.Universe filled with this DAGMCUniverse and bounded
        with a cell. Defaults to a box cell with a vacuum surface however this
        can be changed using the kwargs which are passed directly to
        DAGMCUniverse.bounding_region().

        Parameters
        ----------
        bounding_cell_id : int
            The cell ID number to use for the bounding cell, defaults to 10000 to reduce
            the chance of overlapping ID numbers with the DAGMC geometry.

        Returns
        -------
        openmc.Universe
            Universe instance
        """
        raise NotImplementedError("Bounded universe is not implemented for XDG")

    @classmethod
    def from_hdf5(cls, group):
        """Create DAGMC universe from HDF5 group

        Parameters
        ----------
        group : h5py.Group
            Group in HDF5 file

        Returns
        -------
        openmc.XDGUniverse
            XDGUniverse instance

        """
        id = int(group.name.split('/')[-1].lstrip('universe '))
        fname = group['filename'][()].decode()
        name = group['name'][()].decode() if 'name' in group else None

        out = cls(fname, universe_id=id, name=name)

        out.auto_geom_ids = bool(group.attrs['auto_geom_ids'])
        out.auto_mat_ids = bool(group.attrs['auto_mat_ids'])

        return out

    @classmethod
    def from_xml_element(cls, elem, meshes=None):
        """Generate XDG universe from XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            `<xdg_universe>` element
        meshes : dict
            Dictionary mapping mesh ID strings to :class:`openmc.XDGMesh`
            instances (defined in :meth:`openmc.Geometry.from_xml`)

        Returns
        -------
        openmc.XDGUniverse
            XDGUniverse instance

        """
        id = int(get_text(elem, 'id'))
        mesh_id = int(get_text(elem, 'mesh'))

        out = cls(meshes[mesh_id], universe_id=id)

        name = get_text(elem, 'name')
        if name is not None:
            out.name = name

        out.auto_geom_ids = bool(elem.get('auto_geom_ids'))

        if type := elem.get('type'):
            out.type = type
        else:
            raise ValueError("XDG type is not specified in the XML element")

        if library := elem.get('library'):
            out.library = library

        if background_material := elem.get('background_material'):
            out.background_material = background_material

        return out

    def _partial_deepcopy(self):
        """Clone all of the openmc.DAGMCUniverse object's attributes except for
        its cells, as they are copied within the clone function. This should
        only to be used within the openmc.UniverseBase.clone() context.
        """
        clone = openmc.XDGUniverse(name=self.name, filename=self.filename)
        clone.auto_geom_ids = self.auto_geom_ids
        clone.auto_mat_ids = self.auto_mat_ids
        clone.library = self.library
        clone.type = self.type
        return clone

    def add_cell(self, cell):
        """Add a cell to the universe.

        Parameters
        ----------
        cell : openmc.XDGCell
            Cell to add

        """
        raise NotImplementedError("Add cell is not implemented for XDG")

    def remove_cell(self, cell):
        """Remove a cell from the universe.

        Parameters
        ----------
        cell : openmc.Cell
            Cell to remove

        """
        raise NotImplementedError("Remove cell is not implemented for XDG")

    def sync_dagmc_cells(self, mats: Iterable[openmc.Material]):
        """Synchronize DAGMC cell information between Python and C API

        .. versionadded:: 0.15.1

        Parameters
        ----------
        mats : iterable of openmc.Material
            Iterable of materials to assign to the DAGMC cells

        """
        raise NotImplementedError("Sync cells is not implemented for XDG")

class XDGCell(openmc.Cell):
    """A cell class for XDG-based geometries.

    Parameters
    ----------
    cell_id : int or None, optional
        Unique identifier for the cell. If None, an identifier will be
        automatically assigned.
    name : str, optional
        Name of the cell.
    fill : openmc.Material or None, optional
        Material filling the cell. If None, the cell is filled with vacuum.

    Attributes
    ----------
    id : int
        Unique identifier of the cell
    name : str
        Name of the cell
    fill : openmc.Material or None, optional
        Material filling the cell. If None, the cell is filled with vacuum.
    """
    pass
