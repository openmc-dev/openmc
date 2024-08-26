from collections.abc import Iterable
from numbers import Integral
from pathlib import Path

import h5py
import lxml.etree as ET
import numpy as np

import openmc
import openmc.checkvalue as cv
from ._xml import get_text
from .checkvalue import check_type, check_value
from .surface import _BOUNDARY_TYPES


class DAGMCUniverse(openmc.UniverseBase):
    """A reference to a DAGMC file to be used in the model.

    .. versionadded:: 0.13.0

    Parameters
    ----------
    filename : str
        Path to the DAGMC file used to represent this universe.
    universe_id : int, optional
        Unique identifier of the universe. If not specified, an identifier will
        automatically be assigned.
    name : str, optional
        Name of the universe. If not specified, the name is the empty string.
    auto_geom_ids : bool
        Set IDs automatically on initialization (True) or report overlaps in ID
        space between CSG and DAGMC (False)
    auto_mat_ids : bool
        Set IDs automatically on initialization (True)  or report overlaps in ID
        space between OpenMC and UWUW materials (False)

    Attributes
    ----------
    id : int
        Unique identifier of the universe
    name : str
        Name of the universe
    filename : str
        Path to the DAGMC file used to represent this universe.
    auto_geom_ids : bool
        Set IDs automatically on initialization (True) or report overlaps in ID
        space between CSG and DAGMC (False)
    auto_mat_ids : bool
        Set IDs automatically on initialization (True)  or report overlaps in ID
        space between OpenMC and UWUW materials (False)
    bounding_box : openmc.BoundingBox
        Lower-left and upper-right coordinates of an axis-aligned bounding box
        of the universe.

        .. versionadded:: 0.13.1
    material_names : list of str
        Return a sorted list of materials names that are contained within the
        DAGMC h5m file. This is useful when naming openmc.Material() objects
        as each material name present in the DAGMC h5m file must have a
        matching openmc.Material() with the same name.

        .. versionadded:: 0.13.2
    n_cells : int
        The number of cells in the DAGMC model. This is the number of cells at
        runtime and accounts for the implicit complement whether or not is it
        present in the DAGMC file.

        .. versionadded:: 0.13.2
    n_surfaces : int
        The number of surfaces in the model.

        .. versionadded:: 0.15
    mat_overrides : dict
        A dictionary of material overrides. The keys are material names as
        strings and the values are openmc.Material objects. If a material name
        is found in the DAGMC file, the material will be replaced with the
        openmc.Material object in the value.

    
    """

    def __init__(self,
                 filename,
                 universe_id=None,
                 name='',
                 auto_geom_ids=False,
                 auto_mat_ids=False,
                 mat_overrides={}):
        super().__init__(universe_id, name)
        # Initialize class attributes
        self.filename = filename
        self.auto_geom_ids = auto_geom_ids
        self.auto_mat_ids = auto_mat_ids
        self.material_overrides = mat_overrides
        self._dagmc_cells = []

    def __repr__(self):
        string = super().__repr__()
        string += '{: <16}=\t{}\n'.format('\tGeom', 'DAGMC')
        string += '{: <16}=\t{}\n'.format('\tFile', self.filename)
        return string

    @property
    def bounding_box(self):
        with h5py.File(self.filename) as dagmc_file:
            coords = dagmc_file['tstt']['nodes']['coordinates'][()]
            lower_left_corner = coords.min(axis=0)
            upper_right_corner = coords.max(axis=0)
            return openmc.BoundingBox(lower_left_corner, upper_right_corner)

    @property
    def filename(self):
        return self._filename

    @filename.setter
    def filename(self, val):
        cv.check_type('DAGMC filename', val, (Path, str))
        self._filename = val

    @property
    def auto_geom_ids(self):
        return self._auto_geom_ids

    @auto_geom_ids.setter
    def auto_geom_ids(self, val):
        cv.check_type('DAGMC automatic geometry ids', val, bool)
        self._auto_geom_ids = val

    @property
    def auto_mat_ids(self):
        return self._auto_mat_ids

    @auto_mat_ids.setter
    def auto_mat_ids(self, val):
        cv.check_type('DAGMC automatic material ids', val, bool)
        self._auto_mat_ids = val

    @property
    def material_names(self):
        dagmc_file_contents = h5py.File(self.filename)
        material_tags_hex = dagmc_file_contents['/tstt/tags/NAME'].get(
            'values')
        material_tags_ascii = []
        for tag in material_tags_hex:
            candidate_tag = tag.tobytes().decode().replace('\x00', '')
            # tags might be for temperature or reflective surfaces
            if candidate_tag.startswith('mat:'):
                # if name ends with _comp remove it, it is not parsed
                if candidate_tag.endswith('_comp'):
                   candidate_tag = candidate_tag[:-5]
                # removes first 4 characters as openmc.Material name should be
                # set without the 'mat:' part of the tag
                material_tags_ascii.append(candidate_tag[4:])

        return sorted(set(material_tags_ascii))

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
        cv.check_value('geometry type', geom_type, ('volume', 'surface'))

        def decode_str_tag(tag_val):
            return tag_val.tobytes().decode().replace('\x00', '')

        dagmc_filepath = Path(self.filename).resolve()
        with h5py.File(dagmc_filepath) as dagmc_file:
            category_data = dagmc_file['tstt/tags/CATEGORY/values']
            category_strs = map(decode_str_tag, category_data)
            n = sum([v == geom_type.capitalize() for v in category_strs])

            # check for presence of an implicit complement in the file and
            # increment the number of cells if it doesn't exist
            if geom_type == 'volume':
                name_data = dagmc_file['tstt/tags/NAME/values']
                name_strs = map(decode_str_tag, name_data)
                if not sum(['impl_complement' in n for n in name_strs]):
                    n += 1
        return n

    @property
    def n_cells(self):
        return self._n_geom_elements('volume')

    @property
    def n_surfaces(self):
        return self._n_geom_elements('surface')

    def create_xml_subelement(self, xml_element, memo=None):
        if memo is None:
            memo = set()

        if self in memo:
            return

        memo.add(self)

        # Set xml element values
        dagmc_element = ET.Element('dagmc_universe')
        dagmc_element.set('id', str(self.id))

        if self.auto_geom_ids:
            dagmc_element.set('auto_geom_ids', 'true')
        if self.auto_mat_ids:
            dagmc_element.set('auto_mat_ids', 'true')
        dagmc_element.set('filename', str(self.filename))
        if len(self.material_overrides) == 0:
            mats = self.get_all_materials()
            for mat in mats.values():
                if mat.name[:-4] in self.material_names:
                    self.material_overrides.setdefault(
                        mat.name[:-4].lower(), []).append(mat.name)
                    print(f"Material {mat.name} found in DAGMC file, ")

        if self.material_overrides:
            mat_element = ET.Element('material_overrides')
            for key in self.material_overrides:
                mat_element.set(key, ' '.join(
                    t for t in self.material_overrides[key]))
            dagmc_element.append(mat_element)
        xml_element.append(dagmc_element)

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

        check_type('boundary type', boundary_type, str)
        check_value('boundary type', boundary_type, _BOUNDARY_TYPES)
        check_type('starting surface id', starting_id, Integral)
        check_type('bounded type', bounded_type, str)
        check_value('bounded type', bounded_type, ('box', 'sphere'))

        bbox = self.bounding_box.expand(padding_distance, True)

        if bounded_type == 'sphere':
            radius = np.linalg.norm(bbox.upper_right - bbox.center)
            bounding_surface = openmc.Sphere(
                surface_id=starting_id,
                x0=bbox.center[0],
                y0=bbox.center[1],
                z0=bbox.center[2],
                boundary_type=boundary_type,
                r=radius,
            )

            return -bounding_surface

        if bounded_type == 'box':
            # defines plane surfaces for all six faces of the bounding box
            lower_x = openmc.XPlane(bbox[0][0], surface_id=starting_id)
            upper_x = openmc.XPlane(bbox[1][0], surface_id=starting_id+1)
            lower_y = openmc.YPlane(bbox[0][1], surface_id=starting_id+2)
            upper_y = openmc.YPlane(bbox[1][1], surface_id=starting_id+3)
            lower_z = openmc.ZPlane(bbox[0][2], surface_id=starting_id+4)
            upper_z = openmc.ZPlane(bbox[1][2], surface_id=starting_id+5)

            region = +lower_x & -upper_x & +lower_y & -upper_y & +lower_z & -upper_z

            for surface in region.get_surfaces().values():
                surface.boundary_type = boundary_type

            return region

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
        bounding_cell = openmc.Cell(
            fill=self, cell_id=bounding_cell_id, region=self.bounding_region(**kwargs))
        return openmc.Universe(cells=[bounding_cell])

    @classmethod
    def from_hdf5(cls, group):
        """Create DAGMC universe from HDF5 group

        Parameters
        ----------
        group : h5py.Group
            Group in HDF5 file

        Returns
        -------
        openmc.DAGMCUniverse
            DAGMCUniverse instance

        """
        id = int(group.name.split('/')[-1].lstrip('universe '))
        fname = group['filename'][()].decode()
        name = group['name'][()].decode() if 'name' in group else None

        out = cls(fname, universe_id=id, name=name)

        out.auto_geom_ids = bool(group.attrs['auto_geom_ids'])
        out.auto_mat_ids = bool(group.attrs['auto_mat_ids'])

        return out

    @classmethod
    def from_xml_element(cls, elem):
        """Generate DAGMC universe from XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            `<dagmc_universe>` element

        Returns
        -------
        openmc.DAGMCUniverse
            DAGMCUniverse instance

        """
        id = int(get_text(elem, 'id'))
        fname = get_text(elem, 'filename')

        out = cls(fname, universe_id=id)

        name = get_text(elem, 'name')
        if name is not None:
            out.name = name

        out.auto_geom_ids = bool(elem.get('auto_geom_ids'))
        out.auto_mat_ids = bool(elem.get('auto_mat_ids'))

        return out

    def _partial_deepcopy(self):
        """Clone all of the openmc.DAGMCUniverse object's attributes except for
        its cells, as they are copied within the clone function. This should
        only to be used within the openmc.UniverseBase.clone() context.
        """
        clone = openmc.DAGMCUniverse(name=self.name, filename=self.filename)
        clone.volume = self.volume
        clone.auto_geom_ids = self.auto_geom_ids
        clone.auto_mat_ids = self.auto_mat_ids
        return clone

    def add_cell(self, cell):
        """Add a cell to the universe.

        Parameters
        ----------
        cell : openmc.DAGMCCell
            Cell to add

        """
        if not isinstance(cell, openmc.DAGMCCell):
            msg = f'Unable to add a DAGMCCell to DAGMCUniverse ID="{self._id}" since ' \
                f'"{cell}" is not a DAGMCCell'
            raise TypeError(msg)

        cell_id = cell.id

        if cell_id not in self._cells:
            self._cells[cell_id] = cell

    def remove_cell(self, cell):
        """Remove a cell from the universe.

        Parameters
        ----------
        cell : openmc.Cell
            Cell to remove

        """

        if not isinstance(cell, openmc.DAGMCCell):
            msg = f'Unable to remove a Cell from Universe ID="{self._id}" ' \
                f'since "{cell}" is not a Cell'
            raise TypeError(msg)

        # If the Cell is in the Universe's list of Cells, delete it
        self._cells.pop(cell.id, None)

    def sync_dagmc_cells(self, mats={}):
        """Synchronize DAGMC cell information between Python and C API

        .. versionadded:: 0.13.0

        """
        import openmc.lib
        if not openmc.lib.is_initialized:
            raise RuntimeError("Model must be initialized via Model.init_lib "
                               "before calling this method.")

        mats_per_id = {mat.id: mat for mat in mats}
        for dag_cell_id in openmc.lib.dagmc.get_dagmc_cell_ids(self.id, self._n_geom_elements('volume')):
            dag_cell = openmc.lib.cells[dag_cell_id]
            if isinstance(dag_cell.fill, Iterable):
                fill = [mats_per_id[mat_id.id]
                        for mat_id in dag_cell.fill]
            else:
                fill = mats_per_id[dag_cell.fill.id]
            dag_pseudo_cell = openmc.DAGMCCell(
                cell_id=dag_cell_id, fill=fill
            )
            self.add_cell(dag_pseudo_cell)


class DAGMCCell(openmc.Cell):
    def __init__(self, cell_id=None, name='', fill=None, region=None):
        super().__init__(cell_id, name, fill, region)

    @property
    def DAG_parent_universe(self):
        """Get the parent universe of the cell."""
        return self._parent_universe

    @DAG_parent_universe.setter
    def DAG_parent_universe(self, universe):
        """Set the parent universe of the cell."""
        self._parent_universe = universe.id

    def boundingbox(self):
        print("Warning: Bounding box is not available for cells in a DAGMC universe.")
        return {}

    def get_all_cells(self, memo=None):
        return {}

    def get_all_universes(self, memo=None):
        return {}

    def clone(self, clone_materials=True, clone_regions=True, memo=None):
        print("Warning: clone is not available for cells in a DAGMC universe.")
        return None

    def plot(self, *args, **kwargs):
        print("Warning: plot is not available for cells in a DAGMC universe.")
        return None

    def create_xml_subelement(self, xml_element, memo=None):
        print(
            "Warning: create_xml_subelement is not available for cells in a DAGMC universe.")
        return None

    @classmethod
    def from_xml_element(cls, elem, surfaces, materials, get_universe):
        print("Warning: from_xml_element is not available for cells in a DAGMC universe.")
        return None
