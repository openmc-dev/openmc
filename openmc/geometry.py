from __future__ import annotations
import os
from collections import defaultdict
from copy import deepcopy
from collections.abc import Iterable
from pathlib import Path
import warnings
import lxml.etree as ET

import numpy as np

import openmc
import openmc._xml as xml
from .checkvalue import check_type, check_less_than, check_greater_than, PathLike


class Geometry:
    """Geometry representing a collection of surfaces, cells, and universes.

    Parameters
    ----------
    root : openmc.UniverseBase or Iterable of openmc.Cell, optional
        Root universe which contains all others, or an iterable of cells that
        should be used to create a root universe.

    Attributes
    ----------
    root_universe : openmc.UniverseBase
        Root universe which contains all others
    bounding_box : openmc.BoundingBox
        Lower-left and upper-right coordinates of an axis-aligned bounding box
        of the universe.
    merge_surfaces : bool
        Whether to remove redundant surfaces when the geometry is exported.
    surface_precision : int
        Number of decimal places to round to for comparing the coefficients of
        surfaces for considering them topologically equivalent.

    """

    def __init__(
        self,
        root: openmc.UniverseBase | Iterable[openmc.Cell] | None = None,
        merge_surfaces: bool = False,
        surface_precision: int = 10
    ):
        self._root_universe = None
        self._offsets = {}
        self.merge_surfaces = merge_surfaces
        self.surface_precision = surface_precision
        if root is not None:
            if isinstance(root, openmc.UniverseBase):
                self.root_universe = root
            else:
                univ = openmc.Universe()
                for cell in root:
                    univ.add_cell(cell)
                self._root_universe = univ

    @property
    def root_universe(self) -> openmc.UniverseBase:
        return self._root_universe

    @root_universe.setter
    def root_universe(self, root_universe):
        check_type('root universe', root_universe, openmc.UniverseBase)
        self._root_universe = root_universe

    @property
    def bounding_box(self) -> np.ndarray:
        return self.root_universe.bounding_box

    @property
    def merge_surfaces(self) -> bool:
        return self._merge_surfaces

    @merge_surfaces.setter
    def merge_surfaces(self, merge_surfaces):
        check_type('merge surfaces', merge_surfaces, bool)
        self._merge_surfaces = merge_surfaces

    @property
    def surface_precision(self) -> int:
        return self._surface_precision

    @surface_precision.setter
    def surface_precision(self, surface_precision):
        check_type('surface precision', surface_precision, int)
        check_less_than('surface_precision', surface_precision, 16)
        check_greater_than('surface_precision', surface_precision, 0)
        self._surface_precision = surface_precision

    def add_volume_information(self, volume_calc):
        """Add volume information from a stochastic volume calculation.

        Parameters
        ----------
        volume_calc : openmc.VolumeCalculation
            Results from a stochastic volume calculation

        """
        if volume_calc.domain_type == 'cell':
            for cell in self.get_all_cells().values():
                if cell.id in volume_calc.volumes:
                    cell.add_volume_information(volume_calc)
        elif volume_calc.domain_type == 'material':
            for material in self.get_all_materials().values():
                if material.id in volume_calc.volumes:
                    material.add_volume_information(volume_calc)
        elif volume_calc.domain_type == 'universe':
            for universe in self.get_all_universes().values():
                if universe.id in volume_calc.volumes:
                    universe.add_volume_information(volume_calc)

    def to_xml_element(self, remove_surfs=False) -> ET.Element:
        """Creates a 'geometry' element to be written to an XML file.

        Parameters
        ----------
        remove_surfs : bool
            Whether or not to remove redundant surfaces from the geometry when
            exporting

        """
        # Find and remove redundant surfaces from the geometry
        if remove_surfs:
            warnings.warn("remove_surfs kwarg will be deprecated soon, please "
                          "set the Geometry.merge_surfaces attribute instead.")
            self.merge_surfaces = True

        if self.merge_surfaces:
            self.remove_redundant_surfaces()

        # Create XML representation
        element = ET.Element("geometry")
        self.root_universe.create_xml_subelement(element)

        # Sort the elements in the file
        element[:] = sorted(element, key=lambda x: (
            x.tag, int(x.get('id'))))

        # Clean the indentation in the file to be user-readable
        xml.clean_indentation(element)
        xml.reorder_attributes(element)  # TODO: Remove when support is Python 3.8+

        return element

    def export_to_xml(self, path='geometry.xml', remove_surfs=False):
        """Export geometry to an XML file.

        Parameters
        ----------
        path : str
            Path to file to write. Defaults to 'geometry.xml'.
        remove_surfs : bool
            Whether or not to remove redundant surfaces from the geometry when
            exporting

            .. versionadded:: 0.12

        """
        root_element = self.to_xml_element(remove_surfs)

        # Check if path is a directory
        p = Path(path)
        if p.is_dir():
            p /= 'geometry.xml'

        # Write the XML Tree to the geometry.xml file
        tree = ET.ElementTree(root_element)
        tree.write(str(p), xml_declaration=True, encoding='utf-8')

    @classmethod
    def from_xml_element(cls, elem, materials=None) -> Geometry:
        """Generate geometry from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element
        materials : openmc.Materials or None
            Materials used to assign to cells. If None, an attempt is made to
            generate it from the materials.xml file.

        Returns
        -------
        openmc.Geometry
            Geometry object

        """
        mats = dict()
        if materials is not None:
            mats.update({str(m.id): m for m in materials})
        mats['void'] = None

        # Helper function for keeping a cache of Universe instances
        universes = {}
        def get_universe(univ_id):
            if univ_id not in universes:
                univ = openmc.Universe(univ_id)
                universes[univ_id] = univ
            return universes[univ_id]

        # Get surfaces
        surfaces = {}
        periodic = {}
        for surface in elem.findall('surface'):
            s = openmc.Surface.from_xml_element(surface)
            surfaces[s.id] = s

            # Check for periodic surface
            other_id = xml.get_text(surface, 'periodic_surface_id')
            if other_id is not None:
                periodic[s.id] = int(other_id)

        # Apply periodic surfaces
        for s1, s2 in periodic.items():
            surfaces[s1].periodic_surface = surfaces[s2]

        # Add any DAGMC universes
        for e in elem.findall('dagmc_universe'):
            dag_univ = openmc.DAGMCUniverse.from_xml_element(e)
            universes[dag_univ.id] = dag_univ

        # Dictionary that maps each universe to a list of cells/lattices that
        # contain it (needed to determine which universe is the elem)
        child_of = defaultdict(list)

        for e in elem.findall('lattice'):
            lat = openmc.RectLattice.from_xml_element(e, get_universe)
            universes[lat.id] = lat
            if lat.outer is not None:
                child_of[lat.outer].append(lat)
            for u in lat.universes.ravel():
                child_of[u].append(lat)

        for e in elem.findall('hex_lattice'):
            lat = openmc.HexLattice.from_xml_element(e, get_universe)
            universes[lat.id] = lat
            if lat.outer is not None:
                child_of[lat.outer].append(lat)
            if lat.ndim == 2:
                for ring in lat.universes:
                    for u in ring:
                        child_of[u].append(lat)
            else:
                for axial_slice in lat.universes:
                    for ring in axial_slice:
                        for u in ring:
                            child_of[u].append(lat)

        for e in elem.findall('cell'):
            c = openmc.Cell.from_xml_element(e, surfaces, mats, get_universe)
            if c.fill_type in ('universe', 'lattice'):
                child_of[c.fill].append(c)

        # Determine which universe is the root by finding one which is not a
        # child of any other object
        for u in universes.values():
            if not child_of[u]:
                return cls(u)
        else:
            raise ValueError('Error determining root universe.')

    @classmethod
    def from_xml(
        cls,
        path: PathLike = 'geometry.xml',
        materials: PathLike | 'openmc.Materials' | None = 'materials.xml'
    ) -> Geometry:
        """Generate geometry from XML file

        Parameters
        ----------
        path : PathLike, optional
            Path to geometry XML file
        materials : openmc.Materials or PathLike
            Materials used to assign to cells. If PathLike, an attempt is made
            to generate materials from the provided xml file.

        Returns
        -------
        openmc.Geometry
            Geometry object

        """

        # Using str and os.PathLike here to avoid error when using just the imported PathLike
        # TypeError: Subscripted generics cannot be used with class and instance checks
        check_type('materials', materials, (str, os.PathLike, openmc.Materials))

        if isinstance(materials, (str, os.PathLike)):
            materials = openmc.Materials.from_xml(materials)

        parser = ET.XMLParser(huge_tree=True)
        tree = ET.parse(path, parser=parser)
        root = tree.getroot()

        return cls.from_xml_element(root, materials)

    def find(self, point) -> list:
        """Find cells/universes/lattices which contain a given point

        Parameters
        ----------
        point : 3-tuple of float
            Cartesian coordinates of the point

        Returns
        -------
        list
            Sequence of universes, cells, and lattices which are traversed to
            find the given point

        """
        return self.root_universe.find(point)

    def get_instances(self, paths) -> int | list[int]:
        """Return the instance number(s) for a cell/material in a geometry path.

        The instance numbers are used as indices into distributed
        material/temperature arrays and tally distribcell filter arrays.

        Parameters
        ----------
        paths : str or iterable of str
            The path traversed through the CSG tree to reach a cell or material
            instance. For example, 'u0->c10->l20(2,2,1)->u5->c5' would indicate
            the cell instance whose first level is universe 0 and cell 10,
            second level is lattice 20 position (2,2,1), and third level is
            universe 5 and cell 5.

        Returns
        -------
        int or list of int
            Instance number(s) for the given path(s)

        """
        # Make sure we are working with an iterable
        return_list = (isinstance(paths, Iterable) and
                       not isinstance(paths, str))
        path_list = paths if return_list else [paths]

        indices = []
        for p in path_list:
            # Extract the cell id from the path
            last_index = p.rfind('>')
            last_path = p[last_index+1:]
            uid = int(last_path[1:])

            # Get corresponding cell/material
            if last_path[0] == 'c':
                obj = self.get_all_cells()[uid]
            elif last_path[0] == 'm':
                obj = self.get_all_materials()[uid]

            # Determine index in paths array
            try:
                indices.append(obj.paths.index(p))
            except ValueError:
                indices.append(None)

        return indices if return_list else indices[0]

    def get_all_cells(self) -> dict[int, openmc.Cell]:
        """Return all cells in the geometry.

        Returns
        -------
        dict
            Dictionary mapping cell IDs to :class:`openmc.Cell` instances

        """
        if self.root_universe is not None:
            return self.root_universe.get_all_cells()
        else:
            return {}

    def get_all_universes(self) -> dict[int, openmc.Universe]:
        """Return all universes in the geometry.

        Returns
        -------
        dict
            Dictionary mapping universe IDs to :class:`openmc.Universe`
            instances

        """
        universes = {}
        universes[self.root_universe.id] = self.root_universe
        universes.update(self.root_universe.get_all_universes())
        return universes

    def get_all_nuclides(self) -> list[str]:
        """Return all nuclides within the geometry.

        Returns
        -------
        list
            Sorted list of all nuclides in materials appearing in the geometry

        """
        all_nuclides = set()
        for material in self.get_all_materials().values():
            all_nuclides |= set(material.get_nuclides())
        return sorted(all_nuclides)
    
    def get_all_dagmc_universes(self) -> typing.Dict[int, openmc.DAGMCUniverse]:
        """Return all universes in the geometry.

        Returns
        -------
        dict
            Dictionary mapping universe IDs to :class:`openmc.DAGMCUniverse`
            instances

        """
        universes = self.get_all_universes()
        dag_universes = {}
        for id, uni in universes.items():
            if isinstance(uni, openmc.DAGMCUniverse):
                dag_universes[id] = uni
        return dag_universes

    def get_all_materials(self) -> dict[int, openmc.Material]:
        """Return all materials within the geometry.

        Returns
        -------
        dict
            Dictionary mapping material IDs to :class:`openmc.Material`
            instances

        """
        if self.root_universe is not None:
            return self.root_universe.get_all_materials()
        else:
            return {}

    def get_all_material_cells(self) -> dict[int, openmc.Cell]:
        """Return all cells filled by a material

        Returns
        -------
        dict
            Dictionary mapping cell IDs to :class:`openmc.Cell` instances that
            are filled with materials or distributed materials.

        """
        material_cells = {}

        for cell in self.get_all_cells().values():
            if cell.fill_type in ('material', 'distribmat'):
                if cell not in material_cells:
                    material_cells[cell.id] = cell

        return material_cells

    def get_all_material_universes(self) -> dict[int, openmc.Universe]:
        """Return all universes having at least one material-filled cell.

        This method can be used to find universes that have at least one cell
        that is filled with a material or is void.

        Returns
        -------
        dict
            Dictionary mapping universe IDs to :class:`openmc.Universe`
            instances with at least one material-filled cell

        """
        material_universes = {}

        for universe in self.get_all_universes().values():
            for cell in universe.cells.values():
                if cell.fill_type in ('material', 'distribmat', 'void'):
                    if universe not in material_universes:
                        material_universes[universe.id] = universe

        return material_universes

    def get_all_lattices(self) -> dict[int, openmc.Lattice]:
        """Return all lattices defined

        Returns
        -------
        dict
            Dictionary mapping lattice IDs to :class:`openmc.Lattice` instances

        """
        lattices = {}

        for cell in self.get_all_cells().values():
            if cell.fill_type == 'lattice':
                if cell.fill.id not in lattices:
                    lattices[cell.fill.id] = cell.fill

        return lattices

    def get_all_surfaces(self) -> dict[int, openmc.Surface]:
        """
        Return all surfaces used in the geometry

        Returns
        -------
        dict
            Dictionary mapping surface IDs to :class:`openmc.Surface` instances

        """
        surfaces = {}

        for cell in self.get_all_cells().values():
            if cell.region is not None:
                surfaces = cell.region.get_surfaces(surfaces)
        return surfaces

    def _get_domains_by_name(self, name, case_sensitive, matching, domain_type) -> list:
        if not case_sensitive:
            name = name.lower()

        domains = []

        func = getattr(self, f'get_all_{domain_type}s')
        for domain in func().values():
            domain_name = domain.name if case_sensitive else domain.name.lower()
            if domain_name == name:
                domains.append(domain)
            elif not matching and name in domain_name:
                domains.append(domain)

        domains.sort(key=lambda x: x.id)
        return domains

    def get_materials_by_name(
        self, name, case_sensitive=False, matching=False
    ) -> list[openmc.Material]:
        """Return a list of materials with matching names.

        Parameters
        ----------
        name : str
            The name to match
        case_sensitive : bool
            Whether to distinguish upper and lower case letters in each
            material's name (default is False)
        matching : bool
            Whether the names must match completely (default is False)

        Returns
        -------
        list of openmc.Material
            Materials matching the queried name

        """
        return self._get_domains_by_name(name, case_sensitive, matching, 'material')

    def get_cells_by_name(
        self, name, case_sensitive=False, matching=False
    ) -> list[openmc.Cell]:
        """Return a list of cells with matching names.

        Parameters
        ----------
        name : str
            The name to search match
        case_sensitive : bool
            Whether to distinguish upper and lower case letters in each
            cell's name (default is False)
        matching : bool
            Whether the names must match completely (default is False)

        Returns
        -------
        list of openmc.Cell
            Cells matching the queried name

        """
        return self._get_domains_by_name(name, case_sensitive, matching, 'cell')

    def get_surfaces_by_name(
        self, name, case_sensitive=False, matching=False
    ) -> list[openmc.Surface]:
        """Return a list of surfaces with matching names.

        .. versionadded:: 0.13.3

        Parameters
        ----------
        name : str
            The name to search match
        case_sensitive : bool
            Whether to distinguish upper and lower case letters in each
            surface's name (default is False)
        matching : bool
            Whether the names must match completely (default is False)

        Returns
        -------
        list of openmc.Surface
            Surfaces matching the queried name

        """
        return self._get_domains_by_name(name, case_sensitive, matching, 'surface')

    def get_cells_by_fill_name(
        self, name, case_sensitive=False, matching=False
    ) -> list[openmc.Cell]:
        """Return a list of cells with fills with matching names.

        Parameters
        ----------
        name : str
            The name to match
        case_sensitive : bool
            Whether to distinguish upper and lower case letters in each
            cell's name (default is False)
        matching : bool
            Whether the names must match completely (default is False)

        Returns
        -------
        list of openmc.Cell
            Cells with fills matching the queried name

        """

        if not case_sensitive:
            name = name.lower()

        cells = set()

        for cell in self.get_all_cells().values():
            names = []
            if cell.fill_type in ('material', 'universe', 'lattice'):
                names.append(cell.fill.name)
            elif cell.fill_type == 'distribmat':
                for mat in cell.fill:
                    if mat is not None:
                        names.append(mat.name)

            for fill_name in names:
                if not case_sensitive:
                    fill_name = fill_name.lower()

                if fill_name == name:
                    cells.add(cell)
                elif not matching and name in fill_name:
                    cells.add(cell)

        return sorted(cells, key=lambda x: x.id)

    def get_universes_by_name(
        self, name, case_sensitive=False, matching=False
    ) -> list[openmc.Universe]:
        """Return a list of universes with matching names.

        Parameters
        ----------
        name : str
            The name to match
        case_sensitive : bool
            Whether to distinguish upper and lower case letters in each
            universe's name (default is False)
        matching : bool
            Whether the names must match completely (default is False)

        Returns
        -------
        list of openmc.Universe
            Universes matching the queried name

        """
        return self._get_domains_by_name(name, case_sensitive, matching, 'universe')

    def get_lattices_by_name(
        self, name, case_sensitive=False, matching=False
    ) -> list[openmc.Lattice]:
        """Return a list of lattices with matching names.

        Parameters
        ----------
        name : str
            The name to match
        case_sensitive : bool
            Whether to distinguish upper and lower case letters in each
            lattice's name (default is False)
        matching : bool
            Whether the names must match completely (default is False)

        Returns
        -------
        list of openmc.Lattice
            Lattices matching the queried name

        """
        return self._get_domains_by_name(name, case_sensitive, matching, 'lattice')

    def remove_redundant_surfaces(self) -> dict[int, openmc.Surface]:
        """Remove and return all of the redundant surfaces.

        Uses surface_precision attribute of Geometry instance for rounding and
        comparing surface coefficients.

        .. versionadded:: 0.12

        Returns
        -------
        redundant_surfaces
            Dictionary whose keys are the ID of a redundant surface and whose
            values are the topologically equivalent :class:`openmc.Surface`
            that should replace it.

        """
        # Get redundant surfaces
        redundancies = defaultdict(list)
        for surf in self.get_all_surfaces().values():
            coeffs = tuple(round(surf._coefficients[k],
                                 self.surface_precision)
                           for k in surf._coeff_keys)
            key = (surf._type, surf._boundary_type) + coeffs
            redundancies[key].append(surf)

        redundant_surfaces = {replace.id: keep
                              for keep, *redundant in redundancies.values()
                              for replace in redundant}

        if redundant_surfaces:
            # Iterate through all cells contained in the geometry
            for cell in self.get_all_cells().values():
                # Recursively remove redundant surfaces from regions
                if cell.region:
                    cell.region.remove_redundant_surfaces(redundant_surfaces)

        return redundant_surfaces

    def determine_paths(self, instances_only=False):
        """Determine paths through CSG tree for cells and materials.

        This method recursively traverses the CSG tree to determine each unique
        path that reaches every cell and material. The paths are stored in the
        :attr:`Cell.paths` and :attr:`Material.paths` attributes.

        Parameters
        ----------
        instances_only : bool, optional
            If true, this method will only determine the number of instances of
            each cell and material.

        """
        # (Re-)initialize all cell instances to 0
        for cell in self.get_all_cells().values():
            cell._paths = []
            cell._num_instances = 0
        for material in self.get_all_materials().values():
            material._paths = []
            material._num_instances = 0

        # Recursively traverse the CSG tree to count all cell instances
        self.root_universe._determine_paths(instances_only=instances_only)

    def clone(self) -> Geometry:
        """Create a copy of this geometry with new unique IDs for all of its
        enclosed materials, surfaces, cells, universes and lattices."""

        clone = deepcopy(self)
        clone.root_universe = self.root_universe.clone()
        return clone

    def plot(self, *args, **kwargs):
        """Display a slice plot of the geometry.

        .. versionadded:: 0.14.0

        Parameters
        ----------
        origin : iterable of float
            Coordinates at the origin of the plot. If left as None then the
            bounding box center will be used to attempt to ascertain the origin.
            Defaults to (0, 0, 0) if the bounding box is not finite
        width : iterable of float
            Width of the plot in each basis direction. If left as none then the
            bounding box width will be used to attempt to ascertain the plot
            width. Defaults to (10, 10) if the bounding box is not finite
        pixels : Iterable of int or int
            If iterable of ints provided, then this directly sets the number of
            pixels to use in each basis direction. If int provided, then this
            sets the total number of pixels in the plot and the number of pixels
            in each basis direction is calculated from this total and the image
            aspect ratio.
        basis : {'xy', 'xz', 'yz'}
            The basis directions for the plot
        color_by : {'cell', 'material'}
            Indicate whether the plot should be colored by cell or by material
        colors : dict
            Assigns colors to specific materials or cells. Keys are instances of
            :class:`Cell` or :class:`Material` and values are RGB 3-tuples, RGBA
            4-tuples, or strings indicating SVG color names. Red, green, blue,
            and alpha should all be floats in the range [0.0, 1.0], for
            example::

               # Make water blue
               water = openmc.Cell(fill=h2o)
               universe.plot(..., colors={water: (0., 0., 1.))
        seed : int
            Seed for the random number generator
        openmc_exec : str
            Path to OpenMC executable.
        axes : matplotlib.Axes
            Axes to draw to
        legend : bool
            Whether a legend showing material or cell names should be drawn
        legend_kwargs : dict
            Keyword arguments passed to :func:`matplotlib.pyplot.legend`.
        outline : bool
            Whether outlines between color boundaries should be drawn
        axis_units : {'km', 'm', 'cm', 'mm'}
            Units used on the plot axis
        **kwargs
            Keyword arguments passed to :func:`matplotlib.pyplot.imshow`
        Returns
        -------
        matplotlib.axes.Axes
            Axes containing resulting image
        """

        return self.root_universe.plot(*args, **kwargs)
