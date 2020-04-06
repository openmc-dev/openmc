from collections import OrderedDict, defaultdict
from collections.abc import Iterable
from copy import deepcopy
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np

import openmc
import openmc._xml as xml
from openmc.checkvalue import check_type


class Geometry:
    """Geometry representing a collection of surfaces, cells, and universes.

    Parameters
    ----------
    root : openmc.Universe or Iterable of openmc.Cell, optional
        Root universe which contains all others, or an iterable of cells that
        should be used to create a root universe.

    Attributes
    ----------
    root_universe : openmc.Universe
        Root universe which contains all others
    bounding_box : 2-tuple of numpy.array
        Lower-left and upper-right coordinates of an axis-aligned bounding box
        of the universe.

    """

    def __init__(self, root=None):
        self._root_universe = None
        self._offsets = {}
        if root is not None:
            if isinstance(root, openmc.Universe):
                self.root_universe = root
            else:
                univ = openmc.Universe()
                for cell in root:
                    univ.add_cell(cell)
                self._root_universe = univ

    @property
    def root_universe(self):
        return self._root_universe

    @property
    def bounding_box(self):
        return self.root_universe.bounding_box

    @root_universe.setter
    def root_universe(self, root_universe):
        check_type('root universe', root_universe, openmc.Universe)
        self._root_universe = root_universe

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

    def export_to_xml(self, path='geometry.xml', remove_surfs=False):
        """Export geometry to an XML file.

        Parameters
        ----------
        path : str
            Path to file to write. Defaults to 'geometry.xml'.
        remove_surfs : bool
            Whether or not to remove redundant surfaces from the geometry when
            exporting

        """
        # Find and remove redundant surfaces from the geometry
        if remove_surfs:
            self.remove_redundant_surfaces()

        # Create XML representation
        root_element = ET.Element("geometry")
        self.root_universe.create_xml_subelement(root_element, memo=set())

        # Sort the elements in the file
        root_element[:] = sorted(root_element, key=lambda x: (
            x.tag, int(x.get('id'))))

        # Clean the indentation in the file to be user-readable
        xml.clean_indentation(root_element)

        # Check if path is a directory
        p = Path(path)
        if p.is_dir():
            p /= 'geometry.xml'

        # Write the XML Tree to the geometry.xml file
        tree = ET.ElementTree(root_element)
        tree.write(str(p), xml_declaration=True, encoding='utf-8')

    @classmethod
    def from_xml(cls, path='geometry.xml', materials=None):
        """Generate geometry from XML file

        Parameters
        ----------
        path : str, optional
            Path to geometry XML file
        materials : openmc.Materials or None
            Materials used to assign to cells. If None, an attempt is made to
            generate it from the materials.xml file.

        Returns
        -------
        openmc.Geometry
            Geometry object

        """
        # Helper function for keeping a cache of Universe instances
        universes = {}
        def get_universe(univ_id):
            if univ_id not in universes:
                univ = openmc.Universe(univ_id)
                universes[univ_id] = univ
            return universes[univ_id]

        tree = ET.parse(path)
        root = tree.getroot()

        # Get surfaces
        surfaces = {}
        periodic = {}
        for surface in root.findall('surface'):
            s = openmc.Surface.from_xml_element(surface)
            surfaces[s.id] = s

            # Check for periodic surface
            other_id = xml.get_text(surface, 'periodic_surface_id')
            if other_id is not None:
                periodic[s.id] = int(other_id)

        # Apply periodic surfaces
        for s1, s2 in periodic.items():
            surfaces[s1].periodic_surface = surfaces[s2]

        # Dictionary that maps each universe to a list of cells/lattices that
        # contain it (needed to determine which universe is the root)
        child_of = defaultdict(list)

        for elem in root.findall('lattice'):
            lat = openmc.RectLattice.from_xml_element(elem, get_universe)
            universes[lat.id] = lat
            if lat.outer is not None:
                child_of[lat.outer].append(lat)
            for u in lat.universes.ravel():
                child_of[u].append(lat)

        for elem in root.findall('hex_lattice'):
            lat = openmc.HexLattice.from_xml_element(elem, get_universe)
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

        # Create dictionary to easily look up materials
        if materials is None:
            filename = Path(path).parent / 'materials.xml'
            materials = openmc.Materials.from_xml(str(filename))
        mats = {str(m.id): m for m in materials}
        mats['void'] = None

        for elem in root.findall('cell'):
            c = openmc.Cell.from_xml_element(elem, surfaces, mats, get_universe)
            if c.fill_type in ('universe', 'lattice'):
                child_of[c.fill].append(c)

        # Determine which universe is the root by finding one which is not a
        # child of any other object
        for u in universes.values():
            if not child_of[u]:
                return cls(u)
        else:
            raise ValueError('Error determining root universe.')

    def find(self, point):
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

    def get_instances(self, paths):
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

    def get_all_cells(self):
        """Return all cells in the geometry.

        Returns
        -------
        collections.OrderedDict
            Dictionary mapping cell IDs to :class:`openmc.Cell` instances

        """
        if self.root_universe is not None:
            return self.root_universe.get_all_cells(memo=set())
        else:
            return OrderedDict()

    def get_all_universes(self):
        """Return all universes in the geometry.

        Returns
        -------
        collections.OrderedDict
            Dictionary mapping universe IDs to :class:`openmc.Universe`
            instances

        """
        universes = OrderedDict()
        universes[self.root_universe.id] = self.root_universe
        universes.update(self.root_universe.get_all_universes())
        return universes

    def get_all_materials(self):
        """Return all materials within the geometry.

        Returns
        -------
        collections.OrderedDict
            Dictionary mapping material IDs to :class:`openmc.Material`
            instances

        """
        if self.root_universe is not None:
            return self.root_universe.get_all_materials(memo=set())
        else:
            return OrderedDict()

    def get_all_material_cells(self):
        """Return all cells filled by a material

        Returns
        -------
        collections.OrderedDict
            Dictionary mapping cell IDs to :class:`openmc.Cell` instances that
            are filled with materials or distributed materials.

        """
        material_cells = OrderedDict()

        for cell in self.get_all_cells().values():
            if cell.fill_type in ('material', 'distribmat'):
                if cell not in material_cells:
                    material_cells[cell.id] = cell

        return material_cells

    def get_all_material_universes(self):
        """Return all universes having at least one material-filled cell.

        This method can be used to find universes that have at least one cell
        that is filled with a material or is void.

        Returns
        -------
        collections.OrderedDict
            Dictionary mapping universe IDs to :class:`openmc.Universe`
            instances with at least one material-filled cell

        """
        material_universes = OrderedDict()

        for universe in self.get_all_universes().values():
            for cell in universe.cells.values():
                if cell.fill_type in ('material', 'distribmat', 'void'):
                    if universe not in material_universes:
                        material_universes[universe.id] = universe

        return material_universes

    def get_all_lattices(self):
        """Return all lattices defined

        Returns
        -------
        collections.OrderedDict
            Dictionary mapping lattice IDs to :class:`openmc.Lattice` instances

        """
        lattices = OrderedDict()

        for cell in self.get_all_cells().values():
            if cell.fill_type == 'lattice':
                if cell.fill.id not in lattices:
                    lattices[cell.fill.id] = cell.fill

        return lattices

    def get_all_surfaces(self):
        """
        Return all surfaces used in the geometry

        Returns
        -------
        collections.OrderedDict
            Dictionary mapping surface IDs to :class:`openmc.Surface` instances

        """
        surfaces = OrderedDict()

        for cell in self.get_all_cells().values():
            if cell.region is not None:
                surfaces = cell.region.get_surfaces(surfaces)
        return surfaces

    def get_redundant_surfaces(self):
        """Return all of the topologically redundant surface ids

        Returns
        -------
        dict
            Dictionary whose keys are the ID of a redundant surface and whose
            values are the topologically equivalent :class:`openmc.Surface`
            that should replace it.

        """
        tally = defaultdict(list)
        for surf in self.get_all_surfaces().values():
            coeffs = tuple(surf._coefficients[k] for k in surf._coeff_keys)
            key = (surf._type,) + coeffs
            tally[key].append(surf)
        return {replace.id: keep
                for keep, *redundant in tally.values()
                for replace in redundant}

    def _get_domains_by_name(self, name, case_sensitive, matching, domain_type):
        if not case_sensitive:
            name = name.lower()

        domains = []

        func = getattr(self, 'get_all_{}s'.format(domain_type))
        for domain in func().values():
            domain_name = domain.name if case_sensitive else domain.name.lower()
            if domain_name == name:
                domains.append(domain)
            elif not matching and name in domain_name:
                domains.append(domain)

        domains.sort(key=lambda x: x.id)
        return domains

    def get_materials_by_name(self, name, case_sensitive=False, matching=False):
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

    def get_cells_by_name(self, name, case_sensitive=False, matching=False):
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

    def get_cells_by_fill_name(self, name, case_sensitive=False, matching=False):
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

    def get_universes_by_name(self, name, case_sensitive=False, matching=False):
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

    def get_lattices_by_name(self, name, case_sensitive=False, matching=False):
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

    def remove_redundant_surfaces(self):
        """Remove redundant surfaces from the geometry"""

        # Get redundant surfaces
        redundant_surfaces = self.get_redundant_surfaces()

        # Iterate through all cells contained in the geometry
        for cell in self.get_all_cells().values():
            # Recursively remove redundant surfaces from regions
            cell.region.remove_redundant_surfaces(redundant_surfaces)

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

    def clone(self):
        """Create a copy of this geometry with new unique IDs for all of its
        enclosed materials, surfaces, cells, universes and lattices."""

        clone = deepcopy(self)
        clone.root_universe = self.root_universe.clone()
        return clone
