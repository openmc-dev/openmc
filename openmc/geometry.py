from collections import OrderedDict, defaultdict
from collections.abc import Iterable
from copy import deepcopy
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np

import openmc
from openmc.clean_xml import clean_xml_indentation
from openmc.checkvalue import check_type


class Geometry(object):
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

    def export_to_xml(self, path='geometry.xml'):
        """Export geometry to an XML file.

        Parameters
        ----------
        path : str
            Path to file to write. Defaults to 'geometry.xml'.

        """
        # Create XML representation
        root_element = ET.Element("geometry")
        self.root_universe.create_xml_subelement(root_element)

        # Sort the elements in the file
        root_element[:] = sorted(root_element, key=lambda x: (
            x.tag, int(x.get('id'))))

        # Clean the indentation in the file to be user-readable
        clean_xml_indentation(root_element)

        # Write the XML Tree to the geometry.xml file
        tree = ET.ElementTree(root_element)
        tree.write(path, xml_declaration=True, encoding='utf-8')

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
        # Helper function to get value from an attribute/element
        def get(elem, name, default=None):
            if name in elem.attrib:
                return elem.get(name, default)
            else:
                child = elem.find(name)
                return child.text if child is not None else default

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
            other_id = get(surface, 'periodic_surface_id')
            if other_id is not None:
                periodic[s.id] = int(other_id)

        # Apply periodic surfaces
        for s1, s2 in periodic.items():
            surfaces[s1].periodic_surface = surfaces[s2]

        # Dictionary that maps each universe to a list of cells/lattices that
        # contain it (needed to determine which universe is the root)
        child_of = defaultdict(list)

        for elem in root.findall('lattice'):
            lat_id = int(get(elem, 'id'))
            name = get(elem, 'name')
            lat = openmc.RectLattice(lat_id, name)
            lat.lower_left = [float(i) for i in get(elem, 'lower_left').split()]
            lat.pitch = [float(i) for i in get(elem, 'pitch').split()]
            outer = get(elem, 'outer')
            if outer is not None:
                lat.outer = get_universe(int(outer))
            universes[lat_id] = lat

            # Get array of universes
            dimension = get(elem, 'dimension').split()
            shape = np.array(dimension, dtype=int)[::-1]
            uarray = np.array([get_universe(int(i)) for i in
                               get(elem, 'universes').split()])
            for u in uarray:
                child_of[u].append(lat)
            uarray.shape = shape
            lat.universes = uarray

        for elem in root.findall('hex_lattice'):
            lat_id = int(get(elem, 'id'))
            name = get(elem, 'name')
            lat = openmc.HexLattice(lat_id, name)
            lat.center = [float(i) for i in get(elem, 'center').split()]
            lat.pitch = [float(i) for i in get(elem, 'pitch').split()]
            outer = get(elem, 'outer')
            if outer is not None:
                lat.outer = get_universe(int(outer))
            universes[lat_id] = lat

            # Get nested lists of universes
            lat._num_rings = n_rings = int(get(elem, 'n_rings'))
            lat._num_axial = n_axial = int(get(elem, 'n_axial', 1))

            # Create empty nested lists for one axial level
            univs = [[None for _ in range(max(6*(n_rings - 1 - r), 1))]
                     for r in range(n_rings)]
            if n_axial > 1:
                univs = [deepcopy(univs) for i in range(n_axial)]

            # Get flat array of universes numbers
            uarray = np.array([get_universe(int(i)) for i in
                               get(elem, 'universes').split()])
            for u in uarray:
                child_of[u].append(lat)

            # Fill nested lists
            j = 0
            for z in range(n_axial):
                # Get list for a single axial level
                axial_level = univs[z] if n_axial > 1 else univs

                # Start iterating from top
                x, alpha = 0, n_rings - 1
                while True:
                    # Set entry in list based on (x,alpha,z) coordinates
                    _, i_ring, i_within = lat.get_universe_index((x, alpha, z))
                    axial_level[i_ring][i_within] = uarray[j]

                    # Move to the right
                    x += 2
                    alpha -= 1
                    if not lat.is_valid_index((x, alpha, z)):
                        # Move down in y direction
                        alpha += x - 1
                        x = 1 - x
                        if not lat.is_valid_index((x, alpha, z)):
                            # Move to the right
                            x += 2
                            alpha -= 1
                            if not lat.is_valid_index((x, alpha, z)):
                                # Reached the bottom
                                break
                    j += 1
            lat.universes = univs

        # Create dictionary to easily look up materials
        if materials is None:
            filename = Path(path).parent / 'materials.xml'
            materials = openmc.Materials.from_xml(str(filename))
        mats = {str(m.id): m for m in materials}
        mats['void'] = None

        for elem in root.findall('cell'):
            cell_id = int(get(elem, 'id'))
            name = get(elem, 'name')
            c = openmc.Cell(cell_id, name)

            # Assign material/distributed materials or fill
            mat_text = get(elem, 'material')
            if mat_text is not None:
                mat_ids = mat_text.split()
                if len(mat_ids) > 1:
                    c.fill = [mats[i] for i in mat_ids]
                else:
                    c.fill = mats[mat_ids[0]]
            else:
                fill_id = int(get(elem, 'fill'))
                c.fill = get_universe(fill_id)
                child_of[c.fill].append(c)

            # Assign region
            region = get(elem, 'region')
            if region is not None:
                c.region = openmc.Region.from_expression(region, surfaces)

            # Check for other attributes
            t = get(elem, 'temperature')
            if t is not None:
                c.temperature = float(t)
            for key in ('temperature', 'rotation', 'translation'):
                value = get(elem, key)
                if value is not None:
                    setattr(c, key, [float(x) for x in value.split()])

            # Add this cell to appropriate universe
            univ_id = int(get(elem, 'universe', 0))
            get_universe(univ_id).add_cell(c)

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
        return self.root_universe.get_all_cells()

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
        return self.root_universe.get_all_materials()

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
            surfaces = cell.region.get_surfaces(surfaces)
        return surfaces

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

        if not case_sensitive:
            name = name.lower()

        all_materials = self.get_all_materials().values()
        materials = set()

        for material in all_materials:
            material_name = material.name
            if not case_sensitive:
                material_name = material_name.lower()

            if material_name == name:
                materials.add(material)
            elif not matching and name in material_name:
                materials.add(material)

        return sorted(materials, key=lambda x: x.id)

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

        if not case_sensitive:
            name = name.lower()

        all_cells = self.get_all_cells().values()
        cells = set()

        for cell in all_cells:
            cell_name = cell.name
            if not case_sensitive:
                cell_name = cell_name.lower()

            if cell_name == name:
                cells.add(cell)
            elif not matching and name in cell_name:
                cells.add(cell)

        return sorted(cells, key=lambda x: x.id)

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

        if not case_sensitive:
            name = name.lower()

        all_universes = self.get_all_universes().values()
        universes = set()

        for universe in all_universes:
            universe_name = universe.name
            if not case_sensitive:
                universe_name = universe_name.lower()

            if universe_name == name:
                universes.add(universe)
            elif not matching and name in universe_name:
                universes.add(universe)

        return sorted(universes, key=lambda x: x.id)

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

        if not case_sensitive:
            name = name.lower()

        all_lattices = self.get_all_lattices().values()
        lattices = set()

        for lattice in all_lattices:
            lattice_name = lattice.name
            if not case_sensitive:
                lattice_name = lattice_name.lower()

            if lattice_name == name:
                lattices.add(lattice)
            elif not matching and name in lattice_name:
                lattices.add(lattice)

        return sorted(lattices, key=lambda x: x.id)

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
