from collections import Iterable, OrderedDict
from xml.etree import ElementTree as ET

import openmc
from openmc.clean_xml import *
from openmc.checkvalue import check_type

def reset_auto_ids():
    openmc.reset_auto_material_id()
    openmc.reset_auto_surface_id()
    openmc.reset_auto_cell_id()
    openmc.reset_auto_universe_id()


class Geometry(object):
    """Geometry representing a collection of surfaces, cells, and universes.

    Parameters
    ----------
    root_universe : openmc.Universe, optional
        Root universe which contains all others

    Attributes
    ----------
    root_universe : openmc.Universe
        Root universe which contains all others

    """

    def __init__(self, root_universe=None):
        self._root_universe = None
        self._offsets = {}
        if root_universe is not None:
            self.root_universe = root_universe

    @property
    def root_universe(self):
        return self._root_universe

    @root_universe.setter
    def root_universe(self, root_universe):
        check_type('root universe', root_universe, openmc.Universe)
        if root_universe._id != 0:
            msg = 'Unable to add root Universe "{0}" to Geometry since ' \
                  'it has ID="{1}" instead of ' \
                  'ID=0'.format(root_universe, root_universe._id)
            raise ValueError(msg)

        self._root_universe = root_universe

    def export_to_xml(self):
        """Create a geometry.xml file that can be used for a simulation.

        """

        # Clear OpenMC written IDs used to optimize XML generation
        openmc.universe.WRITTEN_IDS = {}

        # Create XML representation
        geometry_file = ET.Element("geometry")
        self.root_universe.create_xml_subelement(geometry_file)

        # Clean the indentation in the file to be user-readable
        sort_xml_elements(geometry_file)
        clean_xml_indentation(geometry_file)

        # Write the XML Tree to the geometry.xml file
        tree = ET.ElementTree(geometry_file)
        tree.write("geometry.xml", xml_declaration=True, encoding='utf-8',
                   method="xml")

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

    def get_cell_instance(self, path):
        """Return the instance number for the final cell in a geometry path.

        The instance is an index into tally distribcell filter arrays.

        Parameters
        ----------
        path : list
            A list of IDs that form the path to the target. It should begin with
            0 for the base universe, and should cover every universe, cell, and
            lattice passed through. For the case of the lattice, a tuple should
            be provided to indicate which coordinates in the lattice should be
            entered. This should be in the form: (lat_id, i_x, i_y, i_z)

        Returns
        -------
        instance : int
            Index in tally results array for distribcell filters

        """

        # Extract the cell id from the path
        last_index = path.rfind('>')
        cell_id = int(path[last_index+1:])

        # Find the distribcell index of the cell.
        cells = self.get_all_cells()
        for cell in cells:
            if cell.id == cell_id:
                distribcell_index = cell.distribcell_index
                break
        else:
            raise RuntimeError('Could not find cell {} specified in a \
                                distribcell filter'.format(cell_id))

        # Return memoize'd offset if possible
        if (path, distribcell_index) in self._offsets:
            offset = self._offsets[(path, distribcell_index)]

        # Begin recursive call to compute offset starting with the base Universe
        else:
            offset = self._root_universe.get_cell_instance(path,
                                                           distribcell_index)
            self._offsets[(path, distribcell_index)] = offset

        # Return the final offset
        return offset

    def get_all_cells(self):
        """Return all cells defined

        Returns
        -------
        list of openmc.Cell
            Cells in the geometry

        """

        all_cells = self.root_universe.get_all_cells()
        cells = list(set(all_cells.values()))
        cells.sort(key=lambda x: x.id)
        return cells

    def get_all_universes(self):
        """Return all universes defined

        Returns
        -------
        list of openmc.Universe
            Universes in the geometry

        """

        all_universes = self._root_universe.get_all_universes()
        universes = list(set(all_universes.values()))
        universes.sort(key=lambda x: x.id)
        return universes

    def get_all_nuclides(self):
        """Return all nuclides assigned to a material in the geometry

        Returns
        -------
        list of openmc.Nuclide
            Nuclides in the geometry

        """

        nuclides = OrderedDict()
        materials = self.get_all_materials()

        for material in materials:
            nuclides.update(material.get_all_nuclides())

        return nuclides

    def get_all_materials(self):
        """Return all materials assigned to a cell

        Returns
        -------
        list of openmc.Material
            Materials in the geometry

        """

        material_cells = self.get_all_material_cells()
        materials = []

        for cell in material_cells:
            if cell.fill_type == 'distribmat':
                for m in cell.fill:
                    if m is not None and m not in materials:
                        materials.append(m)
            elif cell.fill_type == 'material':
                if cell.fill not in materials:
                    materials.append(cell.fill)

        materials.sort(key=lambda x: x.id)
        return materials

    def get_all_material_cells(self):
        """Return all cells filled by a material

        Returns
        -------
        list of openmc.Cell
            Cells filled by Materials in the geometry

        """

        all_cells = self.get_all_cells()
        material_cells = []

        for cell in all_cells:
            if cell.fill_type in ('material', 'distribmat'):
                if cell not in material_cells:
                    material_cells.append(cell)

        material_cells.sort(key=lambda x: x.id)
        return material_cells

    def get_all_material_universes(self):
        """Return all universes composed of at least one non-fill cell

        Returns
        -------
        list of openmc.Universe
            Universes with non-fill cells

        """

        all_universes = self.get_all_universes()
        material_universes = []

        for universe in all_universes:
            cells = universe.cells
            for cell in cells:
                if cell.fill_type in ('material', 'distribmat', 'void'):
                    if universe not in material_universes:
                        material_universes.append(universe)

        material_universes.sort(key=lambda x: x.id)
        return material_universes

    def get_all_lattices(self):
        """Return all lattices defined

        Returns
        -------
        list of openmc.Lattice
            Lattices in the geometry

        """

        cells = self.get_all_cells()
        lattices = []

        for cell in cells:
            if cell.fill_type == 'lattice':
                if cell.fill not in lattices:
                    lattices.append(cell.fill)

        lattices.sort(key=lambda x: x.id)
        return lattices

    def get_materials_by_name(self, name, case_sensitive=False, matching=False):
        """Return a list of materials with matching names.

        Parameters
        ----------
        name : str
            The name to match
        case_sensitive : bool
            Whether to distinguish upper and lower case letters in each
            material's name (default is True)
        matching : bool
            Whether the names must match completely (default is True)

        Returns
        -------
        list of openmc.Material
            Materials matching the queried name

        """

        if not case_sensitive:
            name = name.lower()

        all_materials = self.get_all_materials()
        materials = set()

        for material in all_materials:
            material_name = material.name
            if not case_sensitive:
                material_name = material_name.lower()

            if material_name == name:
                materials.add(material)
            elif not matching and name in material_name:
                materials.add(material)

        materials = list(materials)
        materials.sort(key=lambda x: x.id)
        return materials

    def get_cells_by_name(self, name, case_sensitive=False, matching=False):
        """Return a list of cells with matching names.

        Parameters
        ----------
        name : str
            The name to search match
        case_sensitive : bool
            Whether to distinguish upper and lower case letters in each
            cell's name (default is True)
        matching : bool
            Whether the names must match completely (default is True)

        Returns
        -------
        list of openmc.Cell
            Cells matching the queried name

        """

        if not case_sensitive:
            name = name.lower()

        all_cells = self.get_all_cells()
        cells = set()

        for cell in all_cells:
            cell_name = cell.name
            if not case_sensitive:
                cell_name = cell_name.lower()

            if cell_name == name:
                cells.add(cell)
            elif not matching and name in cell_name:
                cells.add(cell)

        cells = list(cells)
        cells.sort(key=lambda x: x.id)
        return cells

    def get_cells_by_fill_name(self, name, case_sensitive=False, matching=False):
        """Return a list of cells with fills with matching names.

        Parameters
        ----------
        name : str
            The name to match
        case_sensitive : bool
            Whether to distinguish upper and lower case letters in each
            cell's name (default is True)
        matching : bool
            Whether the names must match completely (default is True)

        Returns
        -------
        list of openmc.Cell
            Cells with fills matching the queried name

        """

        if not case_sensitive:
            name = name.lower()

        all_cells = self.get_all_cells()
        cells = set()

        for cell in all_cells:
            cell_fill_name = cell.fill.name
            if not case_sensitive:
                cell_fill_name = cell_fill_name.lower()

            if cell_fill_name == name:
                cells.add(cell)
            elif not matching and name in cell_fill_name:
                cells.add(cell)

        cells = list(cells)
        cells.sort(key=lambda x: x.id)
        return cells

    def get_universes_by_name(self, name, case_sensitive=False, matching=False):
        """Return a list of universes with matching names.

        Parameters
        ----------
        name : str
            The name to match
        case_sensitive : bool
            Whether to distinguish upper and lower case letters in each
            universe's name (default is True)
        matching : bool
            Whether the names must match completely (default is True)

        Returns
        -------
        list of openmc.Universe
            Universes matching the queried name

        """

        if not case_sensitive:
            name = name.lower()

        all_universes = self.get_all_universes()
        universes = set()

        for universe in all_universes:
            universe_name = universe.name
            if not case_sensitive:
                universe_name = universe_name.lower()

            if universe_name == name:
                universes.add(universe)
            elif not matching and name in universe_name:
                universes.add(universe)

        universes = list(universes)
        universes.sort(key=lambda x: x.id)
        return universes

    def get_lattices_by_name(self, name, case_sensitive=False, matching=False):
        """Return a list of lattices with matching names.

        Parameters
        ----------
        name : str
            The name to match
        case_sensitive : bool
            Whether to distinguish upper and lower case letters in each
            lattice's name (default is True)
        matching : bool
            Whether the names must match completely (default is True)

        Returns
        -------
        list of openmc.Lattice
            Lattices matching the queried name

        """

        if not case_sensitive:
            name = name.lower()

        all_lattices = self.get_all_lattices()
        lattices = set()

        for lattice in all_lattices:
            lattice_name = lattice.name
            if not case_sensitive:
                lattice_name = lattice_name.lower()

            if lattice_name == name:
                lattices.add(lattice)
            elif not matching and name in lattice_name:
                lattices.add(lattice)

        lattices = list(lattices)
        lattices.sort(key=lambda x: x.id)
        return lattices
