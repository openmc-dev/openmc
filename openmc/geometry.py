from collections import OrderedDict
from xml.etree import ElementTree as ET

import openmc
from openmc.clean_xml import sort_xml_elements, clean_xml_indentation
from openmc.checkvalue import check_type


def reset_auto_ids():
    """Reset counters for all auto-generated IDs"""
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
        self._root_universe = root_universe

    def add_volume_information(self, volume_calc):
        """Add volume information from a stochastic volume calculation.

        Parameters
        ----------
        volume_calc : openmc.VolumeCalculation
            Results from a stochastic volume calculation

        """
        if volume_calc.domain_type == 'cell':
            for cell in self.get_all_cells():
                if cell.id in volume_calc.volumes:
                    cell.add_volume_information(volume_calc)
        elif volume_calc.domain_type == 'material':
            for material in self.get_all_materials():
                if material.id in volume_calc.volumes:
                    material.add_volume_information(volume_calc)
        elif volume_calc.domain_type == 'universe':
            for universe in self.get_all_universes():
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

        # Clean the indentation in the file to be user-readable
        sort_xml_elements(root_element)
        clean_xml_indentation(root_element)

        # Write the XML Tree to the geometry.xml file
        tree = ET.ElementTree(root_element)
        tree.write(path, xml_declaration=True, encoding='utf-8', method="xml")

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
                if cell.fill not in lattices:
                    lattices[cell.fill.id] = cell.fill

        return lattices

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

        all_cells = self.get_all_cells().values()
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

        lattices = list(lattices)
        lattices.sort(key=lambda x: x.id)
        return lattices

    def count_cell_instances(self):
        """Count the number of instances for each cell in the Geometry, and
        record the count in the :attr:`Cell.num_instances` properties."""

        # (Re-)initialize all cell instances to 0
        for cell in self.get_all_cells().values():
            cell._num_instances = 0

        # Recursively traverse the CSG tree to count all cell instances
        self.root_universe._count_cell_instances()

    def count_material_instances(self):
        """Count the number of instances for each material in the Geometry, and
        record the count in the :attr:`Material.num_instances` properties."""

        # (Re-)initialize all material instances to 0
        for material in self.get_all_materials().values():
            material._num_instances = 0

        # Recursively traverse the CSG tree to count all cell instances
        self.root_universe._count_material_instances()