from collections import OrderedDict, Iterable
from numbers import Integral
from xml.etree import ElementTree as ET
import sys
import warnings

import numpy as np

import openmc
import openmc.checkvalue as cv

if sys.version_info[0] >= 3:
    basestring = str


# A dictionary for storing IDs of cell elements that have already been written,
# used to optimize the writing process
WRITTEN_IDS = {}

# A static variable for auto-generated Lattice (Universe) IDs
AUTO_UNIVERSE_ID = 10000


def reset_auto_universe_id():
    global AUTO_UNIVERSE_ID
    AUTO_UNIVERSE_ID = 10000


class Universe(object):
    """A collection of cells that can be repeated.

    Parameters
    ----------
    universe_id : int, optional
        Unique identifier of the universe. If not specified, an identifier will
        automatically be assigned
    name : str, optional
        Name of the universe. If not specified, the name is the empty string.

    Attributes
    ----------
    id : int
        Unique identifier of the universe
    name : str
        Name of the universe
    cells : collections.OrderedDict
        Dictionary whose keys are cell IDs and values are Cell instances

    """

    def __init__(self, universe_id=None, name=''):
        # Initialize Cell class attributes
        self.id = universe_id
        self.name = name

        # Keys     - Cell IDs
        # Values - Cells
        self._cells = OrderedDict()

        # Keys     - Cell IDs
        # Values - Offsets
        self._cell_offsets = OrderedDict()
        self._num_regions = 0

    def __eq__(self, other):
        if not isinstance(other, Universe):
            return False
        elif self.id != other.id:
            return False
        elif self.name != other.name:
            return False
        elif self.cells != other.cells:
            return False
        else:
            return True

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(repr(self))

    def __repr__(self):
        string = 'Universe\n'
        string += '{0: <16}{1}{2}\n'.format('\tID', '=\t', self._id)
        string += '{0: <16}{1}{2}\n'.format('\tName', '=\t', self._name)
        string += '{0: <16}{1}{2}\n'.format('\tCells', '=\t',
                                            list(self._cells.keys()))
        string += '{0: <16}{1}{2}\n'.format('\t# Regions', '=\t',
                                            self._num_regions)
        return string

    @property
    def id(self):
        return self._id

    @property
    def name(self):
        return self._name

    @property
    def cells(self):
        return self._cells

    @id.setter
    def id(self, universe_id):
        if universe_id is None:
            global AUTO_UNIVERSE_ID
            self._id = AUTO_UNIVERSE_ID
            AUTO_UNIVERSE_ID += 1
        else:
            cv.check_type('universe ID', universe_id, Integral)
            cv.check_greater_than('universe ID', universe_id, 0, equality=True)
            self._id = universe_id

    @name.setter
    def name(self, name):
        if name is not None:
            cv.check_type('universe name', name, basestring)
            self._name = name
        else:
            self._name = ''

    def add_cell(self, cell):
        """Add a cell to the universe.

        Parameters
        ----------
        cell : openmc.Cell
            Cell to add

        """

        if not isinstance(cell, openmc.Cell):
            msg = 'Unable to add a Cell to Universe ID="{0}" since "{1}" is not ' \
                  'a Cell'.format(self._id, cell)
            raise ValueError(msg)

        cell_id = cell._id

        if cell_id not in self._cells:
            self._cells[cell_id] = cell

    def add_cells(self, cells):
        """Add multiple cells to the universe.

        Parameters
        ----------
        cells : Iterable of openmc.Cell
            Cells to add

        """

        if not isinstance(cells, Iterable):
            msg = 'Unable to add Cells to Universe ID="{0}" since "{1}" is not ' \
                  'iterable'.format(self._id, cells)
            raise ValueError(msg)

        for cell in cells:
            self.add_cell(cell)

    def remove_cell(self, cell):
        """Remove a cell from the universe.

        Parameters
        ----------
        cell : openmc.Cell
            Cell to remove

        """

        if not isinstance(cell, openmc.Cell):
            msg = 'Unable to remove a Cell from Universe ID="{0}" since "{1}" is ' \
                  'not a Cell'.format(self._id, cell)
            raise ValueError(msg)

        # If the Cell is in the Universe's list of Cells, delete it
        if cell.id in self._cells:
            del self._cells[cell.id]

    def clear_cells(self):
        """Remove all cells from the universe."""

        self._cells.clear()

    def get_cell_instance(self, path, distribcell_index):

        # Pop off the root Universe ID from the path
        next_index = path.index('-')
        path = path[next_index+2:]

        # Extract the Cell ID from the path
        if '-' in path:
            next_index = path.index('-')
            cell_id = int(path[:next_index])
            path = path[next_index+2:]
        else:
            cell_id = int(path)
            path = ''

        # Make a recursive call to the Cell within this Universe
        offset = self.cells[cell_id].get_cell_instance(path, distribcell_index)

        # Return the offset computed at all nested Universe levels
        return offset

    def get_all_nuclides(self):
        """Return all nuclides contained in the universe

        Returns
        -------
        nuclides : collections.OrderedDict
            Dictionary whose keys are nuclide names and values are 2-tuples of
            (nuclide, density)

        """

        nuclides = OrderedDict()

        # Append all Nuclides in each Cell in the Universe to the dictionary
        for cell_id, cell in self._cells.items():
            nuclides.update(cell.get_all_nuclides())

        return nuclides

    def get_all_cells(self):
        """Return all cells that are contained within the universe

        Returns
        -------
        cells : collections.OrderedDict
            Dictionary whose keys are cell IDs and values are Cell instances

        """

        cells = OrderedDict()

        # Add this Universe's cells to the dictionary
        cells.update(self._cells)

        # Append all Cells in each Cell in the Universe to the dictionary
        for cell_id, cell in self._cells.items():
            cells.update(cell.get_all_cells())

        return cells

    def get_all_materials(self):
        """Return all materials that are contained within the universe

        Returns
        -------
        materials : Collections.OrderedDict
            Dictionary whose keys are material IDs and values are Material instances

        """

        materials = OrderedDict()

        # Append all Cells in each Cell in the Universe to the dictionary
        cells = self.get_all_cells()
        for cell_id, cell in cells.items():
            materials.update(cell.get_all_materials())

        return materials

    def get_all_universes(self):
        """Return all universes that are contained within this one.

        Returns
        -------
        universes : collections.OrderedDict
            Dictionary whose keys are universe IDs and values are Universe
            instances

        """

        # Get all Cells in this Universe
        cells = self.get_all_cells()

        universes = OrderedDict()

        # Append all Universes containing each Cell to the dictionary
        for cell_id, cell in cells.items():
            universes.update(cell.get_all_universes())

        return universes

    def create_xml_subelement(self, xml_element):

        # Iterate over all Cells
        for cell_id, cell in self._cells.items():

            # If the cell was not already written, write it
            if cell_id not in WRITTEN_IDS:
                WRITTEN_IDS[cell_id] = None

                # Create XML subelement for this Cell
                cell_subelement = cell.create_xml_subelement(xml_element)

                # Append the Universe ID to the subelement and add to Element
                cell_subelement.set("universe", str(self._id))
                xml_element.append(cell_subelement)
