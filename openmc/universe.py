from collections import OrderedDict, Iterable
from numbers import Integral
from xml.etree import ElementTree as ET
import random
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
    cells : Iterable of openmc.Cell, optional
        Cells to add to the universe. By default no cells are added.

    Attributes
    ----------
    id : int
        Unique identifier of the universe
    name : str
        Name of the universe
    cells : collections.OrderedDict
        Dictionary whose keys are cell IDs and values are :class:`Cell`
        instances

    """

    def __init__(self, universe_id=None, name='', cells=None):
        # Initialize Cell class attributes
        self.id = universe_id
        self.name = name

        # Keys     - Cell IDs
        # Values - Cells
        self._cells = OrderedDict()

        # Keys     - Cell IDs
        # Values - Offsets
        self._cell_offsets = OrderedDict()

        if cells is not None:
            self.add_cells(cells)

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
        p = np.asarray(point)
        for cell in self._cells.values():
            if p in cell:
                if cell.fill_type in ('material', 'distribmat', 'void'):
                    return [self, cell]
                elif cell.fill_type == 'universe':
                    if cell.translation is not None:
                        p -= cell.translation
                    if cell.rotation is not None:
                        p[:] = cell.rotation_matrix.dot(p)
                    return [self, cell] + cell.fill.find(p)
                else:
                    return [self, cell] + cell.fill.find(p)
        return []

    def plot(self, center=(0., 0., 0.), width=(1., 1.), pixels=(200, 200),
             basis='xy', color_by='cell', seed=None):
        """Display a slice plot of the universe.

        Parameters
        ----------
        center : Iterable of float
            Coordinates at the center of the plot
        width : Iterable of float
            Width of the plot in each basis direction
        pixels : Iterable of int
            Number of pixels to use in each basis direction
        basis : {'xy', 'xz', 'yz'}
            The basis directions for the plot
        color_by : {'cell', 'material'}
            Indicate whether the plot should be colored by cell or by material
        seed : hashable object or None
            Hashable object which is used to seed the random number generator
            used to select colors. If None, the generator is seeded from the
            current time.

        """
        import matplotlib.pyplot as plt

        # Seed the random number generator
        if seed is not None:
            random.seed(seed)

        if basis == 'xy':
            x_min = center[0] - 0.5*width[0]
            x_max = center[0] + 0.5*width[0]
            y_min = center[1] - 0.5*width[1]
            y_max = center[1] + 0.5*width[1]
        elif basis == 'yz':
            # The x-axis will correspond to physical y and the y-axis will correspond to physical z
            x_min = center[1] - 0.5*width[0]
            x_max = center[1] + 0.5*width[0]
            y_min = center[2] - 0.5*width[1]
            y_max = center[2] + 0.5*width[1]
        elif basis == 'xz':
            # The y-axis will correspond to physical z
            x_min = center[0] - 0.5*width[0]
            x_max = center[0] + 0.5*width[0]
            y_min = center[2] - 0.5*width[1]
            y_max = center[2] + 0.5*width[1]

        # Determine locations to determine cells at
        x_coords = np.linspace(x_min, x_max, pixels[0], endpoint=False) + \
                   0.5*(x_max - x_min)/pixels[0]
        y_coords = np.linspace(y_max, y_min, pixels[1], endpoint=False) - \
                   0.5*(y_max - y_min)/pixels[1]

        colors = {}
        img = np.zeros(pixels + (4,))  # Use RGBA form
        for i, x in enumerate(x_coords):
            for j, y in enumerate(y_coords):
                if basis == 'xy':
                    path = self.find((x, y, center[2]))
                elif basis == 'yz':
                    path = self.find((center[0], x, y))
                elif basis == 'xz':
                    path = self.find((x, center[1], y))

                if len(path) > 0:
                    try:
                        if color_by == 'cell':
                            uid = path[-1].id
                        elif color_by == 'material':
                            if path[-1].fill_type == 'material':
                                uid = path[-1].fill.id
                            else:
                                continue
                    except AttributeError:
                        continue
                    if uid not in colors:
                        colors[uid] = (random.random(), random.random(),
                                       random.random(), 1.0)
                    img[j,i,:] = colors[uid]

        plt.imshow(img, extent=(x_min, x_max, y_min, y_max))
        plt.show()

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
            Dictionary whose keys are cell IDs and values are :class:`Cell`
            instances

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
            Dictionary whose keys are material IDs and values are
            :class:`Material` instances

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
            Dictionary whose keys are universe IDs and values are
            :class:`Universe` instances

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
