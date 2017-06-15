from __future__ import division
from collections import OrderedDict, Iterable
from copy import copy, deepcopy
from numbers import Integral, Real
import random
import sys

from six import string_types
import numpy as np

import openmc
import openmc.checkvalue as cv
from openmc.plots import _SVG_COLORS
from openmc.mixin import IDManagerMixin


class Universe(IDManagerMixin):
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
    volume : float
        Volume of the universe in cm^3. This can either be set manually or
        calculated in a stochastic volume calculation and added via the
        :meth:`Universe.add_volume_information` method.
    bounding_box : 2-tuple of numpy.array
        Lower-left and upper-right coordinates of an axis-aligned bounding box
        of the universe.

    """

    next_id = 1
    used_ids = set()

    def __init__(self, universe_id=None, name='', cells=None):
        # Initialize Cell class attributes
        self.id = universe_id
        self.name = name
        self._volume = None
        self._atoms = {}

        # Keys     - Cell IDs
        # Values - Cells
        self._cells = OrderedDict()

        if cells is not None:
            self.add_cells(cells)

    def __eq__(self, other):
        if not isinstance(other, Universe):
            return False
        elif self.id != other.id:
            return False
        elif self.name != other.name:
            return False
        elif dict.__ne__(self.cells, other.cells):
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
    def name(self):
        return self._name

    @property
    def cells(self):
        return self._cells

    @property
    def volume(self):
        return self._volume

    @property
    def bounding_box(self):
        regions = [c.region for c in self.cells.values()
                   if c.region is not None]
        if regions:
            return openmc.Union(regions).bounding_box
        else:
            # Infinite bounding box
            return openmc.Intersection().bounding_box

    @name.setter
    def name(self, name):
        if name is not None:
            cv.check_type('universe name', name, string_types)
            self._name = name
        else:
            self._name = ''

    @volume.setter
    def volume(self, volume):
        if volume is not None:
            cv.check_type('universe volume', volume, Real)
        self._volume = volume

    @classmethod
    def from_hdf5(cls, group, cells):
        """Create universe from HDF5 group

        Parameters
        ----------
        group : h5py.Group
            Group in HDF5 file
        cells : dict
            Dictionary mapping cell IDs to instances of :class:`openmc.Cell`.

        Returns
        -------
        openmc.Universe
            Universe instance

        """
        universe_id = int(group.name.split('/')[-1].lstrip('universe '))
        cell_ids = group['cells'].value

        # Create this Universe
        universe = cls(universe_id)

        # Add each Cell to the Universe
        for cell_id in cell_ids:
            universe.add_cell(cells[cell_id])

        return universe

    def add_volume_information(self, volume_calc):
        """Add volume information to a universe.

        Parameters
        ----------
        volume_calc : openmc.VolumeCalculation
            Results from a stochastic volume calculation

        """
        if volume_calc.domain_type == 'universe':
            if self.id in volume_calc.volumes:
                self._volume = volume_calc.volumes[self.id][0]
                self._atoms = volume_calc.atoms[self.id]
            else:
                raise ValueError('No volume information found for this universe.')
        else:
            raise ValueError('No volume information found for this universe.')

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

    def plot(self, origin=(0., 0., 0.), width=(1., 1.), pixels=(200, 200),
             basis='xy', color_by='cell', colors=None, filename=None, seed=None,
             **kwargs):
        """Display a slice plot of the universe.

        Parameters
        ----------
        origin : Iterable of float
            Coordinates at the origin of the plot
        width : Iterable of float
            Width of the plot in each basis direction
        pixels : Iterable of int
            Number of pixels to use in each basis direction
        basis : {'xy', 'xz', 'yz'}
            The basis directions for the plot
        color_by : {'cell', 'material'}
            Indicate whether the plot should be colored by cell or by material
        colors : dict
            Assigns colors to specific materials or cells. Keys are instances of
            :class:`Cell` or :class:`Material` and values are RGB 3-tuples, RGBA
            4-tuples, or strings indicating SVG color names. Red, green, blue,
            and alpha should all be floats in the range [0.0, 1.0], for example:

            .. code-block:: python

               # Make water blue
               water = openmc.Cell(fill=h2o)
               universe.plot(..., colors={water: (0., 0., 1.))

        filename : str or None
            Filename to save plot to. If no filename is given, the plot will be
            displayed using the currently enabled matplotlib backend.
        seed : hashable object or None
            Hashable object which is used to seed the random number generator
            used to select colors. If None, the generator is seeded from the
            current time.
        **kwargs
            All keyword arguments are passed to
            :func:`matplotlib.pyplot.imshow`.

        """
        import matplotlib.pyplot as plt

        # Seed the random number generator
        if seed is not None:
            random.seed(seed)

        if colors is None:
            # Create default dictionary if none supplied
            colors = {}
        else:
            # Convert to RGBA if necessary
            colors = copy(colors)
            for obj, color in colors.items():
                if isinstance(color, string_types):
                    if color.lower() not in _SVG_COLORS:
                        raise ValueError("'{}' is not a valid color."
                                         .format(color))
                    colors[obj] = [x/255 for x in
                                   _SVG_COLORS[color.lower()]] + [1.0]
                elif len(color) == 3:
                    colors[obj] = list(color) + [1.0]

        if basis == 'xy':
            x_min = origin[0] - 0.5*width[0]
            x_max = origin[0] + 0.5*width[0]
            y_min = origin[1] - 0.5*width[1]
            y_max = origin[1] + 0.5*width[1]
        elif basis == 'yz':
            # The x-axis will correspond to physical y and the y-axis will
            # correspond to physical z
            x_min = origin[1] - 0.5*width[0]
            x_max = origin[1] + 0.5*width[0]
            y_min = origin[2] - 0.5*width[1]
            y_max = origin[2] + 0.5*width[1]
        elif basis == 'xz':
            # The y-axis will correspond to physical z
            x_min = origin[0] - 0.5*width[0]
            x_max = origin[0] + 0.5*width[0]
            y_min = origin[2] - 0.5*width[1]
            y_max = origin[2] + 0.5*width[1]

        # Determine locations to determine cells at
        x_coords = np.linspace(x_min, x_max, pixels[0], endpoint=False) + \
                   0.5*(x_max - x_min)/pixels[0]
        y_coords = np.linspace(y_max, y_min, pixels[1], endpoint=False) - \
                   0.5*(y_max - y_min)/pixels[1]

        # Initialize output image in RGBA format.  Flip the pixels from
        # traditional (x, y) to (y, x) used in graphics.
        img = np.zeros((pixels[1], pixels[0], 4))
        for i, x in enumerate(x_coords):
            for j, y in enumerate(y_coords):
                if basis == 'xy':
                    path = self.find((x, y, origin[2]))
                elif basis == 'yz':
                    path = self.find((origin[0], x, y))
                elif basis == 'xz':
                    path = self.find((x, origin[1], y))

                if len(path) > 0:
                    try:
                        if color_by == 'cell':
                            obj = path[-1]
                        elif color_by == 'material':
                            if path[-1].fill_type == 'material':
                                obj = path[-1].fill
                            else:
                                continue
                    except AttributeError:
                        continue
                    if obj not in colors:
                        colors[obj] = (random.random(), random.random(),
                                       random.random(), 1.0)
                    img[j, i, :] = colors[obj]

        # Display image
        plt.imshow(img, extent=(x_min, x_max, y_min, y_max),
                   interpolation='nearest', **kwargs)

        # Show or save the plot
        if filename is None:
            plt.show()
        else:
            plt.savefig(filename)

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

        cell_id = cell.id

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

    def get_nuclides(self):
        """Returns all nuclides in the universe

        Returns
        -------
        nuclides : list of str
            List of nuclide names

        """

        nuclides = []

        # Append all Nuclides in each Cell in the Universe to the dictionary
        for cell in self.cells.values():
            for nuclide in cell.get_nuclides():
                if nuclide not in nuclides:
                    nuclides.append(nuclide)

        return nuclides

    def get_nuclide_densities(self):
        """Return all nuclides contained in the universe

        Returns
        -------
        nuclides : collections.OrderedDict
            Dictionary whose keys are nuclide names and values are 2-tuples of
            (nuclide, density)

        """
        nuclides = OrderedDict()

        if self._atoms is not None:
            volume = self.volume
            for name, atoms in self._atoms.items():
                nuclide = openmc.Nuclide(name)
                density = 1.0e-24 * atoms[0]/volume  # density in atoms/b-cm
                nuclides[name] = (nuclide, density)
        else:
            raise RuntimeError(
                'Volume information is needed to calculate microscopic cross '
                'sections for universe {}. This can be done by running a '
                'stochastic volume calculation via the '
                'openmc.VolumeCalculation object'.format(self.id))

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
        for cell in self._cells.values():
            cells.update(cell.get_all_cells())

        return cells

    def get_all_materials(self):
        """Return all materials that are contained within the universe

        Returns
        -------
        materials : collections.OrderedDict
            Dictionary whose keys are material IDs and values are
            :class:`Material` instances

        """

        materials = OrderedDict()

        # Append all Cells in each Cell in the Universe to the dictionary
        cells = self.get_all_cells()
        for cell in cells.values():
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
        # Append all Universes within each Cell to the dictionary
        universes = OrderedDict()
        for cell in self.get_all_cells().values():
            universes.update(cell.get_all_universes())

        return universes

    def clone(self, memo=None):
        """Create a copy of this universe with a new unique ID, and clones
        all cells within this universe.

        Parameters
        ----------
        memo : dict or None
            A nested dictionary of previously cloned objects. This parameter
            is used internally and should not be specified by the user.

        Returns
        -------
        clone : openmc.Universe
            The clone of this universe

        """

        if memo is None:
            memo = {}

        # If no nemoize'd clone exists, instantiate one
        if self not in memo:
            clone = deepcopy(self)
            clone.id = None

            # Clone all cells for the universe clone
            clone._cells = OrderedDict()
            for cell in self._cells.values():
                clone.add_cell(cell.clone(memo))

            # Memoize the clone
            memo[self] = clone

        return memo[self]

    def create_xml_subelement(self, xml_element):
        # Iterate over all Cells
        for cell_id, cell in self._cells.items():
            path = "./cell[@id='{}']".format(cell_id)

            # If the cell was not already written, write it
            if xml_element.find(path) is None:
                # Create XML subelement for this Cell
                cell_element = cell.create_xml_subelement(xml_element)

                # Append the Universe ID to the subelement and add to Element
                cell_element.set("universe", str(self._id))
                xml_element.append(cell_element)

    def _determine_paths(self, path='', instances_only=False):
        """Count the number of instances for each cell in the universe, and
        record the count in the :attr:`Cell.num_instances` properties."""

        univ_path = path + 'u{}'.format(self.id)

        for cell in self.cells.values():
            cell_path = '{}->c{}'.format(univ_path, cell.id)

            # If universe-filled, recursively count cells in filling universe
            if cell.fill_type == 'universe':
                cell.fill._determine_paths(cell_path + '->', instances_only)

            # If lattice-filled, recursively call for all universes in lattice
            elif cell.fill_type == 'lattice':
                latt = cell.fill

                # Count instances in each universe in the lattice
                for index in latt._natural_indices:
                    latt_path = '{}->l{}({})->'.format(
                        cell_path, latt.id, ",".join(str(x) for x in index))
                    univ = latt.get_universe(index)
                    univ._determine_paths(latt_path, instances_only)

            else:
                if cell.fill_type == 'material':
                    mat = cell.fill
                elif cell.fill_type == 'distribmat':
                    mat = cell.fill[cell._num_instances]
                else:
                    mat = None

                if mat is not None:
                    mat._num_instances += 1
                    if not instances_only:
                        mat._paths.append('{}->m{}'.format(cell_path, mat.id))

            # Append current path
            cell._num_instances += 1
            if not instances_only:
                cell._paths.append(cell_path)
