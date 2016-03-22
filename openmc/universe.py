import abc
from collections import OrderedDict, Iterable
from numbers import Real, Integral
from xml.etree import ElementTree as ET
import sys
import warnings

import numpy as np

import openmc
import openmc.checkvalue as cv
from openmc.surface import Halfspace
from openmc.region import Region, Intersection, Complement

if sys.version_info[0] >= 3:
    basestring = str


# A static variable for auto-generated Cell IDs
AUTO_CELL_ID = 10000

# A dictionary for storing IDs of cell elements that have already been written,
# used to optimize the writing process
WRITTEN_IDS = {}


def reset_auto_cell_id():
    global AUTO_CELL_ID
    AUTO_CELL_ID = 10000


class Cell(object):
    """A region of space defined as the intersection of half-space created by
    quadric surfaces.

    Parameters
    ----------
    cell_id : int, optional
        Unique identifier for the cell. If not specified, an identifier will
        automatically be assigned.
    name : str, optional
        Name of the cell. If not specified, the name is the empty string.

    Attributes
    ----------
    id : int
        Unique identifier for the cell
    name : str
        Name of the cell
    fill : Material or Universe or Lattice or 'void' or iterable of Material
        Indicates what the region of space is filled with
    region : openmc.region.Region
        Region of space that is assigned to the cell.
    rotation : ndarray
        If the cell is filled with a universe, this array specifies the angles
        in degrees about the x, y, and z axes that the filled universe should be
        rotated.
    translation : ndarray
        If the cell is filled with a universe, this array specifies a vector
        that is used to translate (shift) the universe.
    offsets : ndarray
        Array of offsets used for distributed cell searches
    distribcell_index : int
        Index of this cell in distribcell arrays

    """

    def __init__(self, cell_id=None, name=''):
        # Initialize Cell class attributes
        self.id = cell_id
        self.name = name
        self._fill = None
        self._type = None
        self._region = None
        self._rotation = None
        self._translation = None
        self._offsets = None
        self._distribcell_index = None

    def __eq__(self, other):
        if not isinstance(other, Cell):
            return False
        elif self.id != other.id:
            return False
        elif self.name != other.name:
            return False
        elif self.fill != other.fill:
             return False
        elif self.region != other.region:
            return False
        elif self.rotation != other.rotation:
            return False
        elif self.translation != other.translation:
            return False
        else:
            return True

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(repr(self))

    def __repr__(self):
        string = 'Cell\n'
        string += '{0: <16}{1}{2}\n'.format('\tID', '=\t', self._id)
        string += '{0: <16}{1}{2}\n'.format('\tName', '=\t', self._name)

        if isinstance(self._fill, openmc.Material):
            string += '{0: <16}{1}{2}\n'.format('\tMaterial', '=\t',
                                                self._fill._id)
        elif isinstance(self._fill, Iterable):
            string += '{0: <16}{1}'.format('\tMaterial', '=\t')
            string += '['
            string += ', '.join(['void' if m == 'void' else str(m.id)
                                 for m in self.fill])
            string += ']\n'
        elif isinstance(self._fill, (Universe, Lattice)):
            string += '{0: <16}{1}{2}\n'.format('\tFill', '=\t',
                                                self._fill._id)
        else:
            string += '{0: <16}{1}{2}\n'.format('\tFill', '=\t', self._fill)

        string += '{0: <16}{1}{2}\n'.format('\tRegion', '=\t', self._region)

        string += '{0: <16}{1}{2}\n'.format('\tRotation', '=\t',
                                            self._rotation)
        string += '{0: <16}{1}{2}\n'.format('\tTranslation', '=\t',
                                            self._translation)
        string += '{0: <16}{1}{2}\n'.format('\tOffset', '=\t', self._offsets)
        string += '{0: <16}{1}{2}\n'.format('\tDistribcell index', '=\t',
                                            self._distribcell_index)

        return string

    @property
    def id(self):
        return self._id

    @property
    def name(self):
        return self._name

    @property
    def fill(self):
        return self._fill

    @property
    def fill_type(self):
        if isinstance(self.fill, openmc.Material):
            return 'material'
        elif isinstance(self.fill, openmc.Universe):
            return 'universe'
        elif isinstance(self.fill, openmc.Lattice):
            return 'lattice'
        else:
            return None

    @property
    def region(self):
        return self._region

    @property
    def rotation(self):
        return self._rotation

    @property
    def translation(self):
        return self._translation

    @property
    def offsets(self):
        return self._offsets

    @property
    def distribcell_index(self):
        return self._distribcell_index

    @id.setter
    def id(self, cell_id):
        if cell_id is None:
            global AUTO_CELL_ID
            self._id = AUTO_CELL_ID
            AUTO_CELL_ID += 1
        else:
            cv.check_type('cell ID', cell_id, Integral)
            cv.check_greater_than('cell ID', cell_id, 0, equality=True)
            self._id = cell_id

    @name.setter
    def name(self, name):
        if name is not None:
            cv.check_type('cell name', name, basestring)
            self._name = name
        else:
            self._name = ''

    @fill.setter
    def fill(self, fill):
        if isinstance(fill, basestring):
            if fill.strip().lower() == 'void':
                self._type = 'void'
            else:
                msg = 'Unable to set Cell ID="{0}" to use a non-Material or ' \
                       'Universe fill "{1}"'.format(self._id, fill)
                raise ValueError(msg)

        elif isinstance(fill, openmc.Material):
            self._type = 'normal'

        elif isinstance(fill, Iterable):
            cv.check_type('cell.fill', fill, Iterable,
                          (openmc.Material, basestring))
            self._type = 'normal'

        elif isinstance(fill, Universe):
            self._type = 'fill'

        elif isinstance(fill, Lattice):
            self._type = 'lattice'

        else:
            msg = 'Unable to set Cell ID="{0}" to use a non-Material or ' \
                   'Universe fill "{1}"'.format(self._id, fill)
            raise ValueError(msg)

        self._fill = fill

    @rotation.setter
    def rotation(self, rotation):
        cv.check_type('cell rotation', rotation, Iterable, Real)
        cv.check_length('cell rotation', rotation, 3)
        self._rotation = rotation

    @translation.setter
    def translation(self, translation):
        cv.check_type('cell translation', translation, Iterable, Real)
        cv.check_length('cell translation', translation, 3)
        self._translation = translation

    @offsets.setter
    def offsets(self, offsets):
        cv.check_type('cell offsets', offsets, Iterable)
        self._offsets = offsets

    @region.setter
    def region(self, region):
        cv.check_type('cell region', region, Region)
        self._region = region

    @distribcell_index.setter
    def distribcell_index(self, ind):
        cv.check_type('distribcell index', ind, Integral)
        self._distribcell_index = ind

    def add_surface(self, surface, halfspace):
        """Add a half-space to the list of half-spaces whose intersection defines the
        cell.

        .. deprecated:: 0.7.1
            Use the Cell.region property to directly specify a Region
            expression.

        Parameters
        ----------
        surface : openmc.surface.Surface
            Quadric surface dividing space
        halfspace : {-1, 1}
            Indicate whether the negative or positive half-space is to be used

        """

        warnings.warn("Cell.add_surface(...) has been deprecated and may be "
                      "removed in a future version. The region for a Cell "
                      "should be defined using the region property directly.",
                      DeprecationWarning)

        if not isinstance(surface, openmc.Surface):
            msg = 'Unable to add Surface "{0}" to Cell ID="{1}" since it is ' \
                        'not a Surface object'.format(surface, self._id)
            raise ValueError(msg)

        if halfspace not in [-1, +1]:
            msg = 'Unable to add Surface "{0}" to Cell ID="{1}" with halfspace ' \
                  '"{2}" since it is not +/-1'.format(surface, self._id, halfspace)
            raise ValueError(msg)

        # If no region has been assigned, simply use the half-space. Otherwise,
        # take the intersection of the current region and the half-space
        # specified
        region = +surface if halfspace == 1 else -surface
        if self.region is None:
            self.region = region
        else:
            if isinstance(self.region, Intersection):
                self.region.nodes.append(region)
            else:
                self.region = Intersection(self.region, region)

    def get_cell_instance(self, path, distribcell_index):

        # If the Cell is filled by a Material
        if self._type == 'normal' or self._type == 'void':
            offset = 0

        # If the Cell is filled by a Universe
        elif self._type == 'fill':
            offset = self.offsets[distribcell_index-1]
            offset += self.fill.get_cell_instance(path, distribcell_index)

        # If the Cell is filled by a Lattice
        else:
            offset = self.fill.get_cell_instance(path, distribcell_index)

        return offset

    def get_all_nuclides(self):
        """Return all nuclides contained in the cell

        Returns
        -------
        nuclides : dict
            Dictionary whose keys are nuclide names and values are 2-tuples of
            (nuclide, density)

        """

        nuclides = OrderedDict()

        if self._type != 'void':
            nuclides.update(self._fill.get_all_nuclides())

        return nuclides

    def get_all_cells(self):
        """Return all cells that are contained within this one if it is filled with a
        universe or lattice

        Returns
        -------
        cells : dict
            Dictionary whose keys are cell IDs and values are Cell instances

        """

        cells = OrderedDict()

        if self._type == 'fill' or self._type == 'lattice':
            cells.update(self._fill.get_all_cells())

        return cells

    def get_all_materials(self):
        """Return all materials that are contained within the cell

        Returns
        -------
        materials : dict
            Dictionary whose keys are material IDs and values are Material instances

        """

        materials = OrderedDict()
        if self.fill_type == 'material':
            materials[self.fill.id] = self.fill

        # Append all Cells in each Cell in the Universe to the dictionary
        cells = self.get_all_cells()
        for cell_id, cell in cells.items():
            materials.update(cell.get_all_materials())

        return materials

    def get_all_universes(self):
        """Return all universes that are contained within this one if any of
        its cells are filled with a universe or lattice.

        Returns
        -------
        universes : dict
            Dictionary whose keys are universe IDs and values are Universe
            instances

        """

        universes = OrderedDict()

        if self._type == 'fill':
            universes[self._fill._id] = self._fill
            universes.update(self._fill.get_all_universes())
        elif self._type == 'lattice':
            universes.update(self._fill.get_all_universes())

        return universes

    def create_xml_subelement(self, xml_element):
        element = ET.Element("cell")
        element.set("id", str(self.id))

        if len(self._name) > 0:
            element.set("name", str(self.name))

        if isinstance(self.fill, basestring):
            element.set("material", "void")

        elif isinstance(self.fill, openmc.Material):
            element.set("material", str(self.fill.id))

        elif isinstance(self.fill, Iterable):
            element.set("material", ' '.join([m if m == 'void' else str(m.id)
                                              for m in self.fill]))

        elif isinstance(self.fill, (Universe, Lattice)):
            element.set("fill", str(self.fill.id))
            self.fill.create_xml_subelement(xml_element)

        else:
            element.set("fill", str(self.fill))
            self.fill.create_xml_subelement(xml_element)

        if self.region is not None:
            # Set the region attribute with the region specification
            element.set("region", str(self.region))

            # Only surfaces that appear in a region are added to the geometry
            # file, so the appropriate check is performed here. First we create
            # a function which is called recursively to navigate through the CSG
            # tree. When it reaches a leaf (a Halfspace), it creates a <surface>
            # element for the corresponding surface if none has been created
            # thus far.
            def create_surface_elements(node, element):
                if isinstance(node, Halfspace):
                    path = './surface[@id=\'{0}\']'.format(node.surface.id)
                    if xml_element.find(path) is None:
                        surface_subelement = node.surface.create_xml_subelement()
                        xml_element.append(surface_subelement)
                elif isinstance(node, Complement):
                    create_surface_elements(node.node, element)
                else:
                    for subnode in node.nodes:
                        create_surface_elements(subnode, element)

            # Call the recursive function from the top node
            create_surface_elements(self.region, xml_element)

        if self.translation is not None:
            element.set("translation", ' '.join(map(str, self.translation)))

        if self.rotation is not None:
            element.set("rotation", ' '.join(map(str, self.rotation)))

        return element


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
    cells : dict
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
        cell : Cell
            Cell to add

        """

        if not isinstance(cell, Cell):
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
        cells : array-like of Cell
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
        cell : Cell
            Cell to remove

        """

        if not isinstance(cell, Cell):
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
        nuclides : dict
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
        cells : dict
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
        materials : dict
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
        universes : dict
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


class Lattice(object):
    """A repeating structure wherein each element is a universe.

    Parameters
    ----------
    lattice_id : int, optional
        Unique identifier for the lattice. If not specified, an identifier will
        automatically be assigned.
    name : str, optional
        Name of the lattice. If not specified, the name is the empty string.

    Attributes
    ----------
    id : int
        Unique identifier for the lattice
    name : str
        Name of the lattice
    pitch : float
        Pitch of the lattice in cm
    outer : int
        The unique identifier of a universe to fill all space outside the
        lattice
    universes : ndarray of Universe
        An array of universes filling each element of the lattice

    """

    # This is an abstract class which cannot be instantiated
    __metaclass__ = abc.ABCMeta

    def __init__(self, lattice_id=None, name=''):
        # Initialize Lattice class attributes
        self.id = lattice_id
        self.name = name
        self._pitch = None
        self._outer = None
        self._universes = None

    def __eq__(self, other):
        if not isinstance(other, Lattice):
            return False
        elif self.id != other.id:
            return False
        elif self.name != other.name:
            return False
        elif self.pitch != other.pitch:
            return False
        elif self.outer != other.outer:
            return False
        elif self.universes != other.universes:
            return False
        else:
            return True

    def __ne__(self, other):
        return not self == other

    @property
    def id(self):
        return self._id

    @property
    def name(self):
        return self._name

    @property
    def pitch(self):
        return self._pitch

    @property
    def outer(self):
        return self._outer

    @property
    def universes(self):
        return self._universes

    @id.setter
    def id(self, lattice_id):
        if lattice_id is None:
            global AUTO_UNIVERSE_ID
            self._id = AUTO_UNIVERSE_ID
            AUTO_UNIVERSE_ID += 1
        else:
            cv.check_type('lattice ID', lattice_id, Integral)
            cv.check_greater_than('lattice ID', lattice_id, 0, equality=True)
            self._id = lattice_id

    @name.setter
    def name(self, name):
        if name is not None:
            cv.check_type('lattice name', name, basestring)
            self._name = name
        else:
            self._name = ''

    @outer.setter
    def outer(self, outer):
        cv.check_type('outer universe', outer, Universe)
        self._outer = outer

    @universes.setter
    def universes(self, universes):
        cv.check_iterable_type('lattice universes', universes, Universe,
                               min_depth=2, max_depth=3)
        self._universes = np.asarray(universes)

    def get_unique_universes(self):
        """Determine all unique universes in the lattice

        Returns
        -------
        universes : dict
            Dictionary whose keys are universe IDs and values are Universe
            instances

        """

        univs = OrderedDict()
        for k in range(len(self._universes)):
            for j in range(len(self._universes[k])):
                if isinstance(self._universes[k][j], Universe):
                    u = self._universes[k][j]
                    univs[u._id] = u
                else:
                    for i in range(len(self._universes[k][j])):
                        u = self._universes[k][j][i]
                        assert isinstance(u, Universe)
                        univs[u._id] = u

        if self.outer is not None:
            univs[self.outer._id] = self.outer

        return univs

    def get_all_nuclides(self):
        """Return all nuclides contained in the lattice

        Returns
        -------
        nuclides : dict
            Dictionary whose keys are nuclide names and values are 2-tuples of
            (nuclide, density)

        """

        nuclides = OrderedDict()

        # Get all unique Universes contained in each of the lattice cells
        unique_universes = self.get_unique_universes()

        # Append all Universes containing each cell to the dictionary
        for universe_id, universe in unique_universes.items():
            nuclides.update(universe.get_all_nuclides())

        return nuclides

    def get_all_cells(self):
        """Return all cells that are contained within the lattice

        Returns
        -------
        cells : dict
            Dictionary whose keys are cell IDs and values are Cell instances

        """

        cells = OrderedDict()
        unique_universes = self.get_unique_universes()

        for universe_id, universe in unique_universes.items():
            cells.update(universe.get_all_cells())

        return cells

    def get_all_materials(self):
        """Return all materials that are contained within the lattice

        Returns
        -------
        materials : dict
            Dictionary whose keys are material IDs and values are Material instances

        """

        materials = OrderedDict()

        # Append all Cells in each Cell in the Universe to the dictionary
        cells = self.get_all_cells()
        for cell_id, cell in cells.items():
            materials.update(cell.get_all_materials())

        return materials

    def get_all_universes(self):
        """Return all universes that are contained within the lattice

        Returns
        -------
        universes : dict
            Dictionary whose keys are universe IDs and values are Universe
            instances

        """

        # Initialize a dictionary of all Universes contained by the Lattice
        # in each nested Universe level
        all_universes = OrderedDict()

        # Get all unique Universes contained in each of the lattice cells
        unique_universes = self.get_unique_universes()

        # Add the unique Universes filling each Lattice cell
        all_universes.update(unique_universes)

        # Append all Universes containing each cell to the dictionary
        for universe_id, universe in unique_universes.items():
            all_universes.update(universe.get_all_universes())

        return all_universes


class RectLattice(Lattice):
    """A lattice consisting of rectangular prisms.

    Parameters
    ----------
    lattice_id : int, optional
        Unique identifier for the lattice. If not specified, an identifier will
        automatically be assigned.
    name : str, optional
        Name of the lattice. If not specified, the name is the empty string.

    Attributes
    ----------
    id : int
        Unique identifier for the lattice
    name : str
        Name of the lattice
    dimension : array-like of int
        An array of two or three integers representing the number of lattice
        cells in the x- and y- (and z-) directions, respectively.
    lower_left : array-like of float
        The coordinates of the lower-left corner of the lattice. If the lattice
        is two-dimensional, only the x- and y-coordinates are specified.

    """

    def __init__(self, lattice_id=None, name=''):
        super(RectLattice, self).__init__(lattice_id, name)

        # Initialize Lattice class attributes
        self._dimension = None
        self._lower_left = None
        self._offsets = None

    def __eq__(self, other):
        if not isinstance(other, RectLattice):
            return False
        elif not super(RectLattice, self).__eq__(other):
            return False
        elif self.dimension != other.dimension:
            return False
        elif self.lower_left != other.lower_left:
            return False
        else:
            return True

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(repr(self))

    def __repr__(self):
        string = 'RectLattice\n'
        string += '{0: <16}{1}{2}\n'.format('\tID', '=\t', self._id)
        string += '{0: <16}{1}{2}\n'.format('\tName', '=\t', self._name)
        string += '{0: <16}{1}{2}\n'.format('\tDimension', '=\t',
                                            self._dimension)
        string += '{0: <16}{1}{2}\n'.format('\tLower Left', '=\t',
                                            self._lower_left)
        string += '{0: <16}{1}{2}\n'.format('\tPitch', '=\t', self._pitch)

        if self._outer is not None:
            string += '{0: <16}{1}{2}\n'.format('\tOuter', '=\t',
                                                self._outer._id)
        else:
            string += '{0: <16}{1}{2}\n'.format('\tOuter', '=\t',
                                                self._outer)

        string += '{0: <16}\n'.format('\tUniverses')

        # Lattice nested Universe IDs - column major for Fortran
        for i, universe in enumerate(np.ravel(self._universes)):
            string += '{0} '.format(universe._id)

            # Add a newline character every time we reach end of row of cells
            if (i+1) % self._dimension[-1] == 0:
                string += '\n'

        string = string.rstrip('\n')

        if self._offsets is not None:
            string += '{0: <16}\n'.format('\tOffsets')

            # Lattice cell offsets
            for i, offset in enumerate(np.ravel(self._offsets)):
                string += '{0} '.format(offset)

                # Add a newline character when we reach end of row of cells
                if (i+1) % self._dimension[-1] == 0:
                    string += '\n'

            string = string.rstrip('\n')

        return string

    @property
    def dimension(self):
        return self._dimension

    @property
    def lower_left(self):
        return self._lower_left

    @property
    def offsets(self):
        return self._offsets

    @dimension.setter
    def dimension(self, dimension):
        cv.check_type('lattice dimension', dimension, Iterable, Integral)
        cv.check_length('lattice dimension', dimension, 2, 3)
        for dim in dimension:
            cv.check_greater_than('lattice dimension', dim, 0)
        self._dimension = dimension

    @lower_left.setter
    def lower_left(self, lower_left):
        cv.check_type('lattice lower left corner', lower_left, Iterable, Real)
        cv.check_length('lattice lower left corner', lower_left, 2, 3)
        self._lower_left = lower_left

    @offsets.setter
    def offsets(self, offsets):
        cv.check_type('lattice offsets', offsets, Iterable)
        self._offsets = offsets

    @Lattice.pitch.setter
    def pitch(self, pitch):
        cv.check_type('lattice pitch', pitch, Iterable, Real)
        cv.check_length('lattice pitch', pitch, 2, 3)
        for dim in pitch:
            cv.check_greater_than('lattice pitch', dim, 0.0)
        self._pitch = pitch

    def get_cell_instance(self, path, distribcell_index):

        # Extract the lattice element from the path
        next_index = path.index('-')
        lat_id_indices = path[:next_index]
        path = path[next_index+2:]

        # Extract the lattice cell indices from the path
        i1 = lat_id_indices.index('(')
        i2 = lat_id_indices.index(')')
        i = lat_id_indices[i1+1:i2]
        lat_x = int(i.split(',')[0]) - 1
        lat_y = int(i.split(',')[1]) - 1
        lat_z = int(i.split(',')[2]) - 1

        # For 2D Lattices
        if len(self._dimension) == 2:
            offset = self._offsets[lat_z, lat_y, lat_x, distribcell_index-1]
            offset += self._universes[lat_x][lat_y].get_cell_instance(path,
                                                              distribcell_index)

        # For 3D Lattices
        else:
            offset = self._offsets[lat_z, lat_y, lat_x, distribcell_index-1]
            offset += self._universes[lat_z][lat_y][lat_x].get_cell_instance(
                                                        path, distribcell_index)

        return offset

    def create_xml_subelement(self, xml_element):

        # Determine if XML element already contains subelement for this Lattice
        path = './lattice[@id=\'{0}\']'.format(self._id)
        test = xml_element.find(path)

        # If the element does contain the Lattice subelement, then return
        if test is not None:
            return

        lattice_subelement = ET.Element("lattice")
        lattice_subelement.set("id", str(self._id))

        if len(self._name) > 0:
            lattice_subelement.set("name", str(self._name))

        # Export the Lattice cell pitch
        pitch = ET.SubElement(lattice_subelement, "pitch")
        pitch.text = ' '.join(map(str, self._pitch))

        # Export the Lattice outer Universe (if specified)
        if self._outer is not None:
            outer = ET.SubElement(lattice_subelement, "outer")
            outer.text = '{0}'.format(self._outer._id)
            self._outer.create_xml_subelement(xml_element)

        # Export Lattice cell dimensions
        dimension = ET.SubElement(lattice_subelement, "dimension")
        dimension.text = ' '.join(map(str, self._dimension))

        # Export Lattice lower left
        lower_left = ET.SubElement(lattice_subelement, "lower_left")
        lower_left.text = ' '.join(map(str, self._lower_left))

        # Export the Lattice nested Universe IDs - column major for Fortran
        universe_ids = '\n'

        # 3D Lattices
        if len(self._dimension) == 3:
            for z in range(self._dimension[2]):
                for y in range(self._dimension[1]):
                    for x in range(self._dimension[0]):
                        universe = self._universes[z][y][x]

                        # Append Universe ID to the Lattice XML subelement
                        universe_ids += '{0} '.format(universe._id)

                        # Create XML subelement for this Universe
                        universe.create_xml_subelement(xml_element)

                    # Add newline character when we reach end of row of cells
                    universe_ids += '\n'

                # Add newline character when we reach end of row of cells
                universe_ids += '\n'

        # 2D Lattices
        else:
            for y in range(self._dimension[1]):
                for x in range(self._dimension[0]):
                    universe = self._universes[y][x]

                    # Append Universe ID to Lattice XML subelement
                    universe_ids += '{0} '.format(universe._id)

                    # Create XML subelement for this Universe
                    universe.create_xml_subelement(xml_element)

                # Add newline character when we reach end of row of cells
                universe_ids += '\n'

        # Remove trailing newline character from Universe IDs string
        universe_ids = universe_ids.rstrip('\n')

        universes = ET.SubElement(lattice_subelement, "universes")
        universes.text = universe_ids

        # Append the XML subelement for this Lattice to the XML element
        xml_element.append(lattice_subelement)


class HexLattice(Lattice):
    """A lattice consisting of hexagonal prisms.

    Parameters
    ----------
    lattice_id : int, optional
        Unique identifier for the lattice. If not specified, an identifier will
        automatically be assigned.
    name : str, optional
        Name of the lattice. If not specified, the name is the empty string.

    Attributes
    ----------
    id : int
        Unique identifier for the lattice
    name : str
        Name of the lattice
    num_rings : int
        Number of radial ring positions in the xy-plane
    num_axial : int
        Number of positions along the z-axis.
    center : array-like of float
        Coordinates of the center of the lattice. If the lattice does not have
        axial sections then only the x- and y-coordinates are specified

    """

    def __init__(self, lattice_id=None, name=''):
        super(HexLattice, self).__init__(lattice_id, name)

        # Initialize Lattice class attributes
        self._num_rings = None
        self._num_axial = None
        self._center = None

    def __eq__(self, other):
        if not isinstance(other, HexLattice):
            return False
        elif not super(HexLattice, self).__eq__(other):
            return False
        elif self.num_rings != other.num_rings:
            return False
        elif self.num_axial != other.num_axial:
            return False
        elif self.center != other.center:
            return False
        else:
            return True

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(repr(self))

    def __repr__(self):
        string = 'HexLattice\n'
        string += '{0: <16}{1}{2}\n'.format('\tID', '=\t', self._id)
        string += '{0: <16}{1}{2}\n'.format('\tName', '=\t', self._name)
        string += '{0: <16}{1}{2}\n'.format('\t# Rings', '=\t', self._num_rings)
        string += '{0: <16}{1}{2}\n'.format('\t# Axial', '=\t', self._num_axial)
        string += '{0: <16}{1}{2}\n'.format('\tCenter', '=\t',
                                            self._center)
        string += '{0: <16}{1}{2}\n'.format('\tPitch', '=\t', self._pitch)

        if self._outer is not None:
            string += '{0: <16}{1}{2}\n'.format('\tOuter', '=\t',
                                                self._outer._id)
        else:
            string += '{0: <16}{1}{2}\n'.format('\tOuter', '=\t',
                                                self._outer)

        string += '{0: <16}\n'.format('\tUniverses')

        if self._num_axial is not None:
            slices = [self._repr_axial_slice(x) for x in self._universes]
            string += '\n'.join(slices)

        else:
            string += self._repr_axial_slice(self._universes)

        return string

    @property
    def num_rings(self):
        return self._num_rings

    @property
    def num_axial(self):
        return self._num_axial

    @property
    def center(self):
        return self._center

    @num_rings.setter
    def num_rings(self, num_rings):
        cv.check_type('number of rings', num_rings, Integral)
        cv.check_greater_than('number of rings', num_rings, 0)
        self._num_rings = num_rings

    @num_axial.setter
    def num_axial(self, num_axial):
        cv.check_type('number of axial', num_axial, Integral)
        cv.check_greater_than('number of axial', num_axial, 0)
        self._num_axial = num_axial

    @center.setter
    def center(self, center):
        cv.check_type('lattice center', center, Iterable, Real)
        cv.check_length('lattice center', center, 2, 3)
        self._center = center

    @Lattice.pitch.setter
    def pitch(self, pitch):
        cv.check_type('lattice pitch', pitch, Iterable, Real)
        cv.check_length('lattice pitch', pitch, 1, 2)
        for dim in pitch:
            cv.check_greater_than('lattice pitch', dim, 0)
        self._pitch = pitch

    @Lattice.universes.setter
    def universes(self, universes):
        # Call Lattice.universes parent class setter property
        Lattice.universes.fset(self, universes)

        # NOTE: This routine assumes that the user creates a "ragged" list of
        # lists, where each sub-list corresponds to one ring of Universes.
        # The sub-lists are ordered from outermost ring to innermost ring.
        # The Universes within each sub-list are ordered from the "top" in a
        # clockwise fashion.

        # Check to see if the given universes look like a 2D or a 3D array.
        if isinstance(self._universes[0][0], Universe):
            n_dims = 2

        elif isinstance(self._universes[0][0][0], Universe):
            n_dims = 3

        else:
            msg = 'HexLattice ID={0:d} does not appear to be either 2D or ' \
                  '3D.  Make sure set_universes was given a two-deep or ' \
                  'three-deep iterable of universes.'.format(self._id)
            raise RuntimeError(msg)

        # Set the number of axial positions.
        if n_dims == 3:
            self.num_axial = len(self._universes)
        else:
            self._num_axial = None

        # Set the number of rings and make sure this number is consistent for
        # all axial positions.
        if n_dims == 3:
            self.num_rings = len(self._universes)
            for rings in self._universes:
                if len(rings) != self._num_rings:
                    msg = 'HexLattice ID={0:d} has an inconsistent number of ' \
                          'rings per axial positon'.format(self._id)
                    raise ValueError(msg)

        else:
            self.num_rings = len(self._universes)

        # Make sure there are the correct number of elements in each ring.
        if n_dims == 3:
            for axial_slice in self._universes:
                # Check the center ring.
                if len(axial_slice[-1]) != 1:
                    msg = 'HexLattice ID={0:d} has the wrong number of ' \
                          'elements in the innermost ring.  Only 1 element is ' \
                          'allowed in the innermost ring.'.format(self._id)
                    raise ValueError(msg)

                # Check the outer rings.
                for r in range(self._num_rings-1):
                    if len(axial_slice[r]) != 6*(self._num_rings - 1 - r):
                        msg = 'HexLattice ID={0:d} has the wrong number of ' \
                              'elements in ring number {1:d} (counting from the '\
                              'outermost ring).  This ring should have {2:d} ' \
                              'elements.'.format(self._id, r,
                                                 6*(self._num_rings - 1 - r))
                        raise ValueError(msg)

        else:
            axial_slice = self._universes
            # Check the center ring.
            if len(axial_slice[-1]) != 1:
                msg = 'HexLattice ID={0:d} has the wrong number of ' \
                      'elements in the innermost ring.  Only 1 element is ' \
                      'allowed in the innermost ring.'.format(self._id)
                raise ValueError(msg)

            # Check the outer rings.
            for r in range(self._num_rings-1):
                if len(axial_slice[r]) != 6*(self._num_rings - 1 - r):
                    msg = 'HexLattice ID={0:d} has the wrong number of ' \
                          'elements in ring number {1:d} (counting from the '\
                          'outermost ring).  This ring should have {2:d} ' \
                          'elements.'.format(self._id, r,
                                             6*(self._num_rings - 1 - r))
                    raise ValueError(msg)

    def create_xml_subelement(self, xml_element):
        # Determine if XML element already contains subelement for this Lattice
        path = './hex_lattice[@id=\'{0}\']'.format(self._id)
        test = xml_element.find(path)

        # If the element does contain the Lattice subelement, then return
        if test is not None:
            return

        lattice_subelement = ET.Element("hex_lattice")
        lattice_subelement.set("id", str(self._id))

        if len(self._name) > 0:
            lattice_subelement.set("name", str(self._name))

        # Export the Lattice cell pitch
        pitch = ET.SubElement(lattice_subelement, "pitch")
        pitch.text = ' '.join(map(str, self._pitch))

        # Export the Lattice outer Universe (if specified)
        if self._outer is not None:
            outer = ET.SubElement(lattice_subelement, "outer")
            outer.text = '{0}'.format(self._outer._id)
            self._outer.create_xml_subelement(xml_element)

        lattice_subelement.set("n_rings", str(self._num_rings))

        if self._num_axial is not None:
            lattice_subelement.set("n_axial", str(self._num_axial))

        # Export Lattice cell center
        dimension = ET.SubElement(lattice_subelement, "center")
        dimension.text = ' '.join(map(str, self._center))

        # Export the Lattice nested Universe IDs.

        # 3D Lattices
        if self._num_axial is not None:
            slices = []
            for z in range(self._num_axial):
                # Initialize the center universe.
                universe = self._universes[z][-1][0]
                universe.create_xml_subelement(xml_element)

                # Initialize the remaining universes.
                for r in range(self._num_rings-1):
                    for theta in range(6*(self._num_rings - 1 - r)):
                        universe = self._universes[z][r][theta]
                        universe.create_xml_subelement(xml_element)

                # Get a string representation of the universe IDs.
                slices.append(self._repr_axial_slice(self._universes[z]))

            # Collapse the list of axial slices into a single string.
            universe_ids = '\n'.join(slices)

        # 2D Lattices
        else:
            # Initialize the center universe.
            universe = self._universes[-1][0]
            universe.create_xml_subelement(xml_element)

            # Initialize the remaining universes.
            for r in range(self._num_rings - 1):
                for theta in range(6*(self._num_rings - 1 - r)):
                    universe = self._universes[r][theta]
                    universe.create_xml_subelement(xml_element)

            # Get a string representation of the universe IDs.
            universe_ids = self._repr_axial_slice(self._universes)

        universes = ET.SubElement(lattice_subelement, "universes")
        universes.text = '\n' + universe_ids

        # Append the XML subelement for this Lattice to the XML element
        xml_element.append(lattice_subelement)

    def _repr_axial_slice(self, universes):
        """Return string representation for the given 2D group of universes.

        The 'universes' argument should be a list of lists of universes where
        each sub-list represents a single ring.  The first list should be the
        outer ring.
        """

        # Find the largest universe ID and count the number of digits so we can
        # properly pad the output string later.
        largest_id = max([max([univ._id for univ in ring])
                          for ring in universes])
        n_digits = len(str(largest_id))
        pad = ' '*n_digits
        id_form = '{: ^' + str(n_digits) + 'd}'

        # Initialize the list for each row.
        rows = [ [] for i in range(1 + 4 * (self._num_rings-1)) ]
        middle = 2 * (self._num_rings - 1)

        # Start with the degenerate first ring.
        universe = universes[-1][0]
        rows[middle] = [id_form.format(universe._id)]

        # Add universes one ring at a time.
        for r in range(1, self._num_rings):
            # r_prime increments down while r increments up.
            r_prime = self._num_rings - 1 - r
            theta = 0
            y = middle + 2*r

            # Climb down the top-right.
            for i in range(r):
                # Add the universe.
                universe = universes[r_prime][theta]
                rows[y].append(id_form.format(universe._id))

                # Translate the indices.
                y -= 1
                theta += 1

            # Climb down the right.
            for i in range(r):
                # Add the universe.
                universe = universes[r_prime][theta]
                rows[y].append(id_form.format(universe._id))

                # Translate the indices.
                y -= 2
                theta += 1

            # Climb down the bottom-right.
            for i in range(r):
                # Add the universe.
                universe = universes[r_prime][theta]
                rows[y].append(id_form.format(universe._id))

                # Translate the indices.
                y -= 1
                theta += 1

            # Climb up the bottom-left.
            for i in range(r):
                # Add the universe.
                universe = universes[r_prime][theta]
                rows[y].insert(0, id_form.format(universe._id))

                # Translate the indices.
                y += 1
                theta += 1

            # Climb up the left.
            for i in range(r):
                # Add the universe.
                universe = universes[r_prime][theta]
                rows[y].insert(0, id_form.format(universe._id))

                # Translate the indices.
                y += 2
                theta += 1

            # Climb up the top-left.
            for i in range(r):
                # Add the universe.
                universe = universes[r_prime][theta]
                rows[y].insert(0, id_form.format(universe._id))

                # Translate the indices.
                y += 1
                theta += 1

        # Flip the rows and join each row into a single string.
        rows = [pad.join(x) for x in rows[::-1]]

        # Pad the beginning of the rows so they line up properly.
        for y in range(self._num_rings - 1):
            rows[y] = (self._num_rings - 1 - y)*pad + rows[y]
            rows[-1 - y] = (self._num_rings - 1 - y)*pad + rows[-1 - y]

        for y in range(self._num_rings % 2, self._num_rings, 2):
            rows[middle + y] = pad + rows[middle + y]
            if y != 0:
                rows[middle - y] = pad + rows[middle - y]

        # Join the rows together and return the string.
        universe_ids = '\n'.join(rows)
        return universe_ids
