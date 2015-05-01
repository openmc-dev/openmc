import abc
from collections import OrderedDict
from xml.etree import ElementTree as ET

import openmc
from openmc.checkvalue import *


################################################################################
####################################  Cell  ####################################
################################################################################

# A static variable for auto-generated Cell IDs
AUTO_CELL_ID = 10000

# A dictionary for storing IDs of cell elements that have already been written,
# used to optimize the writing process
WRITTEN_IDS = {}

def reset_auto_cell_id():
    global AUTO_CELL_ID
    AUTO_CELL_ID = 10000


class Cell(object):

    def __init__(self, cell_id=None, name=''):

        # Initialize Cell class attributes
        self.id = cell_id
        self.name= name
        self._fill = None
        self._type = None
        self._surfaces = {}
        self._rotation = None
        self._translation = None
        self._offset = None


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
    def type(self):
        return self._fill


    @property
    def surfaces(self):
        return self._surfaces


    @property
    def rotation(self):
        return self._rotation


    @property
    def translation(self):
        return self._translation


    @property
    def offset(self):
        return self._offset


    @id.setter
    def id(self, cell_id):

        if cell_id is None:
            global AUTO_CELL_ID
            self._id = AUTO_CELL_ID
            AUTO_CELL_ID += 1

        # Check that the ID is an integer and wasn't already used
        elif not is_integer(cell_id):
             msg = 'Unable to set a non-integer Cell ID {0}'.format(cell_id)
             raise ValueError(msg)

        elif cell_id < 0:
            msg = 'Unable to set Cell ID to {0} since it must be a ' \
                        'non-negative integer'.format(cell_id)
            raise ValueError(msg)

        else:
            self._id = cell_id


    @name.setter
    def name(self, name):

        if not isinstance(name, str):
            msg = 'Unable to set name for Cell ID={0} with a non-string ' \
                        'value {1}'.format(self._id, name)
            raise ValueError(msg)

        else:
            self._name = name


    @fill.setter
    def fill(self, fill):

        if isinstance(fill, str):
            if fill.strip().lower() == 'void':
                self._type = 'void'
            else:
                msg = 'Unable to set Cell ID={0} to use a non-Material or ' \
                       'Universe fill {1}'.format(self._id, fill)
                raise ValueError(msg)

        elif isinstance(fill, openmc.Material):
            self._type = 'normal'

        elif isinstance(fill, Universe):
            self._type = 'fill'

        elif isinstance(fill, Lattice):
            self._type = 'lattice'

        else:
            msg = 'Unable to set Cell ID={0} to use a non-Material or ' \
                   'Universe fill {1}'.format(self._id, fill)
            raise ValueError(msg)

        self._fill = fill


    @rotation.setter
    def rotation(self, rotation):

        if not isinstance(rotation, (tuple, list, np.ndarray)):
            msg = 'Unable to add rotation {0} to Cell ID={1} since ' \
                  'it is not a Python tuple/list or NumPy ' \
                  'array'.format(rotation, self._id)
            raise ValueError(msg)

        elif len(rotation) != 3:
            msg = 'Unable to add rotation {0} to Cell ID={1} since ' \
                  'it does not contain 3 values'.format(rotation, self._id)
            raise ValueError(msg)

        for axis in rotation:

            if not is_integer(axis) and not is_float(axis):
                msg = 'Unable to add rotation {0} to Cell ID={1} since ' \
                      'it is not an integer or floating point ' \
                      'value'.format(axis, self._id)
                raise ValueError(msg)

            self._rotation = rotation


    @translation.setter
    def translation(self, translation):

        if not isinstance(translation, (tuple, list, np.ndarray)):
            msg = 'Unable to add translation {0} to Cell ID={1} since ' \
                  'it is not a Python tuple/list or NumPy ' \
                  'array'.format(translation, self._id)
            raise ValueError(msg)

        elif len(translation) != 3:
            msg = 'Unable to add translation {0} to Cell ID={1} since ' \
                  'it does not contain 3 values'.format(translation, self._id)
            raise ValueError(msg)

        for axis in translation:
            if not is_integer(axis) and not is_float(axis):
                msg = 'Unable to add translation {0} to Cell ID={1} since ' \
                      'it is not an integer or floating point ' \
                      'value'.format(axis, self._id)
                raise ValueError(msg)

        self._translation = translation


    @offset.setter
    def offset(self, offset):

        if not isinstance(offset, (tuple, list, np.ndarray)):
            msg = 'Unable to set offset {0} to Cell ID={1} since ' \
                  'it is not a Python tuple/list or NumPy ' \
                  'array'.format(offset, self._id)
            raise ValueError(msg)

        self._offset = offset


    def add_surface(self, surface, halfspace):

        if not isinstance(surface, openmc.Surface):
            msg = 'Unable to add Surface {0} to Cell ID={1} since it is ' \
                        'not a Surface object'.format(surface, self._id)
            raise ValueError(msg)

        if not halfspace in [-1, +1]:
            msg = 'Unable to add Surface {0} to Cell ID={1} with halfspace ' \
                  '{2} since it is not +/-1'.format(surface, self._id, halfspace)
            raise ValueError(msg)

        # If the Cell does not already contain the Surface, add it
        if not surface._id in self._surfaces:
            self._surfaces[surface._id] = (surface, halfspace)


    def remove_surface(self, surface):

        if not isinstance(surface, openmc.Surface):
            msg = 'Unable to remove Surface {0} from Cell ID={1} since it is ' \
                        'not a Surface object'.format(surface, self._id)
            raise ValueError(msg)

        # If the Cell contains the Surface, delete it
        if surface._id in self._surfaces:
            del self._surfaces[surface._id]


    def set_surfaces(self, surfaces, halfspaces):

        if not isinstance(surfaces, (tuple, list, np.ndarray)):
            msg = 'Unable to set Cell ID={0} with Surfaces {1} since ' \
                  'it is not a a Python tuple/list or NumPy ' \
                  'array'.format(self._id, surfaces)
            raise ValueError(msg)

        if not isinstance(halfspaces, (tuple, list, np.ndarray)):
            msg = 'Unable to set Cell ID={0} with Surface halfspaces {1} ' \
                        'since it is not a Python tuple/list or NumPy ' \
                        'array'.format(self._id, halfspaces)
            raise ValueError(msg)

        for surface in surfaces:
            self.add_surface(surface)


    def get_offset(self, path, filter_offset):

        # Get the current element and remove it from the list
        cell_id = path[0]
        path = path[1:]

        # If the Cell is filled by a Material
        if self._type == 'normal' or self._type == 'void':
            offset = 0

        # If the Cell is filled by a Universe
        elif self._type == 'fill':
            offset = self._offset[filter_offset-1]
            offset += self._fill.get_offset(path, filter_offset)

        # If the Cell is filled by a Lattice
        else:
            offset = self._fill.get_offset(path, filter_offset)

        return offset

        # Make a recursive call to the Universe filling this Cell
        offset = self._cells[cell_id].get_offset(path, filter_offset)

        # Return the offset computed at all nested Universe levels
        return offset


    def get_all_nuclides(self):

        nuclides = {}

        if self._type != 'void':
            nuclides.update(self._fill.get_all_nuclides())

        return nuclides


    def get_all_cells(self):

        cells = {}

        if self._type == 'fill' or self._type == 'lattice':
            cells.update(self._fill.get_all_cells())

        return cells


    def get_all_universes(self):

        universes = {}

        if self._type == 'fill':
            universes[self._fill._id] = self._fill
            universes.update(self._fill.get_all_universes())
        elif self._type == 'lattice':
            universes.update(self._fill.get_all_universes())

        return universes


    def __repr__(self):

        string = 'Cell\n'
        string += '{0: <16}{1}{2}\n'.format('\tID', '=\t', self._id)
        string += '{0: <16}{1}{2}\n'.format('\tName', '=\t', self._name)

        if isinstance(self._fill, openmc.Material):
            string += '{0: <16}{1}{2}\n'.format('\tMaterial', '=\t',
                                                self._fill._id)
        elif isinstance(self._fill, (Universe, Lattice)):
            string += '{0: <16}{1}{2}\n'.format('\tFill', '=\t',
                                                self._fill._id)
        else:
            string += '{0: <16}{1}{2}\n'.format('\tFill', '=\t', self._fill)

        string += '{0: <16}{1}\n'.format('\tSurfaces', '=\t')

        for surface_id in self._surfaces:
            halfspace = self._surfaces[surface_id][1]
            string += '{0} '.format(halfspace * surface_id)

        string = string.rstrip(' ') + '\n'

        string += '{0: <16}{1}{2}\n'.format('\tRotation', '=\t',
                                            self._rotation)
        string += '{0: <16}{1}{2}\n'.format('\tTranslation', '=\t',
                                            self._translation)
        string += '{0: <16}{1}{2}\n'.format('\tOffset', '=\t', self._offset)

        return string


    def create_xml_subelement(self, xml_element):

        element = ET.Element("cell")
        element.set("id", str(self._id))

        if len(self._name) > 0:
            element.set("name", str(self._name))

        if isinstance(self._fill, openmc.Material):
            element.set("material", str(self._fill._id))

        elif isinstance(self._fill, (Universe, Lattice)):
            element.set("fill", str(self._fill._id))
            self._fill.create_xml_subelement(xml_element)

        elif self._fill.strip().lower() == "void":
            element.set("material", "void")

        else:
            element.set("fill", str(self._fill))
            self._fill.create_xml_subelement(xml_element)


        if not self._surfaces is None:

            surfaces = ''

            for surface_id in self._surfaces:

                # Determine if XML element already includes this Surface
                path = './surface[@id=\'{0}\']'.format(surface_id)
                test = xml_element.find(path)

                # If the element does not contain the Surface subelement
                if test is None:

                    # Create the XML subelement for this Surface
                    surface = self._surfaces[surface_id][0]
                    surface_subelement = surface.create_xml_subelement()
                    xml_element.append(surface_subelement)

                # Append the halfspace and Surface ID
                halfspace = self._surfaces[surface_id][1]
                surfaces += '{0} '.format(halfspace * surface_id)

            element.set("surfaces", surfaces.rstrip(' '))


        if not self._translation is None:

            translation = ''

            for axis in self._translation:
                translation += '{0} '.format(axis)

            element.set("translation", translation.rstrip(' '))


        if not self._rotation is None:

            rotation = ''

            for axis in self._rotation:
                rotation += '{0} '.format(axis)

            element.set("rotation", rotation.rstrip(' '))

        return element



################################################################################
###################################  Universe  #################################
################################################################################

# A static variable for auto-generated Lattice (Universe) IDs
AUTO_UNIVERSE_ID = 10000

def reset_auto_universe_id():
    global AUTO_UNIVERSE_ID
    AUTO_UNIVERSE_ID = 10000


class Universe(object):

    def __init__(self, universe_id=None, name=''):

        # Initialize Cell class attributes
        self.id = universe_id
        self.name = name

        # Keys     - Cell IDs
        # Values - Cells
        self._cells = {}

        # Keys     - Cell IDs
        # Values - Offsets
        self._cell_offsets = OrderedDict()
        self._num_regions = 0


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

        # Check that the ID is an integer and wasn't already used
        elif not is_integer(universe_id):
            msg = 'Unable to set Universe ID to a non-integer ' \
                  '{0}'.format(universe_id)
            raise ValueError(msg)

        elif universe_id < 0:
            msg = 'Unable to set Universe ID to {0} since it must be a ' \
                  'non-negative integer'.format(universe_id)
            raise ValueError(msg)

        else:
            self._id = universe_id


    @name.setter
    def name(self, name):

        if not is_string(name):
            msg = 'Unable to set name for Universe ID={0} with a non-string ' \
                  'value {1}'.format(self._id, name)
            raise ValueError(msg)

        else:
            self._name = name


    def add_cell(self, cell):

        if not isinstance(cell, Cell):
            msg = 'Unable to add a Cell to Universe ID={0} since {1} is not ' \
                  'a Cell'.format(self._id, cell)
            raise ValueError(msg)

        cell_id = cell._id

        if not cell_id in self._cells:
            self._cells[cell_id] = cell


    def add_cells(self, cells):

        if not isinstance(cells, (list, tuple, np.ndarray)):
            msg = 'Unable to add Cells to Universe ID={0} since {1} is not a ' \
                  'Python tuple/list or NumPy array'.format(self._id, cells)
            raise ValueError(msg)

        for i in range(len(cells)):
            self.add_cell(cells[i])


    def remove_cell(self, cell):

        if not isinstance(cell, Cell):
            msg = 'Unable to remove a Cell from Universe ID={0} since {1} is ' \
                  'not a Cell'.format(self._id, cell)
            raise ValueError(msg)

        cell_id = cell.getId()

        # If the Cell is in the Universe's list of Cells, delete it
        if cell_id in self._cells:
            del self._cells[cell_id]


    def clear_cells(self):
        self._cells.clear()


    def get_offset(self, path, filter_offset):

        # Get the current element and remove it from the list
        path = path[1:]

        # Get the Cell ID
        cell_id = path[0]

        # Make a recursive call to the Cell within this Universe
        offset = self._cells[cell_id].get_offset(path, filter_offset)

        # Return the offset computed at all nested Universe levels
        return offset


    def get_all_nuclides(self):

        nuclides = {}

        # Append all Nuclides in each Cell in the Universe to the dictionary
        for cell_id, cell in self._cells.items():
            nuclides.update(cell.get_all_nuclides())

        return nuclides


    def get_all_cells(self):

        cells = {}

        # Add this Universe's cells to the dictionary
        cells.update(self._cells)

        # Append all Cells in each Cell in the Universe to the dictionary
        for cell_id, cell in self._cells.items():
            cells.update(cell.get_all_cells())

        return cells


    def get_all_universes(self):

        # Get all Cells in this Universe
        cells = self.get_all_cells()

        universes = {}

        # Append all Universes containing each Cell to the dictionary
        for cell_id, cell in cells.items():
            universes.update(cell.get_all_universes())

        return universes


    def __repr__(self):

        string = 'Universe\n'
        string += '{0: <16}{1}{2}\n'.format('\tID', '=\t', self._id)
        string += '{0: <16}{1}{2}\n'.format('\tName', '=\t', self._name)
        string += '{0: <16}{1}{2}\n'.format('\tCells', '=\t',
                                            list(self._cells.keys()))
        string += '{0: <16}{1}{2}\n'.format('\t# Regions', '=\t',
                                            self._num_regions)
        return string


    def create_xml_subelement(self, xml_element):

        # Iterate over all Cells
        for cell_id, cell in self._cells.items():

            # If the cell was not already written, write it
            if not cell_id in WRITTEN_IDS:
                WRITTEN_IDS[cell_id] = None

                # Create XML subelement for this Cell
                cell_subelement = cell.create_xml_subelement(xml_element)

                # Append the Universe ID to the subelement and add to Element
                cell_subelement.set("universe", str(self._id))
                xml_element.append(cell_subelement)



################################################################################
###################################  Lattice  ##################################
################################################################################


class Lattice(object):

    # This is an abstract class which cannot be instantiated
    metaclass__ = abc.ABCMeta

    def __init__(self, lattice_id=None, name=''):

        # Initialize Lattice class attributes
        self.id = lattice_id
        self.name = name
        self._pitch = None
        self._outer = None
        self._universes = None


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

        # Check that the ID is an integer and wasn't already used
        elif not is_integer(lattice_id):
            msg = 'Unable to set non-integer Lattice ID {0}'.format(lattice_id)
            raise ValueError(msg)

        elif lattice_id < 0:
            msg = 'Unable to set Lattice ID to {0} since it must be a ' \
                  'non-negative integer'.format(lattice_id)
            raise ValueError(msg)

        else:
            self._id = lattice_id


    @name.setter
    def name(self, name):

        if not is_string(name):
            msg = 'Unable to set name for Lattice ID={0} with a non-string ' \
                  'value {1}'.format(self._id, name)
            raise ValueError(msg)

        else:
            self._name = name


    @outer.setter
    def outer(self, outer):

        if not isinstance(outer, Universe):
            msg = 'Unable to set Lattice ID={0} outer universe to {1} ' \
                  'since it is not a Universe object'.format(self._id, outer)
            raise ValueError(msg)

        self._outer = outer


    @universes.setter
    def universes(self, universes):

        if not isinstance(universes, (tuple, list, np.ndarray)):
            msg = 'Unable to set Lattice ID={0} universes to {1} since ' \
                  'it is not a Python tuple/list or NumPy ' \
                  'array'.format(self._id, universes)
            raise ValueError(msg)

        self._universes = np.asarray(universes, dtype=Universe)


    def get_unique_universes(self):

        unique_universes = np.unique(self._universes.ravel())
        universes = {}

        for universe in unique_universes:
            universes[universe._id] = universe

        return universes


    def get_all_nuclides(self):

        nuclides = {}

        # Get all unique Universes contained in each of the lattice cells
        unique_universes = self.get_unique_universes()

        # Append all Universes containing each cell to the dictionary
        for universe_id, universe in unique_universes.items():
            nuclides.update(universe.get_all_nuclides())

        return nuclides


    def get_all_cells(self):

        cells = {}
        unique_universes = self.get_unique_universes()

        for universe_id, universe in unique_universes.items():
            cells.update(universe.get_all_cells())

        return cells


    def get_all_universes(self):

        # Initialize a dictionary of all Universes contained by the Lattice
        # in each nested Universe level
        all_universes = {}

        # Get all unique Universes contained in each of the lattice cells
        unique_universes = self.get_unique_universes()

        # Add the unique Universes filling each Lattice cell
        all_universes.update(unique_universes)

        # Append all Universes containing each cell to the dictionary
        for universe_id, universe in unique_universes.items():
            all_universes.update(universe.get_all_universes())

        return all_universes




################################################################################
#################################  RectLattice  ################################
################################################################################


class RectLattice(Lattice):

    def __init__(self, lattice_id=None, name=''):

        super(RectLattice, self).__init__(lattice_id, name)

        # Initialize Lattice class attributes
        self._dimension = None
        self._lower_left = None
        self._offsets = None


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

        if not isinstance(dimension, (tuple, list, np.ndarray)):
            msg = 'Unable to set RectLattice ID={0} dimension to {1} since ' \
                  'it is not a Python tuple/list or NumPy ' \
                  'array'.format(self._id, dimension)
            raise ValueError(msg)

        elif len(dimension) != 2 and len(dimension) != 3:
            msg = 'Unable to set RectLattice ID={0} dimension to {1} since ' \
                  'it does not contain 2 or 3 ' \
                  'coordinates'.format(self._id, dimension)
            raise ValueError(msg)

        for dim in dimension:

            if not is_integer(dim) and not is_float(dim):
                msg = 'Unable to set RectLattice ID={0} dimension to {1} since ' \
                      'it is not an integer or floating point ' \
                      'value'.format(self._id, dim)
                raise ValueError(msg)


            elif dim < 0:
                msg = 'Unable to set RectLattice ID={0} dimension to {1} ' \
                      'since it is a negative value'.format(self._id, dim)
                raise ValueError(msg)

        self._dimension = dimension


    @lower_left.setter
    def lower_left(self, lower_left):

        if not isinstance(lower_left, (tuple, list, np.ndarray)):
            msg = 'Unable to set RectLattice ID={0} lower_left to {1} since ' \
                  'it is not a Python tuple/list or NumPy ' \
                  'array'.format(self._id, lower_left)
            raise ValueError(msg)

        elif len(lower_left) != 2 and len(lower_left) != 3:
            msg = 'Unable to set RectLattice ID={0} lower_left to {1} ' \
                  'since it does not contain 2 or 3 ' \
                  'coordinates'.format(self._id, lower_left)
            raise ValueError(msg)


        for dim in lower_left:

            if not is_integer(dim) and not is_float(dim):
                msg = 'Unable to set RectLattice ID={0} lower_left to {1} since ' \
                      'it is is not an integer or floating point ' \
                      'value'.format(self._id, dim)
                raise ValueError(msg)

        self._lower_left = lower_left


    @offsets.setter
    def offsets(self, offsets):

        if not isinstance(offsets, (tuple, list, np.ndarray)):
            msg = 'Unable to set Lattice ID={0} offsets to {1} since ' \
                  'it is not a Python tuple/list or NumPy ' \
                  'array'.format(self._id, offsets)
            raise ValueError(msg)

        self._offsets = offsets


    @Lattice.pitch.setter
    def pitch(self, pitch):

        if not isinstance(pitch, (tuple, list, np.ndarray)):
            msg = 'Unable to set Lattice ID={0} pitch to {1} since ' \
                  'it is not a Python tuple/list or NumPy ' \
                  'array'.format(self._id, pitch)
            raise ValueError(msg)


        elif len(pitch) != 2 and len(pitch) != 3:
            msg = 'Unable to set Lattice ID={0} pitch to {1} since it does ' \
                  'not contain 2 or 3 coordinates'.format(self._id, pitch)
            raise ValueError(msg)

        for dim in pitch:

            if not is_integer(dim) and not is_float(dim):
                msg = 'Unable to set Lattice ID={0} pitch to {1} since ' \
                      'it is not an an integer or floating point ' \
                      'value'.format(self._id, dim)
                raise ValueError(msg)

            elif dim < 0:
                msg = 'Unable to set Lattice ID={0} pitch to {1} since it ' \
                      'is a negative value'.format(self._id, dim)
                raise ValueError(msg)

        self._pitch = pitch


    def get_offset(self, path, filter_offset):

        # Get the current element and remove it from the list
        i = path[0]
        path = path[1:]

        # For 2D Lattices
        if len(self._dimension) == 2:
            offset = self._offsets[i[1]-1, i[2]-1, 0, filter_offset-1]
            offset += self._universes[i[1]][i[2]].get_offset(path, filter_offset)

        # For 3D Lattices
        else:
            offset = self._offsets[i[1]-1, i[2]-1, i[3]-1, filter_offset-1]
            offset += self._universes[i[1]-1][i[2]-1][i[3]-1].get_offset(path,
                                                                 filter_offset)

        return offset


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


    def create_xml_subelement(self, xml_element):

        # Determine if XML element already contains subelement for this Lattice
        path = './lattice[@id=\'{0}\']'.format(self._id)
        test = xml_element.find(path)

        # If the element does contain the Lattice subelement, then return
        if not test is None:
            return

        lattice_subelement = ET.Element("lattice")
        lattice_subelement.set("id", str(self._id))

        if len(self._name) > 0:
            lattice_subelement.set("name", str(self._name))

        # Export the Lattice cell pitch
        if len(self._pitch) == 3:
            pitch = ET.SubElement(lattice_subelement, "pitch")
            pitch.text = '{0} {1} {2}'.format(self._pitch[0], \
                                              self._pitch[1], \
                                              self._pitch[2])
        else:
            pitch = ET.SubElement(lattice_subelement, "pitch")
            pitch.text = '{0} {1}'.format(self._pitch[0], \
                                          self._pitch[1])

        # Export the Lattice outer Universe (if specified)
        if self._outer is not None:
            outer = ET.SubElement(lattice_subelement, "outer")
            outer.text = '{0}'.format(self._outer._id)

        # Export Lattice cell dimensions
        if len(self._dimension) == 3:
            dimension = ET.SubElement(lattice_subelement, "dimension")
            dimension.text = '{0} {1} {2}'.format(self._dimension[0], \
                                                  self._dimension[1], \
                                                  self._dimension[2])
        else:
            dimension = ET.SubElement(lattice_subelement, "dimension")
            dimension.text = '{0} {1}'.format(self._dimension[0], \
                                              self._dimension[1])

        # Export Lattice lower left
        if len(self._lower_left) == 3:
            lower_left = ET.SubElement(lattice_subelement, "lower_left")
            lower_left.text = '{0} {1} {2}'.format(self._lower_left[0], \
                                                   self._lower_left[1], \
                                                   self._lower_left[2])
        else:
            lower_left = ET.SubElement(lattice_subelement, "lower_left")
            lower_left.text = '{0} {1}'.format(self._lower_left[0], \
                                               self._lower_left[1])

        # Export the Lattice nested Universe IDs - column major for Fortran
        universe_ids = '\n'

        # 3D Lattices
        if len(self._dimension) == 3:
            for z in range(self._dimension[2]):
                for y in range(self._dimension[1]):
                    for x in range(self._dimension[0]):

                        universe = self._universes[x][y][z]

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

                    universe = self._universes[x][y]

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


################################################################################
##################################  HexLattice  ################################
################################################################################


class HexLattice(Lattice):

    def __init__(self, lattice_id=None, name=''):

        super(HexLattice, self).__init__(lattice_id, name)

        # Initialize Lattice class attributes
        self._num_rings = None
        self._num_axial = None
        self._center = None


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

        if not is_integer(num_rings) and num_rings < 1:
            msg = 'Unable to set HexLattice ID={0} number of rings to {1} ' \
                  'since it is not a positive integer'.format(self._id, num_rings)
            raise ValueError(msg)

        self._num_rings = num_rings


    @num_axial.setter
    def num_axial(self, num_axial):

        if not is_integer(num_axial) and num_axial < 1:
            msg = 'Unable to set HexLattice ID={0} number of axial to {1} ' \
                  'since it is not a positive integer'.format(self._id, num_axial)
            raise ValueError(msg)

        self._num_axial = num_axial


    @center.setter
    def center(self, center):

        if not isinstance(center, (tuple, list, np.ndarray)):
            msg = 'Unable to set HexLattice ID={0} dimension to {1} since ' \
                  'it is not a Python tuple/list or NumPy ' \
                  'array'.format(self._id, center)
            raise ValueError(msg)

        elif len(center) != 2 and len(center) != 3:
            msg = 'Unable to set HexLattice ID={0} center to {1} since ' \
                  'it does not contain 2 or 3 ' \
                  'coordinates'.format(self._id, center)
            raise ValueError(msg)

        for dim in center:

            if not is_integer(dim) and not is_float(dim):
                msg = 'Unable to set HexLattice ID={0} center to {1} since ' \
                      'it is not an integer or floating point ' \
                      'value'.format(self._id, dim)
                raise ValueError(msg)

        self._center = center


    @Lattice.pitch.setter
    def pitch(self, pitch):

        if not isinstance(pitch, (tuple, list, np.ndarray)):
            msg = 'Unable to set Lattice ID={0} pitch to {1} since ' \
                  'it is not a Python tuple/list or NumPy ' \
                  'array'.format(self._id, pitch)
            raise ValueError(msg)


        elif len(pitch) != 1 and len(pitch) != 2:
            msg = 'Unable to set Lattice ID={0} pitch to {1} since it does ' \
                  'not contain 2 or 3 coordinates'.format(self._id, pitch)
            raise ValueError(msg)

        for dim in pitch:

            if not is_integer(dim) and not is_float(dim):
                msg = 'Unable to set Lattice ID={0} pitch to {1} since ' \
                      'it is not an an integer or floating point ' \
                      'value'.format(self._id, dim)
                raise ValueError(msg)

            elif dim < 0:
                msg = 'Unable to set Lattice ID={0} pitch to {1} since it ' \
                      'is a negative value'.format(self._id, dim)
                raise ValueError(msg)

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
            self.num_axial = self._universes.shape[0]

        else:
            self._num_axial = None

        # Set the number of rings and make sure this number is consistent for
        # all axial positions.
        if n_dims == 3:
            self.num_rings = len(self._universes[0])
            for rings in self._universes:
                if len(rings) != self._num_rings:
                    msg = 'HexLattice ID={0:d} has an inconsistent number of ' \
                          'rings per axial positon'.format(self._id)
                    raise ValueError(msg)

        else:
            self.num_rings = self._universes.shape[0]

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


    def create_xml_subelement(self, xml_element):

        # Determine if XML element already contains subelement for this Lattice
        path = './hex_lattice[@id=\'{0}\']'.format(self._id)
        test = xml_element.find(path)

        # If the element does contain the Lattice subelement, then return
        if not test is None:
            return

        lattice_subelement = ET.Element("hex_lattice")
        lattice_subelement.set("id", str(self._id))

        if len(self._name) > 0:
            lattice_subelement.set("name", str(self._name))

        # Export the Lattice cell pitch
        if len(self._pitch) == 2:
            pitch = ET.SubElement(lattice_subelement, "pitch")
            pitch.text = '{0} {1}'.format(self._pitch[0], \
                                          self._pitch[1])
        else:
            pitch = ET.SubElement(lattice_subelement, "pitch")
            pitch.text = '{0}'.format(self._pitch[0])

        # Export the Lattice outer Universe (if specified)
        if self._outer is not None:
            outer = ET.SubElement(lattice_subelement, "outer")
            outer.text = '{0}'.format(self._outer._id)

        lattice_subelement.set("n_rings", str(self._num_rings))

        if self._num_axial is not None:
            lattice_subelement.set("n_axial", str(self._num_axial))

        # Export Lattice cell center
        if len(self._center) == 3:
            dimension = ET.SubElement(lattice_subelement, "center")
            dimension.text = '{0} {1} {2}'.format(self._center[0], \
                                                  self._center[1], \
                                                  self._center[2])
        else:
            dimension = ET.SubElement(lattice_subelement, "center")
            dimension.text = '{0} {1}'.format(self._center[0], \
                                              self._center[1])

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
            for r in range(self._num_rings-1):
                for theta in range(2*(self._num_rings - r)):
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
            if y != 0: rows[middle - y] = pad + rows[middle - y]

        # Join the rows together and return the string.
        universe_ids = '\n'.join(rows)
        return universe_ids
