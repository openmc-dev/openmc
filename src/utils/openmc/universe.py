#!/usr/bin/env python

import openmc
from openmc.checkvalue import *
from xml.etree import ElementTree as ET
from collections import OrderedDict
import operator
import numpy as np


################################################################################
####################################  Cell  ####################################
################################################################################

# A static variable for auto-generated Cell IDs
AUTO_CELL_ID = 10000

def reset_auto_cell_id():
  global AUTO_CELL_ID
  AUTO_CELL_ID = 10000


class Cell(object):

  def __init__(self, cell_id=None, name=''):

    # Initialize Cell class attributes
    self._id = None
    self._name = None
    self._fill = None
    self._type = None
    self._surfaces = dict()
    self._rotation = None
    self._translation = None
    self._offset = None

    self.setId(cell_id)
    self.setName(name)


  def getOffset(self, path, filter_offset):

    # Get the current element and remove it from the list
    cell_id = path[0]
    path = path[1:]

    # If the Cell is filled by a Material
    if self._type == 'normal':
      offset = 0

    # If the Cell is filled by a Universe
    elif self._type == 'fill':
      offset = self._offset[filter_offset-1]
      offset += self._fill.getOffset(path, filter_offset)

    # If the Cell is filled by a Lattice
    else:
      offset = self._fill.getOffset(path, filter_offset)

    return offset

    # Make a recursive call to the Universe filling this Cell
    offset = self._cells[cell_id].getOffset(path, filter_offset)

    # Return the offset computed at all nested Universe levels below this one
    return offset


  def setId(self, cell_id=None):

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


  def setName(self, name):

    if not isinstance(name, str):
      msg = 'Unable to set name for Cell ID={0} with a non-string ' \
            'value {1}'.format(self._id, name)
      raise ValueError(msg)

    else:
      self._name = name


  def setFill(self, fill):

    if not isinstance(fill, (openmc.Material, Universe, Lattice)) \
      and fill != 'void':
      msg = 'Unable to set Cell ID={0} to use a a non-Material or ' \
            'Universe fill {1}'.format(self._id, fill)
      raise ValueError(msg)

    self._fill = fill

    if isinstance(fill, Lattice):
      self._type = 'lattice'
    elif isinstance(fill, Universe):
      self._type = 'fill'
    elif fill == 'void':
      self._type = 'normal'
    else:
      self._type = 'normal'


  def addSurface(self, surface, halfspace):

    if not isinstance(surface, openmc.Surface):
      msg = 'Unable to add Surface {0} to Cell ID={1} since it is ' \
            'not a Surface object'.format(surface, self._id)
      raise ValueError(msg)

    if not halfspace in [-1, +1]:
      msg = 'Unable to add Surface {0} to Cell ID={1} with halfspace {2} ' \
            'since it is not +/-1'.format(surface, self._id, halfspace)
      raise ValueError(msg)

    # If the Cell does not already contain the Surface, add it
    if not surface._id in self._surfaces:
      self._surfaces[surface._id] = (surface, halfspace)


  def removeSurface(self, surface):

    if not isinstance(surface, openmc.Surface):
      msg = 'Unable to remove Surface {0} from Cell ID={1} since it is ' \
            'not a Surface object'.format(surface, self._id)
      raise ValueError(msg)

    # If the Cell contains the Surface, delete it
    if surface._id in self._surfaces:
      del self._surfaces[surface._id]


  def setSurfaces(self, surfaces, halfspaces):

    if not isinstance(surfaces, (tuple, list, np.ndarray)):
      msg = 'Unable to set Cell ID={0} with Surfaces {1} since it is not a ' \
            'Python tuple/list or NumPy array'.format(self._id, surfaces)
      raise ValueError(msg)

    if not isinstance(halfspaces, (tuple, list, np.ndarray)):
      msg = 'Unable to set Cell ID={0} with Surface halfspaces {1} ' \
            'since it is not a Python tuple/list or NumPy ' \
            'array'.format(self._id, halfspaces)
      raise ValueError(msg)

    for surface in surfaces:
      self.addSurface(surface)


  def setRotation(self, rotation):

    if not isinstance(rotation, (tuple, list, np.ndarray)):
      msg = 'Unable to add rotation {0} to Cell ID={1} since it is not a ' \
            'Python tuple/list or NumPy array'.format(rotation, self._id)
      raise ValueError(msg)

    elif len(rotation) != 3:
      msg = 'Unable to add rotation {0} to Cell ID={1} since it does not ' \
            'contain 3 values'.format(rotation, self._id)
      raise ValueError(msg)

    for axis in rotation:

      if not is_integer(axis) and not is_float(axis):
        msg = 'Unable to add rotation {0} to Cell ID={1} since it is not ' \
              'an integer or floating point value'.format(axis, self._id)
        raise ValueError(msg)

      self._rotation = rotation


  def setTranslation(self, translation):

    if not isinstance(translation, (tuple, list, np.ndarray)):
      msg = 'Unable to add translation {0} to Cell ID={1} since it is not a ' \
            'Python tuple/list or NumPy array'.format(translation, self._id)
      raise ValueError(msg)

    elif len(translation) != 3:
      msg = 'Unable to add translation {0} to Cell ID={1} since it does not ' \
            'contain 3 values'.format(translation, self._id)
      raise ValueError(msg)

    for axis in translation:
      if not is_integer(axis) and not is_float(axis):
        msg = 'Unable to add translation {0} to Cell ID={1} since it is ' \
              'not an integer or floating point value'.format(axis, self._id)
        raise ValueError(msg)

    self._translation = translation

  def setOffset(self, offset):

    if not isinstance(offset, (tuple, list, np.ndarray)):
      msg = 'Unable to set offset {0} to Cell ID={1} since it is not a ' \
            'Python tuple/list or NumPy array'.format(offset, self._id)
      raise ValueError(msg)

    self._offset = offset


  def getAllNuclides(self):

    nuclides = dict()

    if self._type != 'void':
      nuclides.update(self._fill.getAllNuclides())

    return nuclides


  def getAllCells(self):

    cells = dict()

    if self._type == 'fill' or self._type == 'lattice':
      cells.update(self._fill.getAllCells())

    return cells


  def getAllUniverses(self):

    universes = dict()

    if self._type == 'fill':
      universes[self._fill._id] = self._fill
      universes.update(self._fill.getAllUniverses())
    elif self._type == 'lattice':
      universes.update(self._fill.getAllUniverses())

    return universes


  def __repr__(self):

    string = 'Cell\n'
    string += '{0: <16}{1}{2}\n'.format('\tID', '=\t', self._id)
    string += '{0: <16}{1}{2}\n'.format('\tName', '=\t', self._name)

    if isinstance(self._fill, openmc.Material):
      string += '{0: <16}{1}{2}\n'.format('\tMaterial', '=\t', self._fill._id)
    elif isinstance(self._fill, (Universe, Lattice)):
      string += '{0: <16}{1}{2}\n'.format('\tFill', '=\t', self._fill._id)
    else:
      string += '{0: <16}{1}{2}\n'.format('\tFill', '=\t', self._fill)

    string += '{0: <16}{1}\n'.format('\tSurfaces', '=\t')

    for surface_id in self._surfaces:
      halfspace = self._surfaces[surface_id][1]
      string += '{0} '.format(halfspace * surface_id)

    string = string.rstrip(' ') + '\n'

    string += '{0: <16}{1}{2}\n'.format('\tRotation', '=\t', self._rotation)
    string += '{0: <16}{1}{2}\n'.format('\tTranslation', '=\t', self._translation)
    string += '{0: <16}{1}{2}\n'.format('\tOffset', '=\t', self._offset)

    return string


  def createXMLSubElement(self, xml_element):

    element = ET.Element("cell")
    element.set("id", str(self._id))

    if isinstance(self._fill, openmc.Material):
      element.set("material", str(self._fill._id))

    elif isinstance(self._fill, (Universe, Lattice)):
      element.set("fill", str(self._fill._id))
      self._fill.createXMLSubElement(xml_element)

    else:
      element.set("fill", str(self._fill))
      self._fill.createXMLSubElement(xml_element)


    if not self._surfaces is None:

      surfaces = ''

      for surface_id in self._surfaces:

        # Determine if XML element already includes subelement for this Surface
        path = './surface[@id=\'{0}\']'.format(surface_id)
        test = xml_element.find(path)

        # If the element does not contain the Surface subelement, then add it
        if test is None:

          # Create the XML subelement for this Surface and add it to the element
          surface = self._surfaces[surface_id][0]
          surface_subelement = surface.createXMLSubElement()

          if len(surface._name) > 0:
            xml_element.append(ET.Comment(surface._name))

          xml_element.append(surface_subelement)

        # Append to the Surfaces's XML attribute the halfspace and Surface ID
        halfspace = self._surfaces[surface_id][1]
        surfaces += '{0} '.format(halfspace * surface_id)

      element.set("surfaces", surfaces.rstrip(' '))


    if not self._translation is None:

      translation  = ''

      for axis in self._translation:
        translation += '{0} '.format(axis)

      element.set("translation", translation.rstrip(' '))


    if not self._rotation is None:

      rotation  = ''

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
    self._id = None
    self._name = None

    # Keys   - Cell IDs
    # Values - Cells
    self._cells = dict()

    # Keys   - Cell IDs
    # Values - Offsets
    self._cell_offsets = OrderedDict()
    self._num_regions = 0

    self.setId(universe_id)
    self.setName(name)


  def setId(self, universe_id=None):

    if universe_id is None:
      global AUTO_UNIVERSE_ID
      self._id = AUTO_UNIVERSE_ID
      AUTO_UNIVERSE_ID += 1

    # Check that the ID is an integer and wasn't already used
    elif not is_integer(universe_id):
      msg = 'Unable to set Universe ID to a non-integer {0}'.format(universe_id)
      raise ValueError(msg)

    elif universe_id < 0:
      msg = 'Unable to set Universe ID to {0} since it must be a ' \
            'non-negative integer'.format(universe_id)
      raise ValueError(msg)

    else:
      self._id = universe_id


  def setName(self, name):

    if not is_string(name):
      msg = 'Unable to set name for Universe ID={0} with a non-string ' \
            'value {1}'.format(self._id, name)
      raise ValueError(msg)

    else:
      self._name = name


  def addCell(self, cell):

    if not isinstance(cell, Cell):
      msg = 'Unable to add a Cell to Universe ID={0} since {1} is not ' \
            'a Cell'.format(self._id, cell)
      raise ValueError(msg)

    cell_id = cell._id

    if not cell_id in self._cells:
      self._cells[cell_id] = cell


  def addCells(self, cells):

    if not isinstance(cells, (list, tuple, np.ndarray)):
      msg = 'Unable to add Cells to Universe ID={0} since {1} is not a ' \
            'Python tuple/list or NumPy array'.format(self._id, cells)
      raise ValueError(msg)

    for i in range(len(cells)):
      self.addCell(cells[i])


  def removeCell(self, cell):

    if not isinstance(cell, Cell):
      msg = 'Unable to remove a Cell from Universe ID={0} since {1} is ' \
            'not a Cell'.format(self._id, cell)
      raise ValueError(msg)

    cell_id = cell.getId()

    # If the Cell is in the Universe's list of Cells, delete it
    if cell_id in self._cells:
      del self._cells[cell_id]


  def clearCells(self):
    self._cells.clear()


  def getOffset(self, path, filter_offset):

    # Get the current element and remove it from the list
    path = path[1:]

    # Get the Cell ID
    cell_id = path[0]

    # Make a recursive call to the Cell within this Universe
    offset = self._cells[cell_id].getOffset(path, filter_offset)

    # Return the offset computed at all nested Universe levels below this one
    return offset


  def getAllNuclides(self):

    nuclides = dict()

    # Append all Nuclides in each Cell in the Universe to the dictionary
    for cell_id, cell in self._cells.items():
      nuclides.update(cell.getAllNuclides())

    return nuclides


  def getAllCells(self):

    cells = dict()

    # Add this Universe's cells to the dictionary
    cells.update(self._cells)

    # Append all Cells in each Cell in the Universe to the dictionary
    for cell_id, cell in self._cells.items():
      cells.update(cell.getAllCells())

    return cells


  def getAllUniverses(self):

    # Get all Cells in this Universe
    cells = self.getAllCells()

    universes = dict()

    # Append all Universes containing each Cell to the dictionary
    for cell_id, cell in cells.items():
      universes.update(cell.getAllUniverses())

    return universes


  def __repr__(self):

    string = 'Universe\n'
    string += '{0: <16}{1}{2}\n'.format('\tID', '=\t', self._id)
    string += '{0: <16}{1}{2}\n'.format('\tName', '=\t', self._name)
    string += '{0: <16}{1}{2}\n'.format('\tCells', '=\t',
                                        list(self._cells.keys()))
    string += '{0: <16}{1}{2}\n'.format('\t# Regions', '=\t', self._num_regions)
    return string


  def createXMLSubElement(self, xml_element):

    # Iterate over all Cells
    for cell_id, cell in self._cells.items():

      # Determine if XML element already contains subelement for this Cell
      path = './cell[@id=\'{0}\']'.format(cell_id)
      test = xml_element.find(path)

      # If the element does not contain the Cell subelement, then add it
      if test is None:

        # Create the XML subelement for this Cell and everything beneath it
        cell_subelement = cell.createXMLSubElement(xml_element)

        # Append the Universe ID to the subelement and add it to the element
        cell_subelement.set("universe", str(self._id))

        if len(cell._name) > 0:
          xml_element.append(ET.Comment(cell._name))

        xml_element.append(cell_subelement)



################################################################################
###################################  Lattice  ##################################
################################################################################


class Lattice(object):

  def __init__(self, lattice_id=None, name='', type='rectangular'):

    # Initialize Lattice class attributes
    self._id = None
    self._name = None
    self._type = ''
    self._dimension = None
    self._lower_left = None
    self._width = None
    self._outside = None
    self._universes = None
    self._offsets = None


    self.setId(lattice_id)
    self.setName(name)
    self.setType(type)


  def setId(self, lattice_id=None):

    if lattice_id is None:
      global AUTO_UNIVERSE_ID
      self._id = AUTO_UNIVERSE_ID
      AUTO_UNIVERSE_ID += 1

    # Check that the ID is an integer and wasn't already used
    elif not is_integer(lattice_id):
      msg = 'Unable to set a non-integer Lattice ID {0}'.format(lattice_id)
      raise ValueError(msg)

    elif lattice_id < 0:
      msg = 'Unable to set Lattice ID to {0} since it must be a ' \
           'non-negative integer'.format(lattice_id)
      raise ValueError(msg)

    else:
      self._id = lattice_id


  def setName(self, name):

    if not is_string(name):
      msg = 'Unable to set name for Lattice ID={0} with a non-string ' \
            'value {1}'.format(self._id, name)
      raise ValueError(msg)

    else:
      self._name = name


  def setType(self, type):

    if not is_string(type):
      msg = 'Unable to set type for Lattice ID={0} with a non-string ' \
           'value {1}'.format(self._id, type)
      raise ValueError(msg)

    elif not type in ['rectangular', 'hexagonal']:
      msg = 'Unable to set type for Lattice ID={0} as {1} since it is not ' \
           'rectangular or hexagonal'.format(self._id, type)
      raise ValueError(msg)

    self._type = type


  def setDimension(self, dimension):

    if not isinstance(dimension, (tuple, list, np.ndarray)):
      msg = 'Unable to set Lattice ID={0} dimension to {1} since it is not ' \
            'a Python tuple/list or NumPy array'.format(self._id, dimension)
      raise ValueError(msg)

    elif len(dimension) != 2 and len(dimension) != 3:
      msg = 'Unable to set Lattice ID={0} dimension to {1} since it does ' \
            'not contain 2 or 3 coordinates'.format(self._id, dimension)
      raise ValueError(msg)

    for dim in dimension:

      if not is_integer(dim) and not is_float(dim):
        msg = 'Unable to set Lattice ID={0} dimension to {1} since it is ' \
             'not an integer or floating point value'.format(self._id, dim)
        raise ValueError(msg)


      elif dim < 0:
        msg = 'Unable to set Lattice ID={0} dimension to {1} since it ' \
              'is a negative value'.format(self._id, dim)
        raise ValueError(msg)

    self._dimension = dimension


  def setLowerLeft(self, lower_left):

    if not isinstance(lower_left, (tuple, list, np.ndarray)):
      msg = 'Unable to set Lattice ID={0} lower_left to {1} since it is not ' \
            'a Python tuple/list or NumPy array'.format(self._id, lower_left)
      raise ValueError(msg)

    elif len(lower_left) != 2 and len(lower_left) != 3:
      msg = 'Unable to set Lattice ID={0} lower_left to {1} since it does ' \
            'not contain 2 or 3 coordinates'.format(self._id, lower_left)
      raise ValueError(msg)


    for dim in lower_left:

      if not is_integer(dim) and not is_float(dim):
        msg = 'Unable to set Lattice ID={0} lower_left to {1} since it is ' \
              'is not an integer or floating point value'.format(self._id, dim)
        raise ValueError(msg)

    self._lower_left = lower_left


  def setWidth(self, width):

    if not isinstance(width, (tuple, list, np.ndarray)):
      msg = 'Unable to set Lattice ID={0} width to {1} since it is not a ' \
             'Python tuple/list or NumPy array'.format(self._id, width)
      raise ValueError(msg)


    elif len(width) != 2 and len(width) != 3:
      msg = 'Unable to set Lattice ID={0} width to {1} since it does ' \
            'not contain 2 or 3 coordinates'.format(self._id, width)
      raise ValueError(msg)

    for dim in width:

      if not is_integer(dim) and not is_float(dim):
        msg = 'Unable to set Lattice ID={0} width to {1} since it is not an ' \
              'an integer or floating point value'.format(self._id, dim)
        raise ValueError(msg)

      elif dim < 0:
        msg = 'Unable to set Lattice ID={0} width to {1} since it ' \
              'is a negative value'.format(self._id, dim)
        raise ValueError(msg)

    self._width = width


  def setOutside(self, outside):

    if not isinstance(outside, Universe):
      msg = 'Unable to set Lattice ID={0} outside universe to {1} ' \
            'since it is not a Universe object'.format(self._id, outside)
      raise ValueError(msg)

    self._outside = outside


  def setUniverses(self, universes):

    if not isinstance(universes, (tuple, list, np.ndarray)):
      msg = 'Unable to set Lattice ID={0} universes to {1} since it is not ' \
            'a Python tuple/list or NumPy array'.format(self._id, universes)
      raise ValueError(msg)

    self._universes = np.asarray(universes, dtype=Universe)


  def setOffsets(self, offsets):

    if not isinstance(offsets, (tuple, list, np.ndarray)):
      msg = 'Unable to set Lattice ID={0} offsets to {1} since it is not ' \
            'a Python tuple/list or NumPy array'.format(self._id, offsets)
      raise ValueError(msg)

    self._offsets = offsets


  def getOffset(self, path, filter_offset):

    # Get the current element and remove it from the list
    i = path[0]
    path = path[1:]

    # For 2D Lattices
    if len(self._dimension) == 2:
      offset = self._offsets[i[1]-1, i[2]-1, 0, filter_offset-1]
      offset += self._universes[i[1]][i[2]].getOffset(path, filter_offset)

    # For 3D Lattices
    else:
      offset = self._offsets[i[1]-1, i[2]-1, i[3]-1, filter_offset-1]
      offset += self._universes[i[1]-1][i[2]-1][i[3]-1].getOffset(path, \
                                                                  filter_offset)

    return offset


  def getUniqueUniverses(self):

    unique_universes = np.unique(self._universes.ravel())
    universes = dict()

    for universe in unique_universes:
      universes[universe._id] = universe

    return universes


  def getAllNuclides(self):

    nuclides = dict()

    # Get all unique Universes contained in each of the lattice cells
    unique_universes = self.getUniqueUniverses()

    # Append all Universes containing each cell to the dictionary
    for universe_id, universe in unique_universes.items():
      nuclides.update(universe.getAllNuclides())

    return nuclides


  def getAllCells(self):

    cells = dict()
    unique_universes = self.getUniqueUniverses()

    for universe_id, universe in unique_universes.items():
      cells.update(universe.getAllCells())

    return cells


  def getAllUniverses(self):

    # Initialize a dictionary of all Universes contained by the Lattice
    # in each nested Universe level
    all_universes = dict()

    # Get all unique Universes contained in each of the lattice cells
    unique_universes = self.getUniqueUniverses()

    # Add the unique Universes filling each Lattice cell
    all_universes.update(unique_universes)

    # Append all Universes containing each cell to the dictionary
    for universe_id, universe in unique_universes.items():
      all_universes.update(universe.getAllUniverses())

    return all_universes


  def __repr__(self):

    string = 'Lattice\n'
    string += '{0: <16}{1}{2}\n'.format('\tID', '=\t', self._id)
    string += '{0: <16}{1}{2}\n'.format('\tName', '=\t', self._name)
    string += '{0: <16}{1}{2}\n'.format('\tType', '=\t', self._type)
    string += '{0: <16}{1}{2}\n'.format('\tDimension', '=\t', self._dimension)
    string += '{0: <16}{1}{2}\n'.format('\tLower Left', '=\t', self._lower_left)
    string += '{0: <16}{1}{2}\n'.format('\tWidth', '=\t', self._width)

    if self._outside is not None:
      string += '{0: <16}{1}{2}\n'.format('\tOutside', '=\t', self._outside._id)
    else:
      string += '{0: <16}{1}{2}\n'.format('\tOutside', '=\t', self._outside)

    string += '{0: <16}\n'.format('\tUniverses')

    # Lattice nested Universe IDs - column major for Fortran
    for i, universe in enumerate(np.ravel(self._universes)):
      string += '{0} '.format(universe._id)

      # Add a newline character every time we reach the end of a row of cells
      if (i+1) % self._dimension[-1] == 0:
        string += '\n'

    string = string.rstrip('\n')

    if self._offsets is not None:
      string += '{0: <16}\n'.format('\tOffsets')

      # Lattice cell offsets
      for i, offset in enumerate(np.ravel(self._offsets)):

        string += '{0} '.format(offset)

        # Add a newline character every time we reach the end of a row of cells
        if (i+1) % self._dimension[-1] == 0:
          string += '\n'

      string = string.rstrip('\n')

    return string


  def createXMLSubElement(self, xml_element):

    # Determine if XML element already contains subelement for this Lattice
    path = './lattice[@id=\'{0}\']'.format(self._id)
    test = xml_element.find(path)

    # If the element does contain the Lattice subelement, then return
    if not test is None:
      return

    lattice_subelement = ET.Element("lattice")
    lattice_subelement.set("id", str(self._id))
    lattice_subelement.set("type", self._type)

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

    # Export the Lattice cell width/height
    if len(self._width) == 3:
      width = ET.SubElement(lattice_subelement, "width")
      width.text = '{0} {1} {2}'.format(self._width[0], \
                                        self._width[1], \
                                        self._width[2])
    else:
      width = ET.SubElement(lattice_subelement, "width")
      width.text = '{0} {1}'.format(self._width[0], \
                                    self._width[1])

    # Export the Lattice outside Universe (if specified)
    if self._outside is not None:
      outside = ET.SubElement(lattice_subelement, "outside")
      outside.text = '{0}'.format(self._outside._id)

    # Export the Lattice nested Universe IDs - column major for Fortran
    universe_ids = '\n'

    # 3D Lattices
    if len(self._dimension) == 3:
      for z in range(self._dimension[2]):
        for y in range(self._dimension[1]):
          for x in range(self._dimension[0]):

            universe = self._universes[x][y][z]

            # Append Universe ID to the string for the Lattice XML subelement
            universe_ids += '{0} '.format(universe._id)

            # Create XML subelement for this Universe and everything beneath it
            universe.createXMLSubElement(xml_element)

          # Add newline character every time we reach the end of row of cells
          universe_ids += '\n'

        # Add newline character every time we reach the end of row of cells
        universe_ids += '\n'

    # 2D Lattices
    else:
      for y in range(self._dimension[1]):
        for x in range(self._dimension[0]):

          universe = self._universes[x][y]

          # Append Universe ID to the string for the Lattice XML subelement
          universe_ids += '{0} '.format(universe._id)

          # Create XML subelement for this Universe and everything beneath it
          universe.createXMLSubElement(xml_element)

        # Add newline character every time we reach the end of row of cells
        universe_ids += '\n'

    # Remove trailing newline character from Universe IDs string
    universe_ids = universe_ids.rstrip('\n')

    universes = ET.SubElement(lattice_subelement, "universes")
    universes.text = universe_ids

    if len(self._name) > 0:
      xml_element.append(ET.Comment(self._name))

    # Append the XML subelement for this Lattice to the XML element
    xml_element.append(lattice_subelement)
