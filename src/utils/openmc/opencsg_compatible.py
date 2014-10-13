#!/usr/bin/env python

import openmc
import opencsg
import copy
import numpy as np


# A dictionary of all OpenMC Materials created
# Keys    - Material IDs
# Values  - Materials
OPENMC_MATERIALS = dict()

# A dictionary of all OpenCSG Materials created
# Keys    - Material IDs
# Values  - Materials
OPENCSG_MATERIALS = dict()

# A dictionary of all OpenMC Surfaces created
# Keys    - Surface IDs
# Values  - Surfaces
OPENMC_SURFACES = dict()

# A dictionary of all OpenCSG Surfaces created
# Keys    - Surface IDs
# Values  - Surfaces
OPENCSG_SURFACES = dict()

# A dictionary of all OpenMC Cells created
# Keys    - Cell IDs
# Values  - Cells
OPENMC_CELLS = dict()

# A dictionary of all OpenCSG Cells created
# Keys    - Cell IDs
# Values  - Cells
OPENCSG_CELLS = dict()

# A dictionary of all OpenMC Universes created
# Keys    - Universes IDs
# Values  - Universes
OPENMC_UNIVERSES = dict()

# A dictionary of all OpenCSG Universes created
# Keys    - Universes IDs
# Values  - Universes
OPENCSG_UNIVERSES = dict()

# A dictionary of all OpenMC Lattices created
# Keys    - Lattice IDs
# Values  - Lattices
OPENMC_LATTICES = dict()

# A dictionary of all OpenCSG Lattices created
# Keys    - Lattice IDs
# Values  - Lattices
OPENCSG_LATTICES = dict()



def get_opencsg_material(openmc_material):

  if not isinstance(openmc_material, openmc.Material):
    msg = 'Unable to create an OpenCSG Material from {0} ' \
          'which is not an OpenMC Material'.format(openmc_material)
    raise ValueError(msg)

  global OPENCSG_MATERIALS
  material_id = openmc_material._id

  # If this Material was already created, use it
  if material_id in OPENCSG_MATERIALS.keys():
    return OPENCSG_MATERIALS[material_id]

  # Create an OpenCSG Material to represent this OpenMC Material
  name = openmc_material._name
  opencsg_material = opencsg.Material(material_id=material_id, name=name)

  # Add the OpenMC Material to the global collection of all OpenMC Materials
  OPENMC_MATERIALS[material_id] = openmc_material

  # Add the OpenCSG Material to the global collection of all OpenCSG Materials
  OPENCSG_MATERIALS[material_id] = opencsg_material

  return opencsg_material


def get_openmc_material(opencsg_material):

  if not isinstance(opencsg_material, opencsg.Material):
    msg = 'Unable to create an OpenMC Material from {0} ' \
          'which is not an OpenCSG Material'.format(opencsg_material)
    raise ValueError(msg)

  global OPENMC_MATERIALS
  material_id = opencsg_material._id

  # If this Material was already created, use it
  if material_id in OPENMC_MATERIALS.keys():
    return OPENMC_MATERIALS[material_id]

  # Create an OpenMC Material to represent this OpenCSG Material
  name = opencsg_material._name
  openmc_material = openmc.Material(material_id=material_id, name=name)

  # Add the OpenMC Material to the global collection of all OpenMC Materials
  OPENMC_MATERIALS[material_id] = openmc_material

  # Add the OpenCSG Material to the global collection of all OpenCSG Materials
  OPENCSG_MATERIALS[material_id] = opencsg_material

  return openmc_material


def is_opencsg_surface_compatible(opencsg_surface):

  if not isinstance(opencsg_surface, opencsg.Surface):
    msg = 'Unable to check if OpenCSG Surface is compatible' \
          'since {0} is not a Surface'.format(opencsg_surface)
    raise ValueError(msg)

  if opencsg_surface._type in ['x-squareprism', 'y-squareprism', 'z-squareprism']:
    return False
  else:
    return True


def get_opencsg_surface(openmc_surface):

  if not isinstance(openmc_surface, openmc.Surface):
    msg = 'Unable to create an OpenCSG Surface from {0} ' \
          'which is not an OpenMC Surface'.format(openmc_surface)
    raise ValueError(msg)

  global OPENCSG_SURFACES
  surface_id = openmc_surface._id

  # If this Material was already created, use it
  if surface_id in OPENCSG_SURFACES.keys():
    return OPENCSG_SURFACES[surface_id]

  # Create an OpenCSG Surface to represent this OpenMC Surface
  name = openmc_surface._name

  # Correct for OpenMC's syntax for Surfaces dividing Cells
  boundary = openmc_surface._bc_type
  if boundary == 'transmission':
    boundary = 'interface'

  opencsg_surface = None

  if openmc_surface._type == 'plane':
    A = openmc_surface._coeffs['A']
    B = openmc_surface._coeffs['B']
    C = openmc_surface._coeffs['C']
    D = openmc_surface._coeffs['D']
    opencsg_surface = opencsg.Plane(surface_id, name, boundary, A, B, C, D)

  elif openmc_surface._type == 'x-plane':
    x0 = openmc_surface._coeffs['x0']
    opencsg_surface = opencsg.XPlane(surface_id, name, boundary, x0)

  elif openmc_surface._type == 'y-plane':
    y0 = openmc_surface._coeffs['y0']
    opencsg_surface = opencsg.YPlane(surface_id, name, boundary, y0)

  elif openmc_surface._type == 'z-plane':
    z0 = openmc_surface._coeffs['z0']
    opencsg_surface = opencsg.ZPlane(surface_id, name, boundary, z0)

  elif openmc_surface._type == 'x-cylinder':
    y0 = openmc_surface._coeffs['y0']
    z0 = openmc_surface._coeffs['z0']
    R = openmc_surface._coeffs['R']
    opencsg_surface = opencsg.XCylinder(surface_id, name, boundary, y0, z0, R)

  elif openmc_surface._type == 'y-cylinder':
    x0 = openmc_surface._coeffs['x0']
    z0 = openmc_surface._coeffs['z0']
    R = openmc_surface._coeffs['R']
    opencsg_surface = opencsg.YCylinder(surface_id, name, boundary, x0, z0, R)

  elif openmc_surface._type == 'z-cylinder':
    x0 = openmc_surface._coeffs['x0']
    y0 = openmc_surface._coeffs['y0']
    R = openmc_surface._coeffs['R']
    opencsg_surface = opencsg.ZCylinder(surface_id, name, boundary, x0, y0, R)

  # Add the OpenMC Surface to the global collection of all OpenMC Surfaces
  OPENMC_SURFACES[surface_id] = openmc_surface

  # Add the OpenCSG Surface to the global collection of all OpenCSG Surfaces
  OPENCSG_SURFACES[surface_id] = opencsg_surface

  return opencsg_surface


def get_openmc_surface(opencsg_surface):

  if not isinstance(opencsg_surface, opencsg.Surface):
    msg = 'Unable to create an OpenMC Surface from {0} which ' \
          'is not an OpenCSG Surface'.format(opencsg_surface)
    raise ValueError(msg)

  global openmc_surface
  surface_id = opencsg_surface._id

  # If this Surface was already created, use it
  if surface_id in OPENMC_SURFACES.keys():
    return OPENMC_SURFACES[surface_id]

  # Create an OpenMC Surface to represent this OpenCSG Surface
  name = opencsg_surface._name

  # Correct for OpenMC's syntax for Surfaces dividing Cells
  boundary = opencsg_surface._boundary_type
  if boundary == 'interface':
    boundary = 'transmission'

  if opencsg_surface._type == 'plane':
    A = opencsg_surface._coeffs['A']
    B = opencsg_surface._coeffs['B']
    C = opencsg_surface._coeffs['C']
    D = opencsg_surface._coeffs['D']
    openmc_surface = openmc.Plane(surface_id, boundary, A, B, C, D, name)

  elif opencsg_surface._type == 'x-plane':
    x0 = opencsg_surface._coeffs['x0']
    openmc_surface = openmc.XPlane(surface_id, boundary, x0, name)

  elif opencsg_surface._type == 'y-plane':
    y0 = opencsg_surface._coeffs['y0']
    openmc_surface = openmc.YPlane(surface_id, boundary, y0, name)

  elif opencsg_surface._type == 'z-plane':
    z0 = opencsg_surface._coeffs['z0']
    openmc_surface = openmc.ZPlane(surface_id, boundary, z0, name)

  elif opencsg_surface._type == 'x-cylinder':
    y0 = opencsg_surface._coeffs['y0']
    z0 = opencsg_surface._coeffs['z0']
    R = opencsg_surface._coeffs['R']
    openmc_surface = openmc.XCylinder(surface_id, boundary, y0, z0, R, name)

  elif opencsg_surface._type == 'y-cylinder':
    x0 = opencsg_surface._coeffs['x0']
    z0 = opencsg_surface._coeffs['z0']
    R = opencsg_surface._coeffs['R']
    openmc_surface = openmc.YCylinder(surface_id, boundary, x0, z0, R, name)

  elif opencsg_surface._type == 'z-cylinder':
    x0 = opencsg_surface._coeffs['x0']
    y0 = opencsg_surface._coeffs['y0']
    R = opencsg_surface._coeffs['R']
    openmc_surface = openmc.ZCylinder(surface_id, boundary, x0, y0, R, name)

  else:
    msg = 'Unable to create an OpenMC Surface from an OpenCSG ' \
          'Surface of type {0} since it is not a compatible ' \
          'Surface type in OpenMC'.format(opencsg_surface._type)
    raise ValueError(msg)


  # Add the OpenMC Surface to the global collection of all OpenMC Surfaces
  OPENMC_SURFACES[surface_id] = openmc_surface

  # Add the OpenCSG Surface to the global collection of all OpenCSG Surfaces
  OPENCSG_SURFACES[surface_id] = opencsg_surface

  return openmc_surface


def get_compatible_opencsg_surfaces(opencsg_surface):

  if not isinstance(opencsg_surface, opencsg.Surface):
    msg = 'Unable to create an OpenMC Surface from {0} which ' \
          'is not an OpenCSG Surface'.format(opencsg_surface)
    raise ValueError(msg)

  global OPENMC_SURFACES
  surface_id = opencsg_surface._id

  # If this Surface was already created, use it
  if surface_id in OPENMC_SURFACES.keys():
    return OPENMC_SURFACES[surface_id]

  # Create an OpenMC Surface to represent this OpenCSG Surface
  name = opencsg_surface._name
  boundary = opencsg_surface._boundary_type

  if opencsg_surface._type == 'x-squareprism':
    y0 = opencsg_surface._coeffs['y0']
    z0 = opencsg_surface._coeffs['z0']
    R = opencsg_surface._coeffs['R']

    # Create a list of the four planes we need
    left = opencsg.YPlane(name=name, boundary=boundary, y0=y0-R)
    right = opencsg.YPlane(name=name, boundary=boundary, y0=y0+R)
    bottom = opencsg.ZPlane(name=name, boundary=boundary, z0=z0-R)
    top = opencsg.ZPlane(name=name, boundary=boundary, z0=z0+R)
    surfaces = [left, right, bottom, top]

  elif opencsg_surface._type == 'y-squareprism':
    x0 = opencsg_surface._coeffs['x0']
    z0 = opencsg_surface._coeffs['z0']
    R = opencsg_surface._coeffs['R']

    # Create a list of the four planes we need
    left = opencsg.XPlane(name=name, boundary=boundary, x0=x0-R)
    right = opencsg.XPlane(name=name, boundary=boundary, x0=x0+R)
    bottom = opencsg.ZPlane(name=name, boundary=boundary, z0=z0-R)
    top = opencsg.ZPlane(name=name, boundary=boundary, z0=z0+R)
    surfaces = [left, right, bottom, top]

  elif opencsg_surface._type == 'z-squareprism':
    x0 = opencsg_surface._coeffs['x0']
    y0 = opencsg_surface._coeffs['y0']
    R = opencsg_surface._coeffs['R']

    # Create a list of the four planes we need
    left = opencsg.XPlane(name=name, boundary=boundary, x0=x0-R)
    right = opencsg.XPlane(name=name, boundary=boundary, x0=x0+R)
    bottom = opencsg.YPlane(name=name, boundary=boundary, y0=y0-R)
    top = opencsg.YPlane(name=name, boundary=boundary, y0=y0+R)
    surfaces = [left, right, bottom, top]

  else:
    msg = 'Unable to create a compatible OpenMC Surface an OpenCSG ' \
          'Surface of type {0} since it already a compatible ' \
          'Surface type in OpenMC'.format(opencsg_surface._type)
    raise ValueError(msg)

  # Add the OpenMC Surface(s) to the global collection of all OpenMC Surfaces
  OPENMC_SURFACES[surface_id] = surfaces

  # Add the OpenCSG Surface to the global collection of all OpenCSG Surfaces
  OPENCSG_SURFACES[surface_id] = opencsg_surface

  return surfaces


def get_opencsg_cell(openmc_cell):

  if not isinstance(openmc_cell, openmc.Cell):
    msg = 'Unable to create an OpenCSG Cell from {0} which ' \
          'is not an OpenMC Cell'.format(openmc_cell)
    raise ValueError(msg)

  global OPENCSG_CELLS
  cell_id = openmc_cell._id

  # If this Cell was already created, use it
  if cell_id in OPENCSG_CELLS.keys():
    return OPENCSG_CELLS[cell_id]

  # Create an OpenCSG Cell to represent this OpenMC Cell
  name = openmc_cell._name
  opencsg_cell = opencsg.Cell(cell_id, name)

  fill = openmc_cell._fill

  if (openmc_cell._type == 'normal'):
    opencsg_cell.setFill(get_opencsg_material(fill))
  elif (openmc_cell._type == 'fill'):
    opencsg_cell.setFill(get_opencsg_universe(fill))
  else:
    opencsg_cell.setFill(get_opencsg_lattice(fill))

  surfaces = openmc_cell._surfaces

  for surface_id in surfaces.keys():
    surface = surfaces[surface_id][0]
    halfspace = surfaces[surface_id][1]
    opencsg_cell.addSurface(get_opencsg_surface(surface), halfspace)

  # Add the OpenMC Cell to the global collection of all OpenMC Cells
  OPENMC_CELLS[cell_id] = openmc_cell

  # Add the OpenCSG Cell to the global collection of all OpenCSG Cells
  OPENCSG_CELLS[cell_id] = opencsg_cell

  return opencsg_cell


def get_compatible_opencsg_cells(opencsg_cell, opencsg_surface, halfspace):

  if not isinstance(opencsg_cell, opencsg.Cell):
    msg = 'Unable to create compatible OpenMC Cell from {0} which ' \
          'is not an OpenCSG Cell'.format(opencsg_cell)
    raise ValueError(msg)

  elif not isinstance(opencsg_surface, opencsg.Surface):
    msg = 'Unable to create compatible OpenMC Cell since {0} is ' \
          'not an OpenCSG Surface'.format(opencsg_surface)
    raise ValueError(msg)

  elif not halfspace in [-1, +1]:
    msg = 'Unable to create compatible Cell since {0}' \
          'is not a +/-1 halfspace'.format(halfspace)
    raise ValueError(msg)

  # Initialize an empty list for the new compatible cells
  compatible_cells = list()

  # SquarePrism Surfaces
  if opencsg_surface._type in ['x-squareprism', 'y-squareprism', 'z-squareprism']:

    # Get the compatible Surfaces (XPlanes and YPlanes)
    compatible_surfaces = get_compatible_opencsg_surfaces(opencsg_surface)

    opencsg_cell.removeSurface(opencsg_surface)

    # If Cell is inside SquarePrism, add "inside" of Surface halfspaces
    if halfspace == -1:
      opencsg_cell.addSurface(compatible_surfaces[0], +1)
      opencsg_cell.addSurface(compatible_surfaces[1], -1)
      opencsg_cell.addSurface(compatible_surfaces[2], +1)
      opencsg_cell.addSurface(compatible_surfaces[3], -1)
      compatible_cells.append(opencsg_cell)

    # If Cell is outside SquarePrism, add "outside" of Surface halfspaces
    else:

      # Create 8 Cell clones to represent each of the disjoint planar
      # Surface halfspace intersections
      num_clones = 8

      for clone_id in range(num_clones):

        # Create a cloned OpenCSG Cell with Surfaces compatible with OpenMC
        clone = opencsg_cell.clone()
        compatible_cells.append(clone)

        # Top left subcell - add left XPlane, top YPlane
        if clone_id == 0:
          clone.addSurface(compatible_surfaces[0], -1)
          clone.addSurface(compatible_surfaces[3], +1)

        # Top center subcell - add top YPlane, left/right XPlanes
        elif clone_id == 1:
          clone.addSurface(compatible_surfaces[0], +1)
          clone.addSurface(compatible_surfaces[1], -1)
          clone.addSurface(compatible_surfaces[3], +1)

        # Top right subcell - add top YPlane, right XPlane
        elif clone_id == 2:
          clone.addSurface(compatible_surfaces[1], +1)
          clone.addSurface(compatible_surfaces[3], +1)

        # Right center subcell - add right XPlane, top/bottom YPlanes
        elif clone_id == 3:
          clone.addSurface(compatible_surfaces[1], +1)
          clone.addSurface(compatible_surfaces[3], -1)
          clone.addSurface(compatible_surfaces[2], +1)

        # Bottom right subcell - add right XPlane, bottom YPlane
        elif clone_id == 4:
          clone.addSurface(compatible_surfaces[1], +1)
          clone.addSurface(compatible_surfaces[2], -1)

        # Bottom center subcell - add bottom YPlane, left/right XPlanes
        elif clone_id == 5:
          clone.addSurface(compatible_surfaces[0], +1)
          clone.addSurface(compatible_surfaces[1], -1)
          clone.addSurface(compatible_surfaces[2], -1)

        # Bottom left subcell - add bottom YPlane, left XPlane
        elif clone_id == 6:
          clone.addSurface(compatible_surfaces[0], -1)
          clone.addSurface(compatible_surfaces[2], -1)

        # Left center subcell - add left XPlane, top/bottom YPlanes
        elif clone_id == 7:
          clone.addSurface(compatible_surfaces[0], -1)
          clone.addSurface(compatible_surfaces[3], -1)
          clone.addSurface(compatible_surfaces[2], +1)

  # Remove redundant Surfaces from the Cells
  for cell in compatible_cells:
    cell.removeRedundantSurfaces()

  # Return the list of OpenMC compatible OpenCSG Cells
  return compatible_cells


def make_opencsg_cells_compatible(opencsg_universe):

  if not isinstance(opencsg_universe, opencsg.Universe):
    msg = 'Unable to make compatible OpenCSG Cells for {0} which ' \
          'is not an OpenCSG Universe'.format(opencsg_universe)
    raise ValueError(msg)

  # Check all OpenCSG Cells in this Universe for compatibility with OpenMC
  opencsg_cells = opencsg_universe._cells

  for cell_id, opencsg_cell in opencsg_cells.items():

    # Check each of the OpenCSG Surfaces for OpenMC compatibility
    surfaces = opencsg_cell._surfaces

    for surface_id in surfaces.keys():
      surface = surfaces[surface_id][0]
      halfspace = surfaces[surface_id][1]

      # If this Surface is not compatible with OpenMC, create compatible
      # OpenCSG cells with a compatible version of this OpenCSG Surface
      if not is_opencsg_surface_compatible(surface):

        # Get the one or more OpenCSG Cells that are compatible with OpenMC
        # NOTE: This does not necessarily make the OpenCSG fully compatible.
        #       It only removes the incompatible Surface and replaces it with
        #       compatible OpenCSG Surface(s). The recursive call at the end
        #       of this block is necessary in the event that there are more
        #       incompatible Surfaces in this Cell that are not accounted for.
        cells = get_compatible_opencsg_cells(opencsg_cell, surface, halfspace)

        # Remove the non-compatible OpenCSG Cell from the Universe
        opencsg_universe.removeCell(opencsg_cell)

        # Add the compatible OpenCSG Cells to the Universe
        opencsg_universe.addCells(cells)

        # Make recursive call to look at the updated state of the
        # OpenCSG Universe and return
        return make_opencsg_cells_compatible(opencsg_universe)

  # If all OpenCSG Cells in the OpenCSG Universe are compatible, return
  return



def get_openmc_cell(opencsg_cell):

  if not isinstance(opencsg_cell, opencsg.Cell):
    msg = 'Unable to create an OpenMC Cell from {0} which ' \
          'is not an OpenCSG Cell'.format(opencsg_cell)
    raise ValueError(msg)

  global OPENMC_CELLS
  cell_id = opencsg_cell._id

  # If this Cell was already created, use it
  if cell_id in OPENMC_CELLS.keys():
    return OPENMC_CELLS[cell_id]

  # Create an OpenCSG Cell to represent this OpenMC Cell
  name = opencsg_cell._name
  openmc_cell = openmc.Cell(cell_id, name)

  fill = opencsg_cell._fill

  if (opencsg_cell._type == 'universe'):
    openmc_cell.setFill(get_openmc_universe(fill))
  elif (opencsg_cell._type == 'lattice'):
    openmc_cell.setFill(get_openmc_lattice(fill))
  else:
    openmc_cell.setFill(get_openmc_material(fill))

  surfaces = opencsg_cell._surfaces

  for surface_id in surfaces.keys():
    surface = surfaces[surface_id][0]
    halfspace = surfaces[surface_id][1]
    openmc_cell.addSurface(get_openmc_surface(surface), halfspace)

  # Add the OpenMC Cell to the global collection of all OpenMC Cells
  OPENMC_CELLS[cell_id] = openmc_cell

  # Add the OpenCSG Cell to the global collection of all OpenCSG Cells
  OPENCSG_CELLS[cell_id] = opencsg_cell

  return openmc_cell



def get_opencsg_universe(openmc_universe):

  if not isinstance(openmc_universe, openmc.Universe):
    msg = 'Unable to create an OpenCSG Universe from {0} which ' \
          'is not an OpenMC Universe'.format(openmc_universe)
    raise ValueError(msg)

  global OPENCSG_UNIVERSES
  universe_id = openmc_universe._id

  # If this Universe was already created, use it
  if universe_id in OPENCSG_UNIVERSES.keys():
    return OPENCSG_UNIVERSES[universe_id]

  # Create an OpenCSG Universe to represent this OpenMC Universe
  name = openmc_universe._name
  opencsg_universe = opencsg.Universe(universe_id, name)

  # Convert all OpenMC Cells in this Universe to OpenCSG Cells
  openmc_cells = openmc_universe._cells

  for cell_id, openmc_cell in openmc_cells.items():
    opencsg_cell = get_opencsg_cell(openmc_cell)
    opencsg_universe.addCell(opencsg_cell)

  # Add the OpenMC Universe to the global collection of all OpenMC Universes
  OPENMC_UNIVERSES[universe_id] = openmc_universe

  # Add the OpenCSG Universe to the global collection of all OpenCSG Universes
  OPENCSG_UNIVERSES[universe_id] = opencsg_universe

  return opencsg_universe


def get_openmc_universe(opencsg_universe):

  if not isinstance(opencsg_universe, opencsg.Universe):
    msg = 'Unable to create an OpenMC Universe from {0} which ' \
          'is not an OpenCSG Universe'.format(opencsg_universe)
    raise ValueError(msg)

  global OPENMC_UNIVERSES
  universe_id = opencsg_universe._id

  # If this Universe was already created, use it
  if universe_id in OPENMC_UNIVERSES.keys():
    return OPENMC_UNIVERSES[universe_id]

  # Make all OpenCSG Cells and Surfaces in this Universe compatible with OpenMC
  make_opencsg_cells_compatible(opencsg_universe)

  # Create an OpenMC Universe to represent this OpenCSg Universe
  name = opencsg_universe._name
  openmc_universe = openmc.Universe(universe_id, name)

  # Convert all OpenCSG Cells in this Universe to OpenMC Cells
  opencsg_cells = opencsg_universe._cells

  for cell_id, opencsg_cell in opencsg_cells.items():
    openmc_cell = get_openmc_cell(opencsg_cell)
    openmc_universe.addCell(openmc_cell)

  # Add the OpenMC Universe to the global collection of all OpenMC Universes
  OPENMC_UNIVERSES[universe_id] = openmc_universe

  # Add the OpenCSG Universe to the global collection of all OpenCSG Universes
  OPENCSG_UNIVERSES[universe_id] = opencsg_universe

  return openmc_universe


def get_opencsg_lattice(openmc_lattice):

  if not isinstance(openmc_lattice, openmc.Lattice):
    msg = 'Unable to create an OpenCSG Lattice from {0} which ' \
          'is not an OpenMC Lattice'.format(openmc_lattice)
    raise ValueError(msg)

  global OPENCSG_LATTICES
  lattice_id = openmc_lattice._id

  # If this Lattice was already created, use it
  if lattice_id in OPENCSG_LATTICES.keys():
    return OPENCSG_LATTICES[lattice_id]

  # Create an OpenCSG Lattice to represent this OpenMC Lattice
  name = openmc_lattice._name
  dimension = openmc_lattice._dimension
  width = openmc_lattice._width
  lower_left = openmc_lattice._lower_left
  universes = openmc_lattice._universes

  # Initialize an empty array for the OpenCSG nested Universes in this Lattice
  universe_array = np.ndarray(tuple(np.array(dimension)[::-1]), \
                              dtype=opencsg.Universe)

  # Create OpenCSG Universes for each unique nested Universe in this Lattice
  unique_universes = openmc_lattice.getUniqueUniverses()

  for universe_id, universe in unique_universes.items():
    unique_universes[universe_id] = get_opencsg_universe(universe)

  # Build the nested Universe array
  for z in range(dimension[2]):
    for y in range(dimension[1]):
      for x in range(dimension[0]):
        universe_id = universes[x][y][z]._id
        universe_array[z][y][x] = unique_universes[universe_id]

  opencsg_lattice = opencsg.Lattice(lattice_id, name)
  opencsg_lattice.setDimension(dimension)
  opencsg_lattice.setWidth(width)
  opencsg_lattice.setUniverses(universe_array)

  offset = np.array(lower_left, dtype=np.float64) - \
           ((np.array(width, dtype=np.float64) * \
             np.array(dimension, dtype=np.float64))) / -2.0
  opencsg_lattice.setOffset(offset)

  # Add the OpenMC Lattice to the global collection of all OpenMC Lattices
  OPENMC_LATTICES[lattice_id] = openmc_lattice

  # Add the OpenCSG Lattice to the global collection of all OpenCSG Lattices
  OPENCSG_LATTICES[lattice_id] = opencsg_lattice

  return opencsg_lattice


def get_openmc_lattice(opencsg_lattice):

  if not isinstance(opencsg_lattice, opencsg.Lattice):
    msg = 'Unable to create an OpenMC Lattice from {0} which ' \
          'is not an OpenCSG Lattice'.format(opencsg_lattice)
    raise ValueError(msg)

  global OPENMC_LATTICES
  lattice_id = opencsg_lattice._id

  # If this Lattice was already created, use it
  if lattice_id in OPENMC_LATTICES.keys():
    return OPENMC_LATTICES[lattice_id]

  dimension = opencsg_lattice._dimension
  width = opencsg_lattice._width
  offset = opencsg_lattice._offset
  universes = opencsg_lattice._universes

  # Initialize an empty array for the OpenMC nested Universes in this Lattice
  universe_array = np.ndarray(tuple(np.array(dimension)), \
                              dtype=openmc.Universe)

  # Create OpenMC Universes for each unique nested Universe in this Lattice
  unique_universes = opencsg_lattice.getUniqueUniverses()

  for universe_id, universe in unique_universes.items():
    unique_universes[universe_id] = get_openmc_universe(universe)

  # Build the nested Universe array
  for z in range(dimension[2]):
    for y in range(dimension[1]):
      for x in range(dimension[0]):
        universe_id = universes[z][y][x]._id
        universe_array[x][y][z] = unique_universes[universe_id]

  openmc_lattice = openmc.Lattice(lattice_id=lattice_id)
  openmc_lattice.setDimension(dimension)
  openmc_lattice.setWidth(width)
  openmc_lattice.setUniverses(universe_array)

  lower_left = np.array(offset, dtype=np.float64) + \
               ((np.array(width, dtype=np.float64) * \
                 np.array(dimension, dtype=np.float64))) / -2.0
  openmc_lattice.setLowerLeft(lower_left)

  # Add the OpenMC Lattice to the global collection of all OpenMC Lattices
  OPENMC_LATTICES[lattice_id] = openmc_lattice

  # Add the OpenCSG Lattice to the global collection of all OpenCSG Lattices
  OPENCSG_LATTICES[lattice_id] = opencsg_lattice

  return openmc_lattice


def get_opencsg_geometry(openmc_geometry):

  if not isinstance(openmc_geometry, openmc.Geometry):
    msg = 'Unable to get OpenCSG geometry from {0} which is ' \
          'not an OpenMC Geometry object'.format(openmc_geometry)
    raise ValueError(msg)

  # Clear dictionaries and auto-generated IDs
  OPENMC_SURFACES.clear()
  OPENCSG_SURFACES.clear()
  OPENMC_CELLS.clear()
  OPENCSG_CELLS.clear()
  OPENMC_UNIVERSES.clear()
  OPENCSG_UNIVERSES.clear()
  OPENMC_LATTICES.clear()
  OPENCSG_LATTICES.clear()

  opencsg.geometry.reset_auto_ids()

  openmc_root_universe = openmc_geometry._root_universe
  opencsg_root_universe = get_opencsg_universe(openmc_root_universe)

  opencsg_geometry = opencsg.Geometry()
  opencsg_geometry.setRootUniverse(opencsg_root_universe)
  opencsg_geometry.initializeCellOffsets()

  return opencsg_geometry


def get_openmc_geometry(opencsg_geometry):

  if not isinstance(opencsg_geometry, opencsg.Geometry):
    msg = 'Unable to get OpenMC geometry from {0} which is ' \
          'not an OpenCSG Geometry object'.format(opencsg_geometry)
    raise ValueError(msg)

  # Deep copy the goemetry since it may be modified to make all Surfaces
  # compatible with OpenMC's specifications
  opencsg_geometry = copy.deepcopy(opencsg_geometry)

  # Update Cell bounding boxes in Geometry
  opencsg_geometry.updateBoundingBoxes()

  # Clear dictionaries and auto-generated ID
  OPENMC_SURFACES.clear()
  OPENCSG_SURFACES.clear()
  OPENMC_CELLS.clear()
  OPENCSG_CELLS.clear()
  OPENMC_UNIVERSES.clear()
  OPENCSG_UNIVERSES.clear()
  OPENMC_LATTICES.clear()
  OPENCSG_LATTICES.clear()

  openmc.reset_auto_ids()

  opencsg_root_universe = opencsg_geometry._root_universe
  openmc_root_universe = get_openmc_universe(opencsg_root_universe)

  openmc_geometry = openmc.Geometry()
  openmc_geometry.setRootUniverse(openmc_root_universe)

  return openmc_geometry
