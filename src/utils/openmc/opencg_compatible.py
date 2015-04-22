import copy

import numpy as np
import opencg

import openmc


# A dictionary of all OpenMC Materials created
# Keys   - Material IDs
# Values - Materials
OPENMC_MATERIALS = {}

# A dictionary of all OpenCG Materials created
# Keys   - Material IDs
# Values - Materials
OPENCG_MATERIALS = {}

# A dictionary of all OpenMC Surfaces created
# Keys   - Surface IDs
# Values - Surfaces
OPENMC_SURFACES = {}

# A dictionary of all OpenCG Surfaces created
# Keys   - Surface IDs
# Values - Surfaces
OPENCG_SURFACES = {}

# A dictionary of all OpenMC Cells created
# Keys   - Cell IDs
# Values - Cells
OPENMC_CELLS = {}

# A dictionary of all OpenCG Cells created
# Keys   - Cell IDs
# Values - Cells
OPENCG_CELLS = {}

# A dictionary of all OpenMC Universes created
# Keys   - Universes IDs
# Values - Universes
OPENMC_UNIVERSES = {}

# A dictionary of all OpenCG Universes created
# Keys   - Universes IDs
# Values - Universes
OPENCG_UNIVERSES = {}

# A dictionary of all OpenMC Lattices created
# Keys   - Lattice IDs
# Values - Lattices
OPENMC_LATTICES = {}

# A dictionary of all OpenCG Lattices created
# Keys   - Lattice IDs
# Values - Lattices
OPENCG_LATTICES = {}



def get_opencg_material(openmc_material):

    if not isinstance(openmc_material, openmc.Material):
        msg = 'Unable to create an OpenCG Material from {0} ' \
              'which is not an OpenMC Material'.format(openmc_material)
        raise ValueError(msg)

    global OPENCG_MATERIALS
    material_id = openmc_material._id

    # If this Material was already created, use it
    if material_id in OPENCG_MATERIALS:
        return OPENCG_MATERIALS[material_id]

    # Create an OpenCG Material to represent this OpenMC Material
    name = openmc_material._name
    opencg_material = opencg.Material(material_id=material_id, name=name)

    # Add the OpenMC Material to the global collection of all OpenMC Materials
    OPENMC_MATERIALS[material_id] = openmc_material

    # Add the OpenCG Material to the global collection of all OpenCG Materials
    OPENCG_MATERIALS[material_id] = opencg_material

    return opencg_material


def get_openmc_material(opencg_material):

    if not isinstance(opencg_material, opencg.Material):
        msg = 'Unable to create an OpenMC Material from {0} ' \
              'which is not an OpenCG Material'.format(opencg_material)
        raise ValueError(msg)

    global OPENMC_MATERIALS
    material_id = opencg_material._id

    # If this Material was already created, use it
    if material_id in OPENMC_MATERIALS:
        return OPENMC_MATERIALS[material_id]

    # Create an OpenMC Material to represent this OpenCG Material
    name = opencg_material._name
    openmc_material = openmc.Material(material_id=material_id, name=name)

    # Add the OpenMC Material to the global collection of all OpenMC Materials
    OPENMC_MATERIALS[material_id] = openmc_material

    # Add the OpenCG Material to the global collection of all OpenCG Materials
    OPENCG_MATERIALS[material_id] = opencg_material

    return openmc_material


def is_opencg_surface_compatible(opencg_surface):

    if not isinstance(opencg_surface, opencg.Surface):
        msg = 'Unable to check if OpenCG Surface is compatible' \
              'since {0} is not a Surface'.format(opencg_surface)
        raise ValueError(msg)

    if opencg_surface._type in ['x-squareprism',
                                 'y-squareprism', 'z-squareprism']:
        return False
    else:
        return True


def get_opencg_surface(openmc_surface):

    if not isinstance(openmc_surface, openmc.Surface):
        msg = 'Unable to create an OpenCG Surface from {0} ' \
              'which is not an OpenMC Surface'.format(openmc_surface)
        raise ValueError(msg)

    global OPENCG_SURFACES
    surface_id = openmc_surface._id

    # If this Material was already created, use it
    if surface_id in OPENCG_SURFACES:
        return OPENCG_SURFACES[surface_id]

    # Create an OpenCG Surface to represent this OpenMC Surface
    name = openmc_surface._name

    # Correct for OpenMC's syntax for Surfaces dividing Cells
    boundary = openmc_surface._boundary_type
    if boundary == 'transmission':
        boundary = 'interface'

    opencg_surface = None

    if openmc_surface._type == 'plane':
        A = openmc_surface._coeffs['A']
        B = openmc_surface._coeffs['B']
        C = openmc_surface._coeffs['C']
        D = openmc_surface._coeffs['D']
        opencg_surface = opencg.Plane(surface_id, name, boundary, A, B, C, D)

    elif openmc_surface._type == 'x-plane':
        x0 = openmc_surface._coeffs['x0']
        opencg_surface = opencg.XPlane(surface_id, name, boundary, x0)

    elif openmc_surface._type == 'y-plane':
        y0 = openmc_surface._coeffs['y0']
        opencg_surface = opencg.YPlane(surface_id, name, boundary, y0)

    elif openmc_surface._type == 'z-plane':
        z0 = openmc_surface._coeffs['z0']
        opencg_surface = opencg.ZPlane(surface_id, name, boundary, z0)

    elif openmc_surface._type == 'x-cylinder':
        y0 = openmc_surface._coeffs['y0']
        z0 = openmc_surface._coeffs['z0']
        R = openmc_surface._coeffs['R']
        opencg_surface = opencg.XCylinder(surface_id, name,
                                            boundary, y0, z0, R)

    elif openmc_surface._type == 'y-cylinder':
        x0 = openmc_surface._coeffs['x0']
        z0 = openmc_surface._coeffs['z0']
        R = openmc_surface._coeffs['R']
        opencg_surface = opencg.YCylinder(surface_id, name,
                                            boundary, x0, z0, R)

    elif openmc_surface._type == 'z-cylinder':
        x0 = openmc_surface._coeffs['x0']
        y0 = openmc_surface._coeffs['y0']
        R = openmc_surface._coeffs['R']
        opencg_surface = opencg.ZCylinder(surface_id, name,
                                            boundary, x0, y0, R)

    # Add the OpenMC Surface to the global collection of all OpenMC Surfaces
    OPENMC_SURFACES[surface_id] = openmc_surface

    # Add the OpenCG Surface to the global collection of all OpenCG Surfaces
    OPENCG_SURFACES[surface_id] = opencg_surface

    return opencg_surface


def get_openmc_surface(opencg_surface):

    if not isinstance(opencg_surface, opencg.Surface):
        msg = 'Unable to create an OpenMC Surface from {0} which ' \
              'is not an OpenCG Surface'.format(opencg_surface)
        raise ValueError(msg)

    global openmc_surface
    surface_id = opencg_surface._id

    # If this Surface was already created, use it
    if surface_id in OPENMC_SURFACES:
        return OPENMC_SURFACES[surface_id]

    # Create an OpenMC Surface to represent this OpenCG Surface
    name = opencg_surface._name

    # Correct for OpenMC's syntax for Surfaces dividing Cells
    boundary = opencg_surface._boundary_type
    if boundary == 'interface':
        boundary = 'transmission'

    if opencg_surface._type == 'plane':
        A = opencg_surface._coeffs['A']
        B = opencg_surface._coeffs['B']
        C = opencg_surface._coeffs['C']
        D = opencg_surface._coeffs['D']
        openmc_surface = openmc.Plane(surface_id, boundary, A, B, C, D, name)

    elif opencg_surface._type == 'x-plane':
        x0 = opencg_surface._coeffs['x0']
        openmc_surface = openmc.XPlane(surface_id, boundary, x0, name)

    elif opencg_surface._type == 'y-plane':
        y0 = opencg_surface._coeffs['y0']
        openmc_surface = openmc.YPlane(surface_id, boundary, y0, name)

    elif opencg_surface._type == 'z-plane':
        z0 = opencg_surface._coeffs['z0']
        openmc_surface = openmc.ZPlane(surface_id, boundary, z0, name)

    elif opencg_surface._type == 'x-cylinder':
        y0 = opencg_surface._coeffs['y0']
        z0 = opencg_surface._coeffs['z0']
        R = opencg_surface._coeffs['R']
        openmc_surface = openmc.XCylinder(surface_id, boundary, y0, z0, R, name)

    elif opencg_surface._type == 'y-cylinder':
        x0 = opencg_surface._coeffs['x0']
        z0 = opencg_surface._coeffs['z0']
        R = opencg_surface._coeffs['R']
        openmc_surface = openmc.YCylinder(surface_id, boundary, x0, z0, R, name)

    elif opencg_surface._type == 'z-cylinder':
        x0 = opencg_surface._coeffs['x0']
        y0 = opencg_surface._coeffs['y0']
        R = opencg_surface._coeffs['R']
        openmc_surface = openmc.ZCylinder(surface_id, boundary, x0, y0, R, name)

    else:
        msg = 'Unable to create an OpenMC Surface from an OpenCG ' \
              'Surface of type {0} since it is not a compatible ' \
              'Surface type in OpenMC'.format(opencg_surface._type)
        raise ValueError(msg)


    # Add the OpenMC Surface to the global collection of all OpenMC Surfaces
    OPENMC_SURFACES[surface_id] = openmc_surface

    # Add the OpenCG Surface to the global collection of all OpenCG Surfaces
    OPENCG_SURFACES[surface_id] = opencg_surface

    return openmc_surface


def get_compatible_opencg_surfaces(opencg_surface):

    if not isinstance(opencg_surface, opencg.Surface):
        msg = 'Unable to create an OpenMC Surface from {0} which ' \
              'is not an OpenCG Surface'.format(opencg_surface)
        raise ValueError(msg)

    global OPENMC_SURFACES
    surface_id = opencg_surface._id

    # If this Surface was already created, use it
    if surface_id in OPENMC_SURFACES:
        return OPENMC_SURFACES[surface_id]

    # Create an OpenMC Surface to represent this OpenCG Surface
    name = opencg_surface._name
    boundary = opencg_surface._boundary_type

    if opencg_surface._type == 'x-squareprism':
        y0 = opencg_surface._coeffs['y0']
        z0 = opencg_surface._coeffs['z0']
        R = opencg_surface._coeffs['R']

        # Create a list of the four planes we need
        left = opencg.YPlane(name=name, boundary=boundary, y0=y0-R)
        right = opencg.YPlane(name=name, boundary=boundary, y0=y0+R)
        bottom = opencg.ZPlane(name=name, boundary=boundary, z0=z0-R)
        top = opencg.ZPlane(name=name, boundary=boundary, z0=z0+R)
        surfaces = [left, right, bottom, top]

    elif opencg_surface._type == 'y-squareprism':
        x0 = opencg_surface._coeffs['x0']
        z0 = opencg_surface._coeffs['z0']
        R = opencg_surface._coeffs['R']

        # Create a list of the four planes we need
        left = opencg.XPlane(name=name, boundary=boundary, x0=x0-R)
        right = opencg.XPlane(name=name, boundary=boundary, x0=x0+R)
        bottom = opencg.ZPlane(name=name, boundary=boundary, z0=z0-R)
        top = opencg.ZPlane(name=name, boundary=boundary, z0=z0+R)
        surfaces = [left, right, bottom, top]

    elif opencg_surface._type == 'z-squareprism':
        x0 = opencg_surface._coeffs['x0']
        y0 = opencg_surface._coeffs['y0']
        R = opencg_surface._coeffs['R']

        # Create a list of the four planes we need
        left = opencg.XPlane(name=name, boundary=boundary, x0=x0-R)
        right = opencg.XPlane(name=name, boundary=boundary, x0=x0+R)
        bottom = opencg.YPlane(name=name, boundary=boundary, y0=y0-R)
        top = opencg.YPlane(name=name, boundary=boundary, y0=y0+R)
        surfaces = [left, right, bottom, top]

    else:
        msg = 'Unable to create a compatible OpenMC Surface an OpenCG ' \
              'Surface of type {0} since it already a compatible ' \
              'Surface type in OpenMC'.format(opencg_surface._type)
        raise ValueError(msg)

    # Add the OpenMC Surface(s) to the global collection of all OpenMC Surfaces
    OPENMC_SURFACES[surface_id] = surfaces

    # Add the OpenCG Surface to the global collection of all OpenCG Surfaces
    OPENCG_SURFACES[surface_id] = opencg_surface

    return surfaces


def get_opencg_cell(openmc_cell):

    if not isinstance(openmc_cell, openmc.Cell):
        msg = 'Unable to create an OpenCG Cell from {0} which ' \
              'is not an OpenMC Cell'.format(openmc_cell)
        raise ValueError(msg)

    global OPENCG_CELLS
    cell_id = openmc_cell._id

    # If this Cell was already created, use it
    if cell_id in OPENCG_CELLS:
        return OPENCG_CELLS[cell_id]

    # Create an OpenCG Cell to represent this OpenMC Cell
    name = openmc_cell._name
    opencg_cell = opencg.Cell(cell_id, name)

    fill = openmc_cell._fill

    if (openmc_cell._type == 'normal'):
        opencg_cell.setFill(get_opencg_material(fill))
    elif (openmc_cell._type == 'fill'):
        opencg_cell.setFill(get_opencg_universe(fill))
    else:
        opencg_cell.setFill(get_opencg_lattice(fill))

    if openmc_cell._rotation is not None:
        opencg_cell.setRotation(openmc_cell._rotation)

    if openmc_cell._translation is not None:
        opencg_cell.setTranslation(openmc_cell._translation)

    surfaces = openmc_cell._surfaces

    for surface_id in surfaces:
        surface = surfaces[surface_id][0]
        halfspace = surfaces[surface_id][1]
        opencg_cell.addSurface(get_opencg_surface(surface), halfspace)

    # Add the OpenMC Cell to the global collection of all OpenMC Cells
    OPENMC_CELLS[cell_id] = openmc_cell

    # Add the OpenCG Cell to the global collection of all OpenCG Cells
    OPENCG_CELLS[cell_id] = opencg_cell

    return opencg_cell


def get_compatible_opencg_cells(opencg_cell, opencg_surface, halfspace):

    if not isinstance(opencg_cell, opencg.Cell):
        msg = 'Unable to create compatible OpenMC Cell from {0} which ' \
              'is not an OpenCG Cell'.format(opencg_cell)
        raise ValueError(msg)

    elif not isinstance(opencg_surface, opencg.Surface):
        msg = 'Unable to create compatible OpenMC Cell since {0} is ' \
              'not an OpenCG Surface'.format(opencg_surface)
        raise ValueError(msg)

    elif not halfspace in [-1, +1]:
        msg = 'Unable to create compatible Cell since {0}' \
              'is not a +/-1 halfspace'.format(halfspace)
        raise ValueError(msg)

    # Initialize an empty list for the new compatible cells
    compatible_cells = []

    # SquarePrism Surfaces
    if opencg_surface._type in ['x-squareprism',
                                 'y-squareprism', 'z-squareprism']:

        # Get the compatible Surfaces (XPlanes and YPlanes)
        compatible_surfaces = get_compatible_opencg_surfaces(opencg_surface)

        opencg_cell.removeSurface(opencg_surface)

        # If Cell is inside SquarePrism, add "inside" of Surface halfspaces
        if halfspace == -1:
            opencg_cell.addSurface(compatible_surfaces[0], +1)
            opencg_cell.addSurface(compatible_surfaces[1], -1)
            opencg_cell.addSurface(compatible_surfaces[2], +1)
            opencg_cell.addSurface(compatible_surfaces[3], -1)
            compatible_cells.append(opencg_cell)

        # If Cell is outside SquarePrism, add "outside" of Surface halfspaces
        else:

            # Create 8 Cell clones to represent each of the disjoint planar
            # Surface halfspace intersections
            num_clones = 8

            for clone_id in range(num_clones):

                # Create a cloned OpenCG Cell with Surfaces compatible with OpenMC
                clone = opencg_cell.clone()
                compatible_cells.append(clone)

                # Top left subcell - add left XPlane, top YPlane
                if clone_id == 0:
                    clone.add_surface(compatible_surfaces[0], -1)
                    clone.add_surface(compatible_surfaces[3], +1)

                # Top center subcell - add top YPlane, left/right XPlanes
                elif clone_id == 1:
                    clone.add_surface(compatible_surfaces[0], +1)
                    clone.add_surface(compatible_surfaces[1], -1)
                    clone.add_surface(compatible_surfaces[3], +1)

                # Top right subcell - add top YPlane, right XPlane
                elif clone_id == 2:
                    clone.add_surface(compatible_surfaces[1], +1)
                    clone.add_surface(compatible_surfaces[3], +1)

                # Right center subcell - add right XPlane, top/bottom YPlanes
                elif clone_id == 3:
                    clone.add_surface(compatible_surfaces[1], +1)
                    clone.add_surface(compatible_surfaces[3], -1)
                    clone.add_surface(compatible_surfaces[2], +1)

                # Bottom right subcell - add right XPlane, bottom YPlane
                elif clone_id == 4:
                    clone.add_surface(compatible_surfaces[1], +1)
                    clone.add_surface(compatible_surfaces[2], -1)

                # Bottom center subcell - add bottom YPlane, left/right XPlanes
                elif clone_id == 5:
                    clone.add_surface(compatible_surfaces[0], +1)
                    clone.add_surface(compatible_surfaces[1], -1)
                    clone.add_surface(compatible_surfaces[2], -1)

                # Bottom left subcell - add bottom YPlane, left XPlane
                elif clone_id == 6:
                    clone.add_surface(compatible_surfaces[0], -1)
                    clone.add_surface(compatible_surfaces[2], -1)

                # Left center subcell - add left XPlane, top/bottom YPlanes
                elif clone_id == 7:
                    clone.add_surface(compatible_surfaces[0], -1)
                    clone.add_surface(compatible_surfaces[3], -1)
                    clone.add_surface(compatible_surfaces[2], +1)

    # Remove redundant Surfaces from the Cells
    for cell in compatible_cells:
        cell.removeRedundantSurfaces()

    # Return the list of OpenMC compatible OpenCG Cells
    return compatible_cells


def make_opencg_cells_compatible(opencg_universe):

    if not isinstance(opencg_universe, opencg.Universe):
        msg = 'Unable to make compatible OpenCG Cells for {0} which ' \
              'is not an OpenCG Universe'.format(opencg_universe)
        raise ValueError(msg)

    # Check all OpenCG Cells in this Universe for compatibility with OpenMC
    opencg_cells = opencg_universe._cells

    for cell_id, opencg_cell in opencg_cells.items():

        # Check each of the OpenCG Surfaces for OpenMC compatibility
        surfaces = opencg_cell._surfaces

        for surface_id in surfaces:
            surface = surfaces[surface_id][0]
            halfspace = surfaces[surface_id][1]

            # If this Surface is not compatible with OpenMC, create compatible
            # OpenCG cells with a compatible version of this OpenCG Surface
            if not is_opencg_surface_compatible(surface):

                # Get one or more OpenCG Cells that are compatible with OpenMC
                # NOTE: This does not necessarily make OpenCG fully compatible.
                # It only removes the incompatible Surface and replaces it with
                # compatible OpenCG Surface(s). The recursive call at the end
                # of this block is necessary in the event that there are more
                # incompatible Surfaces in this Cell that are not accounted for.
                cells = get_compatible_opencg_cells(opencg_cell,
                                                     surface, halfspace)

                # Remove the non-compatible OpenCG Cell from the Universe
                opencg_universe.removeCell(opencg_cell)

                # Add the compatible OpenCG Cells to the Universe
                opencg_universe.addCells(cells)

                # Make recursive call to look at the updated state of the
                # OpenCG Universe and return
                return make_opencg_cells_compatible(opencg_universe)

    # If all OpenCG Cells in the OpenCG Universe are compatible, return
    return



def get_openmc_cell(opencg_cell):

    if not isinstance(opencg_cell, opencg.Cell):
        msg = 'Unable to create an OpenMC Cell from {0} which ' \
              'is not an OpenCG Cell'.format(opencg_cell)
        raise ValueError(msg)

    global OPENMC_CELLS
    cell_id = opencg_cell._id

    # If this Cell was already created, use it
    if cell_id in OPENMC_CELLS:
        return OPENMC_CELLS[cell_id]

    # Create an OpenCG Cell to represent this OpenMC Cell
    name = opencg_cell._name
    openmc_cell = openmc.Cell(cell_id, name)

    fill = opencg_cell._fill

    if (opencg_cell._type == 'universe'):
        openmc_cell.fill = get_openmc_universe(fill)
    elif (opencg_cell._type == 'lattice'):
        openmc_cell.fill = get_openmc_lattice(fill)
    else:
        openmc_cell.fill = get_openmc_material(fill)

    if opencg_cell._rotation:
        rotation = np.asarray(opencg_cell._rotation, dtype=np.int)
        openmc_cell.rotation = rotation

    if opencg_cell._translation:
        translation = np.asarray(opencg_cell._translation, dtype=np.float64)
        openmc_cell.setTranslation(translation)

    surfaces = opencg_cell._surfaces

    for surface_id in surfaces:
        surface = surfaces[surface_id][0]
        halfspace = surfaces[surface_id][1]
        openmc_cell.add_surface(get_openmc_surface(surface), halfspace)

    # Add the OpenMC Cell to the global collection of all OpenMC Cells
    OPENMC_CELLS[cell_id] = openmc_cell

    # Add the OpenCG Cell to the global collection of all OpenCG Cells
    OPENCG_CELLS[cell_id] = opencg_cell

    return openmc_cell



def get_opencg_universe(openmc_universe):

    if not isinstance(openmc_universe, openmc.Universe):
        msg = 'Unable to create an OpenCG Universe from {0} which ' \
              'is not an OpenMC Universe'.format(openmc_universe)
        raise ValueError(msg)

    global OPENCG_UNIVERSES
    universe_id = openmc_universe._id

    # If this Universe was already created, use it
    if universe_id in OPENCG_UNIVERSES:
        return OPENCG_UNIVERSES[universe_id]

    # Create an OpenCG Universe to represent this OpenMC Universe
    name = openmc_universe._name
    opencg_universe = opencg.Universe(universe_id, name)

    # Convert all OpenMC Cells in this Universe to OpenCG Cells
    openmc_cells = openmc_universe._cells

    for cell_id, openmc_cell in openmc_cells.items():
        opencg_cell = get_opencg_cell(openmc_cell)
        opencg_universe.addCell(opencg_cell)

    # Add the OpenMC Universe to the global collection of all OpenMC Universes
    OPENMC_UNIVERSES[universe_id] = openmc_universe

    # Add the OpenCG Universe to the global collection of all OpenCG Universes
    OPENCG_UNIVERSES[universe_id] = opencg_universe

    return opencg_universe


def get_openmc_universe(opencg_universe):

    if not isinstance(opencg_universe, opencg.Universe):
        msg = 'Unable to create an OpenMC Universe from {0} which ' \
              'is not an OpenCG Universe'.format(opencg_universe)
        raise ValueError(msg)

    global OPENMC_UNIVERSES
    universe_id = opencg_universe._id

    # If this Universe was already created, use it
    if universe_id in OPENMC_UNIVERSES:
        return OPENMC_UNIVERSES[universe_id]

    # Make all OpenCG Cells and Surfaces in this Universe compatible with OpenMC
    make_opencg_cells_compatible(opencg_universe)

    # Create an OpenMC Universe to represent this OpenCSg Universe
    name = opencg_universe._name
    openmc_universe = openmc.Universe(universe_id, name)

    # Convert all OpenCG Cells in this Universe to OpenMC Cells
    opencg_cells = opencg_universe._cells

    for cell_id, opencg_cell in opencg_cells.items():
        openmc_cell = get_openmc_cell(opencg_cell)
        openmc_universe.add_cell(openmc_cell)

    # Add the OpenMC Universe to the global collection of all OpenMC Universes
    OPENMC_UNIVERSES[universe_id] = openmc_universe

    # Add the OpenCG Universe to the global collection of all OpenCG Universes
    OPENCG_UNIVERSES[universe_id] = opencg_universe

    return openmc_universe


def get_opencg_lattice(openmc_lattice):

    if not isinstance(openmc_lattice, openmc.Lattice):
        msg = 'Unable to create an OpenCG Lattice from {0} which ' \
              'is not an OpenMC Lattice'.format(openmc_lattice)
        raise ValueError(msg)

    global OPENCG_LATTICES
    lattice_id = openmc_lattice._id

    # If this Lattice was already created, use it
    if lattice_id in OPENCG_LATTICES:
        return OPENCG_LATTICES[lattice_id]

    # Create an OpenCG Lattice to represent this OpenMC Lattice
    name = openmc_lattice.name
    dimension = openmc_lattice.dimension
    pitch = openmc_lattice.pitch
    lower_left = openmc_lattice.lower_left
    universes = openmc_lattice.universes

    if len(pitch) == 2:
        new_pitch = np.ones(3, dtype=np.float64)
        new_pitch[:2] = pitch
        pitch = new_pitch

    if len(lower_left) == 2:
        new_lower_left = np.ones(3, dtype=np.float64)
        new_lower_left[:2] = lower_left
        lower_left = new_lower_left

    # Initialize an empty array for the OpenCG nested Universes in this Lattice
    universe_array = np.ndarray(tuple(np.array(dimension)[::-1]), \
                                dtype=opencg.Universe)

    # Create OpenCG Universes for each unique nested Universe in this Lattice
    unique_universes = openmc_lattice.get_unique_universes()

    for universe_id, universe in unique_universes.items():
        unique_universes[universe_id] = get_opencg_universe(universe)

    # Build the nested Universe array
    for z in range(dimension[2]):
        for y in range(dimension[1]):
            for x in range(dimension[0]):
                universe_id = universes[x][dimension[1]-y-1][z]._id
                universe_array[z][y][x] = unique_universes[universe_id]

    opencg_lattice = opencg.Lattice(lattice_id, name)
    opencg_lattice.setDimension(dimension)
    opencg_lattice.setWidth(pitch)
    opencg_lattice.setUniverses(universe_array)

    offset = np.array(lower_left, dtype=np.float64) - \
                     ((np.array(pitch, dtype=np.float64) * \
                       np.array(dimension, dtype=np.float64))) / -2.0
    opencg_lattice.setOffset(offset)

    # Add the OpenMC Lattice to the global collection of all OpenMC Lattices
    OPENMC_LATTICES[lattice_id] = openmc_lattice

    # Add the OpenCG Lattice to the global collection of all OpenCG Lattices
    OPENCG_LATTICES[lattice_id] = opencg_lattice

    return opencg_lattice


def get_openmc_lattice(opencg_lattice):

    if not isinstance(opencg_lattice, opencg.Lattice):
        msg = 'Unable to create an OpenMC Lattice from {0} which ' \
              'is not an OpenCG Lattice'.format(opencg_lattice)
        raise ValueError(msg)

    global OPENMC_LATTICES
    lattice_id = opencg_lattice._id

    # If this Lattice was already created, use it
    if lattice_id in OPENMC_LATTICES:
        return OPENMC_LATTICES[lattice_id]

    dimension = opencg_lattice._dimension
    width = opencg_lattice._width
    offset = opencg_lattice._offset
    universes = opencg_lattice._universes

    # Initialize an empty array for the OpenMC nested Universes in this Lattice
    universe_array = np.ndarray(tuple(np.array(dimension)), \
                                dtype=openmc.Universe)

    # Create OpenMC Universes for each unique nested Universe in this Lattice
    unique_universes = opencg_lattice.getUniqueUniverses()

    for universe_id, universe in unique_universes.items():
        unique_universes[universe_id] = get_openmc_universe(universe)

    # Build the nested Universe array
    for z in range(dimension[2]):
        for y in range(dimension[1]):
            for x in range(dimension[0]):
                universe_id = universes[z][y][x]._id
                universe_array[x][y][z] = unique_universes[universe_id]

    # Reverse y-dimension in array to match ordering in OpenCG
    universe_array = universe_array[:,::-1,:]

    lower_left = np.array(offset, dtype=np.float64) + \
                             ((np.array(width, dtype=np.float64) * \
                               np.array(dimension, dtype=np.float64))) / -2.0

    openmc_lattice = openmc.RectLattice(lattice_id=lattice_id)
    openmc_lattice.dimension = dimension
    openmc_lattice.pitch = width
    openmc_lattice.universes = universe_array
    openmc_lattice.lower_left = lower_left

    # Add the OpenMC Lattice to the global collection of all OpenMC Lattices
    OPENMC_LATTICES[lattice_id] = openmc_lattice

    # Add the OpenCG Lattice to the global collection of all OpenCG Lattices
    OPENCG_LATTICES[lattice_id] = opencg_lattice

    return openmc_lattice


def get_opencg_geometry(openmc_geometry):

    if not isinstance(openmc_geometry, openmc.Geometry):
        msg = 'Unable to get OpenCG geometry from {0} which is ' \
              'not an OpenMC Geometry object'.format(openmc_geometry)
        raise ValueError(msg)

    # Clear dictionaries and auto-generated IDs
    OPENMC_SURFACES.clear()
    OPENCG_SURFACES.clear()
    OPENMC_CELLS.clear()
    OPENCG_CELLS.clear()
    OPENMC_UNIVERSES.clear()
    OPENCG_UNIVERSES.clear()
    OPENMC_LATTICES.clear()
    OPENCG_LATTICES.clear()

    openmc_root_universe = openmc_geometry._root_universe
    opencg_root_universe = get_opencg_universe(openmc_root_universe)

    opencg_geometry = opencg.Geometry()
    opencg_geometry.setRootUniverse(opencg_root_universe)
    opencg_geometry.initializeCellOffsets()

    return opencg_geometry


def get_openmc_geometry(opencg_geometry):

    if not isinstance(opencg_geometry, opencg.Geometry):
        msg = 'Unable to get OpenMC geometry from {0} which is ' \
              'not an OpenCG Geometry object'.format(opencg_geometry)
        raise ValueError(msg)

    # Deep copy the goemetry since it may be modified to make all Surfaces
    # compatible with OpenMC's specifications
    opencg_geometry.assignAutoIds()
    opencg_geometry = copy.deepcopy(opencg_geometry)

    # Update Cell bounding boxes in Geometry
    opencg_geometry.updateBoundingBoxes()

    # Clear dictionaries and auto-generated ID
    OPENMC_SURFACES.clear()
    OPENCG_SURFACES.clear()
    OPENMC_CELLS.clear()
    OPENCG_CELLS.clear()
    OPENMC_UNIVERSES.clear()
    OPENCG_UNIVERSES.clear()
    OPENMC_LATTICES.clear()
    OPENCG_LATTICES.clear()

    # Make the entire geometry "compatible" before assigning auto IDs
    universes = opencg_geometry.getAllUniverses()
    for universe_id, universe in universes.items():
      if not isinstance(universe, opencg.Lattice):
        make_opencg_cells_compatible(universe)

    opencg_geometry.assignAutoIds()

    opencg_root_universe = opencg_geometry._root_universe
    openmc_root_universe = get_openmc_universe(opencg_root_universe)

    openmc_geometry = openmc.Geometry()
    openmc_geometry.root_universe = openmc_root_universe

    return openmc_geometry
