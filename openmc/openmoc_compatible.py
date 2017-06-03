import copy
import operator

import numpy as np

import openmoc
import openmc
import openmc.checkvalue as cv


# A dictionary of all OpenMC Materials created
# Keys   - Material IDs
# Values - Materials
OPENMC_MATERIALS = {}

# A dictionary of all OpenMOC Materials created
# Keys   - Material IDs
# Values - Materials
OPENMOC_MATERIALS = {}

# A dictionary of all OpenMC Surfaces created
# Keys   - Surface IDs
# Values - Surfaces
OPENMC_SURFACES = {}

# A dictionary of all OpenMOC Surfaces created
# Keys   - Surface IDs
# Values - Surfaces
OPENMOC_SURFACES = {}

# A dictionary of all OpenMC Cells created
# Keys   - Cell IDs
# Values - Cells
OPENMC_CELLS = {}

# A dictionary of all OpenMOC Cells created
# Keys   - Cell IDs
# Values - Cells
OPENMOC_CELLS = {}

# A dictionary of all OpenMC Universes created
# Keys   - Universes IDs
# Values - Universes
OPENMC_UNIVERSES = {}

# A dictionary of all OpenMOC Universes created
# Keys   - Universes IDs
# Values - Universes
OPENMOC_UNIVERSES = {}

# A dictionary of all OpenMC Lattices created
# Keys   - Lattice IDs
# Values - Lattices
OPENMC_LATTICES = {}

# A dictionary of all OpenMOC Lattices created
# Keys   - Lattice IDs
# Values - Lattices
OPENMOC_LATTICES = {}


def get_openmoc_material(openmc_material):
    """Return an OpenMOC material corresponding to an OpenMC material.

    Parameters
    ----------
    openmc_material : openmc.Material
        OpenMC material

    Returns
    -------
    openmoc_material : openmoc.Material
        Equivalent OpenMOC material

    """

    cv.check_type('openmc_material', openmc_material, openmc.Material)

    material_id = openmc_material.id

    # If this Material was already created, use it
    if material_id in OPENMOC_MATERIALS:
        return OPENMOC_MATERIALS[material_id]

    # Create an OpenMOC Material to represent this OpenMC Material
    name = str(openmc_material.name)
    openmoc_material = openmoc.Material(id=material_id, name=name)

    # Add the OpenMC Material to the global collection of all OpenMC Materials
    OPENMC_MATERIALS[material_id] = openmc_material

    # Add the OpenMOC Material to the global collection of all OpenMOC Materials
    OPENMOC_MATERIALS[material_id] = openmoc_material

    return openmoc_material


def get_openmc_material(openmoc_material):
    """Return an OpenMC material corresponding to an OpenMOC material.

    Parameters
    ----------
    openmoc_material : openmoc.Material
        OpenMOC material

    Returns
    -------
    openmc_material : openmc.Material
        Equivalent OpenMC material

    """

    cv.check_type('openmoc_material', openmoc_material, openmoc.Material)

    material_id = openmoc_material.getId()

    # If this Material was already created, use it
    if material_id in OPENMC_MATERIALS:
        return OPENMC_MATERIALS[material_id]

    # Create an OpenMC Material to represent this OpenMOC Material
    name = openmoc_material.getName()
    openmc_material = openmc.Material(material_id=material_id, name=name)

    # Add the OpenMOC Material to the global collection of all OpenMOC Materials
    OPENMOC_MATERIALS[material_id] = openmoc_material

    # Add the OpenMC Material to the global collection of all OpenMC Materials
    OPENMC_MATERIALS[material_id] = openmc_material

    return openmc_material


def get_openmoc_surface(openmc_surface):
    """Return an OpenMOC surface corresponding to an OpenMC surface.

    Parameters
    ----------
    openmc_surface : openmc.Surface
        OpenMC surface

    Returns
    -------
    openmoc_surface : openmoc.Surface
        Equivalent OpenMOC surface

    """

    cv.check_type('openmc_surface', openmc_surface, openmc.Surface)

    surface_id = openmc_surface.id

    # If this Material was already created, use it
    if surface_id in OPENMOC_SURFACES:
        return OPENMOC_SURFACES[surface_id]

    # Create an OpenMOC Surface to represent this OpenMC Surface
    name = openmc_surface.name

    # Determine the type of boundary conditions applied to the Surface
    if openmc_surface.boundary_type == 'vacuum':
        boundary = openmoc.VACUUM
    elif openmc_surface.boundary_type == 'reflective':
        boundary = openmoc.REFLECTIVE
    elif openmc_surface.boundary_type == 'periodic':
        boundary = openmoc.PERIODIC
    else:
        boundary = openmoc.BOUNDARY_NONE

    if openmc_surface.type == 'plane':
        A = openmc_surface.a
        B = openmc_surface.b
        C = openmc_surface.c
        D = openmc_surface.d

        # OpenMOC uses the opposite sign on D
        openmoc_surface = openmoc.Plane(A, B, C, -D, surface_id, name)

    elif openmc_surface.type == 'x-plane':
        x0 = openmc_surface.x0
        openmoc_surface = openmoc.XPlane(x0, surface_id, name)

    elif openmc_surface.type == 'y-plane':
        y0 = openmc_surface.y0
        openmoc_surface = openmoc.YPlane(y0, surface_id, name)

    elif openmc_surface.type == 'z-plane':
        z0 = openmc_surface.z0
        openmoc_surface = openmoc.ZPlane(z0, surface_id, name)

    elif openmc_surface.type == 'z-cylinder':
        x0 = openmc_surface.x0
        y0 = openmc_surface.y0
        R = openmc_surface.r
        openmoc_surface = openmoc.ZCylinder(x0, y0, R, surface_id, name)

    else:
        msg = 'Unable to create an OpenMOC Surface from an OpenMC ' \
              'Surface of type "{}" since it is not a compatible ' \
              'Surface type in OpenMOC'.format(type(openmc_surface))
        raise ValueError(msg)

    # Set the boundary condition for this Surface
    openmoc_surface.setBoundaryType(boundary)

    # Add the OpenMC Surface to the global collection of all OpenMC Surfaces
    OPENMC_SURFACES[surface_id] = openmc_surface

    # Add the OpenMOC Surface to the global collection of all OpenMOC Surfaces
    OPENMOC_SURFACES[surface_id] = openmoc_surface

    return openmoc_surface


def get_openmc_surface(openmoc_surface):
    """Return an OpenMC surface corresponding to an OpenMOC surface.

    Parameters
    ----------
    openmoc_surface : openmoc.Surface
        OpenMOC surface

    Returns
    -------
    openmc_surface : openmc.Surface
        Equivalent OpenMC surface

    """

    cv.check_type('openmoc_surface', openmoc_surface, openmoc.Surface)

    surface_id = openmoc_surface.getId()

    # If this Surface was already created, use it
    if surface_id in OPENMC_SURFACES:
        return OPENMC_SURFACES[surface_id]

    # Create an OpenMC Surface to represent this OpenMOC Surface
    name = openmoc_surface.getName()

    # Correct for OpenMC's syntax for Surfaces dividing Cells
    boundary = openmoc_surface.getBoundaryType()
    if boundary == openmoc.VACUUM:
        boundary = 'vacuum'
    elif boundary == openmoc.REFLECTIVE:
        boundary = 'reflective'
    elif boundary == openmoc.PERIODIC:
        boundary = 'periodic'
    else:
        boundary = 'transmission'

    if openmoc_surface.getSurfaceType() == openmoc.PLANE:
        openmoc_surface = openmoc.castSurfaceToPlane(openmoc_surface)
        A = openmoc_surface.getA()
        B = openmoc_surface.getB()
        C = openmoc_surface.getC()
        D = openmoc_surface.getD()

        # OpenMOC uses the opposite sign on D
        openmc_surface = openmc.Plane(surface_id, boundary, A, B, C, -D, name)

    elif openmoc_surface.getSurfaceType() == openmoc.XPLANE:
        openmoc_surface = openmoc.castSurfaceToXPlane(openmoc_surface)
        x0 = openmoc_surface.getX()
        openmc_surface = openmc.XPlane(surface_id, boundary, x0, name)

    elif openmoc_surface.getSurfaceType() == openmoc.YPLANE:
        openmoc_surface = openmoc.castSurfaceToYPlane(openmoc_surface)
        y0 = openmoc_surface.getY()
        openmc_surface = openmc.YPlane(surface_id, boundary, y0, name)

    elif openmoc_surface.getSurfaceType() == openmoc.ZPLANE:
        openmoc_surface = openmoc.castSurfaceToZPlane(openmoc_surface)
        z0 = openmoc_surface.getZ()
        openmc_surface = openmc.ZPlane(surface_id, boundary, z0, name)

    elif openmoc_surface.getSurfaceType() == openmoc.ZCYLINDER:
        openmoc_surface = openmoc.castSurfaceToZCylinder(openmoc_surface)
        x0 = openmoc_surface.getX0()
        y0 = openmoc_surface.getY0()
        R = openmoc_surface.getRadius()
        openmc_surface = openmc.ZCylinder(surface_id, boundary, x0, y0, R, name)

    # Add the OpenMC Surface to the global collection of all OpenMC Surfaces
    OPENMC_SURFACES[surface_id] = openmc_surface

    # Add the OpenMOC Surface to the global collection of all OpenMOC Surfaces
    OPENMOC_SURFACES[surface_id] = openmoc_surface

    return openmc_surface


def get_openmoc_cell(openmc_cell):
    """Return an OpenMOC cell corresponding to an OpenMC cell.

    Parameters
    ----------
    openmc_cell : openmc.Cell
        OpenMC cell

    Returns
    -------
    openmoc_cell : openmoc.Cell
        Equivalent OpenMOC cell

    """

    cv.check_type('openmc_cell', openmc_cell, openmc.Cell)

    cell_id = openmc_cell.id

    # If this Cell was already created, use it
    if cell_id in OPENMOC_CELLS:
        return OPENMOC_CELLS[cell_id]

    # Create an OpenMOC Cell to represent this OpenMC Cell
    name = openmc_cell.name
    openmoc_cell = openmoc.Cell(cell_id, name)

    fill = openmc_cell.fill

    if openmc_cell.fill_type == 'material':
        openmoc_cell.setFill(get_openmoc_material(fill))
    elif openmc_cell.fill_type == 'universe':
        openmoc_cell.setFill(get_openmoc_universe(fill))
    else:
        openmoc_cell.setFill(get_openmoc_lattice(fill))

    if openmc_cell.rotation is not None:
        rotation = np.asarray(openmc_cell.rotation, dtype=np.float64)
        openmoc_cell.setRotation(rotation)
    if openmc_cell.translation is not None:
        translation = np.asarray(openmc_cell.translation, dtype=np.float64)
        openmoc_cell.setTranslation(translation)

    # Convert OpenMC's cell region to an equivalent OpenMOC region
    if openmc_cell.region is not None:
        openmoc_cell.setRegion(get_openmoc_region(openmc_cell.region))

    # Add the OpenMC Cell to the global collection of all OpenMC Cells
    OPENMC_CELLS[cell_id] = openmc_cell

    # Add the OpenMOC Cell to the global collection of all OpenMOC Cells
    OPENMOC_CELLS[cell_id] = openmoc_cell

    return openmoc_cell


def get_openmoc_region(openmc_region):
    """Return an OpenMOC region corresponding to an OpenMC region.

    Parameters
    ----------
    openmc_region : openmc.Region
        OpenMC region

    Returns
    -------
    openmoc_region : openmoc.Region
        Equivalent OpenMOC region

    """

    cv.check_type('openmc_region', openmc_region, openmc.Region)

    # Recursively instantiate a region of the appropriate type
    if isinstance(openmc_region, openmc.Halfspace):
        surface = openmc_region.surface
        halfspace = -1 if openmc_region.side == '-' else 1
        openmoc_region = \
            openmoc.Halfspace(halfspace, get_openmoc_surface(surface))
    elif isinstance(openmc_region, openmc.Intersection):
        openmoc_region = openmoc.Intersection()
        for openmc_node in openmc_region:
            openmoc_region.addNode(get_openmoc_region(openmc_node))
    elif isinstance(openmc_region, openmc.Union):
        openmoc_region = openmoc.Union()
        for openmc_node in openmc_region:
            openmoc_region.addNode(get_openmoc_region(openmc_node))
    elif isinstance(openmc_region, openmc.Complement):
        openmoc_region = openmoc.Complement()
        openmoc_region.addNode(get_openmoc_region(openmc_region.node))

    return openmoc_region


def get_openmc_region(openmoc_region):
    """Return an OpenMC region corresponding to an OpenMOC region.

    Parameters
    ----------
    openmoc_region : openmoc.Region
        OpenMOC region

    Returns
    -------
    openmc_region : openmc.Region
        Equivalent OpenMC region

    """

    cv.check_type('openmoc_region', openmoc_region, openmoc.Region)

    # Recursively instantiate a region of the appropriate type
    if openmoc_region.getRegionType() == openmoc.HALFSPACE:
        openmoc_region = openmoc.castRegionToHalfspace(openmoc_region)
        surface = get_openmc_surface(openmoc_region.getSurface())
        side = '-' if openmoc_region.getHalfspace() == -1 else '+'
        openmc_region = openmc.Halfspace(surface, side)
    elif openmoc_region.getRegionType() == openmoc.INTERSECTION:
        openmc_region = openmc.Intersection()
        for openmoc_node in openmoc_region.getNodes():
            openmc_node = get_openmc_region(openmoc_node)
            openmc_region.append(openmc_node)
    elif openmoc_region.getRegionType() == openmoc.UNION:
        openmc_region = openmc.Union()
        for openmoc_node in openmoc_region.getNodes():
            openmc_node = get_openmc_region(openmoc_node)
            openmc_region.append(openmc_node)
    elif openmoc_region.getRegionType() == openmoc.COMPLEMENT:
        openmoc_nodes = openmoc_region.getNodes()
        openmc_node = get_openmc_region(openmoc_nodes[0])
        openmc_region = openmc.Complement(openmc_node)

    return openmc_region


def get_openmc_cell(openmoc_cell):
    """Return an OpenMC cell corresponding to an OpenMOC cell.

    Parameters
    ----------
    openmoc_cell : openmoc.Cell
        OpenMOC cell

    Returns
    -------
    openmc_cell : openmc.Cell
        Equivalent OpenMC cell

    """

    cv.check_type('openmoc_cell', openmoc_cell, openmoc.Cell)

    cell_id = openmoc_cell.getId()

    # If this Cell was already created, use it
    if cell_id in OPENMC_CELLS:
        return OPENMC_CELLS[cell_id]

    # Create an OpenMOC Cell to represent this OpenMC Cell
    name = openmoc_cell.getName()
    openmc_cell = openmc.Cell(cell_id, name)

    if (openmoc_cell.getType() == openmoc.MATERIAL):
        fill = openmoc_cell.getFillMaterial()
        openmc_cell.fill = get_openmc_material(fill)
    elif (openmoc_cell.getType() == openmoc.FILL):
        fill = openmoc_cell.getFillUniverse()
        if fill.getType() == openmoc.LATTICE:
            fill = openmoc.castUniverseToLattice(fill)
            openmc_cell.fill = get_openmc_lattice(fill)
        else:
            openmc_cell.fill = get_openmc_universe(fill)

    if openmoc_cell.isRotated():
        rotation = openmoc_cell.retrieveRotation(3)
        openmc_cell.rotation = rotation
    if openmoc_cell.isTranslated():
        translation = openmoc_cell.retrieveTranslation(3)
        openmc_cell.translation = translation

    # Convert OpenMC's cell region to an equivalent OpenMOC region
    openmoc_region = openmoc_cell.getRegion()
    if openmoc_region is not None:
        openmc_cell.region = get_openmc_region(openmoc_region)

    # Add the OpenMC Cell to the global collection of all OpenMC Cells
    OPENMC_CELLS[cell_id] = openmc_cell

    # Add the OpenMOC Cell to the global collection of all OpenMOC Cells
    OPENMOC_CELLS[cell_id] = openmoc_cell

    return openmc_cell


def get_openmoc_universe(openmc_universe):
    """Return an OpenMOC universe corresponding to an OpenMC universe.

    Parameters
    ----------
    openmc_universe : openmc.Universe
        OpenMC universe

    Returns
    -------
    openmoc_universe : openmoc.Universe
        Equivalent OpenMOC universe

    """

    cv.check_type('openmc_universe', openmc_universe, openmc.Universe)

    universe_id = openmc_universe.id

    # If this Universe was already created, use it
    if universe_id in OPENMOC_UNIVERSES:
        return OPENMOC_UNIVERSES[universe_id]

    # Create an OpenMOC Universe to represent this OpenMC Universe
    name = openmc_universe.name
    openmoc_universe = openmoc.Universe(universe_id, name)

    # Convert all OpenMC Cells in this Universe to OpenMOC Cells
    openmc_cells = openmc_universe.cells

    for openmc_cell in openmc_cells.values():
        openmoc_cell = get_openmoc_cell(openmc_cell)
        openmoc_universe.addCell(openmoc_cell)

    # Add the OpenMC Universe to the global collection of all OpenMC Universes
    OPENMC_UNIVERSES[universe_id] = openmc_universe

    # Add the OpenMOC Universe to the global collection of all OpenMOC Universes
    OPENMOC_UNIVERSES[universe_id] = openmoc_universe

    return openmoc_universe


def get_openmc_universe(openmoc_universe):
    """Return an OpenMC universe corresponding to an OpenMOC universe.

    Parameters
    ----------
    openmoc_universe : openmoc.Universe
        OpenMOC universe

    Returns
    -------
    openmc_universe : openmc.Universe
        Equivalent OpenMC universe

    """

    cv.check_type('openmoc_universe', openmoc_universe, openmoc.Universe)

    universe_id = openmoc_universe.getId()

    # If this Universe was already created, use it
    if universe_id in OPENMC_UNIVERSES:
        return OPENMC_UNIVERSES[universe_id]

    # Create an OpenMC Universe to represent this OpenMOC Universe
    name = openmoc_universe.getName()
    openmc_universe = openmc.Universe(universe_id, name)

    # Convert all OpenMOC Cells in this Universe to OpenMC Cells
    for openmoc_cell in openmoc_universe.getCells().values():
        openmc_cell = get_openmc_cell(openmoc_cell)
        openmc_universe.add_cell(openmc_cell)

    # Add the OpenMC Universe to the global collection of all OpenMC Universes
    OPENMC_UNIVERSES[universe_id] = openmc_universe

    # Add the OpenMOC Universe to the global collection of all OpenMOC Universes
    OPENMOC_UNIVERSES[universe_id] = openmoc_universe

    return openmc_universe


def get_openmoc_lattice(openmc_lattice):
    """Return an OpenMOC lattice corresponding to an OpenMOC lattice.

    Parameters
    ----------
    openmc_lattice : openmc.RectLattice
        OpenMC lattice

    Returns
    -------
    openmoc_lattice : openmoc.Lattice
        Equivalent OpenMOC lattice

    """

    cv.check_type('openmc_lattice', openmc_lattice, openmc.RectLattice)

    lattice_id = openmc_lattice.id

    # If this Lattice was already created, use it
    if lattice_id in OPENMOC_LATTICES:
        return OPENMOC_LATTICES[lattice_id]

    # Create an OpenMOC Lattice to represent this OpenMC Lattice
    name = openmc_lattice.name
    dimension = openmc_lattice.shape
    pitch = openmc_lattice.pitch
    lower_left = openmc_lattice.lower_left
    universes = openmc_lattice.universes

    # Convert 2D dimension to 3D for OpenMOC
    if len(dimension) == 2:
        new_dimension = np.ones(3, dtype=int)
        new_dimension[:2] = dimension
        dimension = new_dimension

    # Convert 2D pitch to 3D for OpenMOC
    if len(pitch) == 2:
        new_pitch = np.ones(3, dtype=np.float64) * np.finfo(np.float64).max
        new_pitch[:2] = pitch
        pitch = new_pitch

    # Convert 2D lower left to 3D for OpenMOC
    if len(lower_left) == 2:
        new_lower_left = np.ones(3, dtype=np.float64)
        new_lower_left *= np.finfo(np.float64).min / 2.
        new_lower_left[:2] = lower_left
        lower_left = new_lower_left

        # Convert 2D universes array to 3D for OpenMOC
        if len(universes.shape) == 2:
            new_universes = universes.copy()
            new_universes.shape = (1,) + universes.shape
            universes = new_universes

    # Initialize an empty array for the OpenMOC nested Universes in this Lattice
    universe_array = np.ndarray(tuple(dimension[::-1]), dtype=openmoc.Universe)

    # Create OpenMOC Universes for each unique nested Universe in this Lattice
    unique_universes = openmc_lattice.get_unique_universes()

    for universe_id, universe in unique_universes.items():
        unique_universes[universe_id] = get_openmoc_universe(universe)

    # Build the nested Universe array
    for z in range(dimension[2]):
        for y in range(dimension[1]):
            for x in range(dimension[0]):
                universe_id = universes[z][y][x].id
                universe_array[z][y][x] = unique_universes[universe_id]

    openmoc_lattice = openmoc.Lattice(lattice_id, name)
    openmoc_lattice.setWidth(pitch[0], pitch[1], pitch[2])
    openmoc_lattice.setUniverses(universe_array.tolist())

    offset = np.array(lower_left, dtype=np.float64) - \
             ((np.array(pitch, dtype=np.float64) *
               np.array(dimension, dtype=np.float64))) / -2.0
    openmoc_lattice.setOffset(offset[0], offset[1], offset[2])

    # Add the OpenMC Lattice to the global collection of all OpenMC Lattices
    OPENMC_LATTICES[lattice_id] = openmc_lattice

    # Add the OpenMOC Lattice to the global collection of all OpenMOC Lattices
    OPENMOC_LATTICES[lattice_id] = openmoc_lattice

    return openmoc_lattice


def get_openmc_lattice(openmoc_lattice):
    """Return an OpenMC lattice corresponding to an OpenMOC lattice.

    Parameters
    ----------
    openmoc_lattice : openmoc.Lattice
        OpenMOC lattice

    Returns
    -------
    openmc_lattice : openmc.RectLattice
        Equivalent OpenMC lattice

    """

    cv.check_type('openmoc_lattice', openmoc_lattice, openmoc.Lattice)

    lattice_id = openmoc_lattice.getId()

    # If this Lattice was already created, use it
    if lattice_id in OPENMC_LATTICES:
        return OPENMC_LATTICES[lattice_id]

    name = openmoc_lattice.getName()
    dimension = [openmoc_lattice.getNumX(),
                 openmoc_lattice.getNumY(),
                 openmoc_lattice.getNumZ()]
    width = [openmoc_lattice.getWidthX(),
             openmoc_lattice.getWidthY(),
             openmoc_lattice.getWidthZ()]
    offset = openmoc_lattice.getOffset()
    offset = [offset.getX(), offset.getY(), offset.getZ()]
    lower_left = np.array(offset, dtype=np.float64) + \
                 ((np.array(width, dtype=np.float64) *
                   np.array(dimension, dtype=np.float64))) / -2.0

    # Initialize an empty array for the OpenMOC nested Universes in this Lattice
    universe_array = np.ndarray(tuple(np.array(dimension)),
                                dtype=openmoc.Universe)

    # Create OpenMOC Universes for each unique nested Universe in this Lattice
    unique_universes = openmoc_lattice.getUniqueUniverses()

    for universe_id, universe in unique_universes.items():
        unique_universes[universe_id] = get_openmc_universe(universe)

    # Build the nested Universe array
    for x in range(dimension[0]):
        for y in range(dimension[1]):
            for z in range(dimension[2]):
                universe = openmoc_lattice.getUniverse(x, y, z)
                universe_id = universe.getId()
                universe_array[x][y][z] = \
                    unique_universes[universe_id]

    universe_array = np.swapaxes(universe_array, 0, 2)

    # Convert axially infinite 3D OpenMOC lattice to a 2D OpenMC lattice
    if width[2] == np.finfo(np.float64).max:
        dimension = dimension[:2]
        width = width[:2]
        offset = offset[:2]
        lower_left = lower_left[:2]
        universe_array = np.squeeze(universe_array, 2)

    openmc_lattice = openmc.RectLattice(lattice_id=lattice_id, name=name)
    openmc_lattice.pitch = width
    openmc_lattice.lower_left = lower_left
    openmc_lattice.universes = universe_array

    # Add the OpenMC Lattice to the global collection of all OpenMC Lattices
    OPENMC_LATTICES[lattice_id] = openmc_lattice

    # Add the OpenMOC Lattice to the global collection of all OpenMOC Lattices
    OPENMOC_LATTICES[lattice_id] = openmoc_lattice

    return openmc_lattice


def get_openmoc_geometry(openmc_geometry):
    """Return an OpenMC geometry corresponding to an OpenMOC geometry.

    Parameters
    ----------
    openmc_geometry : openmc.Geometry
        OpenMC geometry

    Returns
    -------
    openmoc_geometry : openmoc.Geometry
        Equivalent OpenMOC geometry

    """

    cv.check_type('openmc_geometry', openmc_geometry, openmc.Geometry)

    # Clear dictionaries and auto-generated IDs
    OPENMC_SURFACES.clear()
    OPENMOC_SURFACES.clear()
    OPENMC_CELLS.clear()
    OPENMOC_CELLS.clear()
    OPENMC_UNIVERSES.clear()
    OPENMOC_UNIVERSES.clear()
    OPENMC_LATTICES.clear()
    OPENMOC_LATTICES.clear()

    openmc_root_universe = openmc_geometry.root_universe
    openmoc_root_universe = get_openmoc_universe(openmc_root_universe)

    openmoc_geometry = openmoc.Geometry()
    openmoc_geometry.setRootUniverse(openmoc_root_universe)

    # Update OpenMOC's auto-generated object IDs (e.g., Surface, Material)
    # with the maximum of those created from the OpenMC objects
    all_materials = openmoc_geometry.getAllMaterials()
    all_surfaces = openmoc_geometry.getAllSurfaces()
    all_cells = openmoc_geometry.getAllCells()
    all_universes = openmoc_geometry.getAllUniverses()

    max_material_id = max(all_materials.keys())
    max_surface_id = max(all_surfaces.keys())
    max_cell_id = max(all_cells.keys())
    max_universe_id = max(all_universes.keys())

    openmoc.maximize_material_id(max_material_id+1)
    openmoc.maximize_surface_id(max_surface_id+1)
    openmoc.maximize_cell_id(max_cell_id+1)
    openmoc.maximize_universe_id(max_universe_id+1)

    return openmoc_geometry


def get_openmc_geometry(openmoc_geometry):
    """Return an OpenMC geometry corresponding to an OpenMOC geometry.

    Parameters
    ----------
    openmoc_geometry : openmoc.Geometry
        OpenMOC geometry

    Returns
    -------
    openmc_geometry : openmc.Geometry
        Equivalent OpenMC geometry

    """

    cv.check_type('openmoc_geometry', openmoc_geometry, openmoc.Geometry)

    # Clear dictionaries and auto-generated ID
    OPENMC_SURFACES.clear()
    OPENMOC_SURFACES.clear()
    OPENMC_CELLS.clear()
    OPENMOC_CELLS.clear()
    OPENMC_UNIVERSES.clear()
    OPENMOC_UNIVERSES.clear()
    OPENMC_LATTICES.clear()
    OPENMOC_LATTICES.clear()

    openmoc_root_universe = openmoc_geometry.getRootUniverse()
    openmc_root_universe = get_openmc_universe(openmoc_root_universe)

    openmc_geometry = openmc.Geometry()
    openmc_geometry.root_universe = openmc_root_universe

    return openmc_geometry
