import copy
import operator

import numpy as np

try:
    import opencg
except ImportError:
    msg = 'Unable to import opencg which is needed by openmc.opencg_compatible'
    raise ImportError(msg)

import openmc
import openmc.checkvalue as cv


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
    """Return an OpenCG material corresponding to an OpenMC material.

    Parameters
    ----------
    openmc_material : openmc.material.Material
        OpenMC material

    Returns
    -------
    opencg_material : opencg.Material
        Equivalent OpenCG material

    """

    cv.check_type('openmc_material', openmc_material, openmc.Material)

    global OPENCG_MATERIALS
    material_id = openmc_material.id

    # If this Material was already created, use it
    if material_id in OPENCG_MATERIALS:
        return OPENCG_MATERIALS[material_id]

    # Create an OpenCG Material to represent this OpenMC Material
    name = openmc_material.name
    opencg_material = opencg.Material(material_id=material_id, name=name)

    # Add the OpenMC Material to the global collection of all OpenMC Materials
    OPENMC_MATERIALS[material_id] = openmc_material

    # Add the OpenCG Material to the global collection of all OpenCG Materials
    OPENCG_MATERIALS[material_id] = opencg_material

    return opencg_material


def get_openmc_material(opencg_material):
    """Return an OpenMC material corresponding to an OpenCG material.

    Parameters
    ----------
    opencg_material : opencg.Material
        OpenCG material

    Returns
    -------
    openmc_material : openmc.material.Material
        Equivalent OpenMC material

    """

    cv.check_type('opencg_material', opencg_material, opencg.Material)

    global OPENMC_MATERIALS
    material_id = opencg_material.id

    # If this Material was already created, use it
    if material_id in OPENMC_MATERIALS:
        return OPENMC_MATERIALS[material_id]

    # Create an OpenMC Material to represent this OpenCG Material
    name = opencg_material.name
    openmc_material = openmc.Material(material_id=material_id, name=name)

    # Add the OpenMC Material to the global collection of all OpenMC Materials
    OPENMC_MATERIALS[material_id] = openmc_material

    # Add the OpenCG Material to the global collection of all OpenCG Materials
    OPENCG_MATERIALS[material_id] = opencg_material

    return openmc_material


def is_opencg_surface_compatible(opencg_surface):
    """Determine whether OpenCG surface is compatible with OpenMC geometry.

    A surface is considered compatible if there is a one-to-one correspondence
    between OpenMC and OpenCG surface types. Note that some OpenCG surfaces,
    e.g. SquarePrism, do not have a one-to-one correspondence with OpenMC
    surfaces but can still be converted into an equivalent collection of OpenMC
    surfaces.

    Parameters
    ----------
    opencg_surface : opencg.Surface
        OpenCG surface

    Returns
    -------
    bool
        Whether OpenCG surface is compatible with OpenMC

    """

    cv.check_type('opencg_surface', opencg_surface, opencg.Surface)

    if opencg_surface.type in ['x-squareprism',
                               'y-squareprism', 'z-squareprism']:
        return False
    else:
        return True


def get_opencg_surface(openmc_surface):
    """Return an OpenCG surface corresponding to an OpenMC surface.

    Parameters
    ----------
    openmc_surface : openmc.surface.Surface
        OpenMC surface

    Returns
    -------
    opencg_surface : opencg.Surface
        Equivalent OpenCG surface

    """

    cv.check_type('openmc_surface', openmc_surface, openmc.Surface)

    global OPENCG_SURFACES
    surface_id = openmc_surface.id

    # If this Material was already created, use it
    if surface_id in OPENCG_SURFACES:
        return OPENCG_SURFACES[surface_id]

    # Create an OpenCG Surface to represent this OpenMC Surface
    name = openmc_surface.name

    # Correct for OpenMC's syntax for Surfaces dividing Cells
    boundary = openmc_surface.boundary_type
    if boundary == 'transmission':
        boundary = 'interface'

    opencg_surface = None

    if openmc_surface.type == 'plane':
        A = openmc_surface.a
        B = openmc_surface.b
        C = openmc_surface.c
        D = openmc_surface.d
        opencg_surface = opencg.Plane(surface_id, name, boundary, A, B, C, D)

    elif openmc_surface.type == 'x-plane':
        x0 = openmc_surface.x0
        opencg_surface = opencg.XPlane(surface_id, name, boundary, x0)

    elif openmc_surface.type == 'y-plane':
        y0 = openmc_surface.y0
        opencg_surface = opencg.YPlane(surface_id, name, boundary, y0)

    elif openmc_surface.type == 'z-plane':
        z0 = openmc_surface.z0
        opencg_surface = opencg.ZPlane(surface_id, name, boundary, z0)

    elif openmc_surface.type == 'x-cylinder':
        y0 = openmc_surface.y0
        z0 = openmc_surface.z0
        R = openmc_surface.r
        opencg_surface = opencg.XCylinder(surface_id, name,
                                            boundary, y0, z0, R)

    elif openmc_surface.type == 'y-cylinder':
        x0 = openmc_surface.x0
        z0 = openmc_surface.z0
        R = openmc_surface.r
        opencg_surface = opencg.YCylinder(surface_id, name,
                                            boundary, x0, z0, R)

    elif openmc_surface.type == 'z-cylinder':
        x0 = openmc_surface.x0
        y0 = openmc_surface.y0
        R = openmc_surface.r
        opencg_surface = opencg.ZCylinder(surface_id, name,
                                            boundary, x0, y0, R)

    # Add the OpenMC Surface to the global collection of all OpenMC Surfaces
    OPENMC_SURFACES[surface_id] = openmc_surface

    # Add the OpenCG Surface to the global collection of all OpenCG Surfaces
    OPENCG_SURFACES[surface_id] = opencg_surface

    return opencg_surface


def get_openmc_surface(opencg_surface):
    """Return an OpenMC surface corresponding to an OpenCG surface.

    Parameters
    ----------
    opencg_surface : opencg.Surface
        OpenCG surface

    Returns
    -------
    openmc_surface : openmc.surface.Surface
        Equivalent OpenMC surface

    """

    cv.check_type('opencg_surface', opencg_surface, opencg.Surface)

    global openmc_surface
    surface_id = opencg_surface.id

    # If this Surface was already created, use it
    if surface_id in OPENMC_SURFACES:
        return OPENMC_SURFACES[surface_id]

    # Create an OpenMC Surface to represent this OpenCG Surface
    name = opencg_surface.name

    # Correct for OpenMC's syntax for Surfaces dividing Cells
    boundary = opencg_surface.boundary_type
    if boundary == 'interface':
        boundary = 'transmission'

    if opencg_surface.type == 'plane':
        A = opencg_surface.a
        B = opencg_surface.b
        C = opencg_surface.c
        D = opencg_surface.d
        openmc_surface = openmc.Plane(surface_id, boundary, A, B, C, D, name)

    elif opencg_surface.type == 'x-plane':
        x0 = opencg_surface.x0
        openmc_surface = openmc.XPlane(surface_id, boundary, x0, name)

    elif opencg_surface.type == 'y-plane':
        y0 = opencg_surface.y0
        openmc_surface = openmc.YPlane(surface_id, boundary, y0, name)

    elif opencg_surface.type == 'z-plane':
        z0 = opencg_surface.z0
        openmc_surface = openmc.ZPlane(surface_id, boundary, z0, name)

    elif opencg_surface.type == 'x-cylinder':
        y0 = opencg_surface.y0
        z0 = opencg_surface.z0
        R = opencg_surface.r
        openmc_surface = openmc.XCylinder(surface_id, boundary, y0, z0, R, name)

    elif opencg_surface.type == 'y-cylinder':
        x0 = opencg_surface.x0
        z0 = opencg_surface.z0
        R = opencg_surface.r
        openmc_surface = openmc.YCylinder(surface_id, boundary, x0, z0, R, name)

    elif opencg_surface.type == 'z-cylinder':
        x0 = opencg_surface.x0
        y0 = opencg_surface.y0
        R = opencg_surface.r
        openmc_surface = openmc.ZCylinder(surface_id, boundary, x0, y0, R, name)

    else:
        msg = 'Unable to create an OpenMC Surface from an OpenCG ' \
              'Surface of type "{0}" since it is not a compatible ' \
              'Surface type in OpenMC'.format(opencg_surface.type)
        raise ValueError(msg)

    # Add the OpenMC Surface to the global collection of all OpenMC Surfaces
    OPENMC_SURFACES[surface_id] = openmc_surface

    # Add the OpenCG Surface to the global collection of all OpenCG Surfaces
    OPENCG_SURFACES[surface_id] = opencg_surface

    return openmc_surface


def get_compatible_opencg_surfaces(opencg_surface):
    """Generate OpenCG surfaces that are compatible with OpenMC equivalent to an
    OpenCG surface that is not compatible. For example, this method may be used
    to convert a ZSquarePrism OpenCG surface into a collection of equivalent
    XPlane and YPlane OpenCG surfaces.

    Parameters
    ----------
    opencg_surface : opencg.Surface
        OpenCG surface that is incompatible with OpenMC

    Returns
    -------
    surfaces : list of opencg.Surface
        Collection of surfaces equivalent to the original one but compatible
        with OpenMC

    """

    cv.check_type('opencg_surface', opencg_surface, opencg.Surface)

    global OPENMC_SURFACES
    surface_id = opencg_surface.id

    # If this Surface was already created, use it
    if surface_id in OPENMC_SURFACES:
        return OPENMC_SURFACES[surface_id]

    # Create an OpenMC Surface to represent this OpenCG Surface
    name = opencg_surface.name
    boundary = opencg_surface.boundary_type

    if opencg_surface.type == 'x-squareprism':
        y0 = opencg_surface.y0
        z0 = opencg_surface.z0
        R = opencg_surface.r

        # Create a list of the four planes we need
        left = opencg.YPlane(name=name, boundary=boundary, y0=y0-R)
        right = opencg.YPlane(name=name, boundary=boundary, y0=y0+R)
        bottom = opencg.ZPlane(name=name, boundary=boundary, z0=z0-R)
        top = opencg.ZPlane(name=name, boundary=boundary, z0=z0+R)
        surfaces = [left, right, bottom, top]

    elif opencg_surface.type == 'y-squareprism':
        x0 = opencg_surface.x0
        z0 = opencg_surface.z0
        R = opencg_surface.r

        # Create a list of the four planes we need
        left = opencg.XPlane(name=name, boundary=boundary, x0=x0-R)
        right = opencg.XPlane(name=name, boundary=boundary, x0=x0+R)
        bottom = opencg.ZPlane(name=name, boundary=boundary, z0=z0-R)
        top = opencg.ZPlane(name=name, boundary=boundary, z0=z0+R)
        surfaces = [left, right, bottom, top]

    elif opencg_surface.type == 'z-squareprism':
        x0 = opencg_surface.x0
        y0 = opencg_surface.y0
        R = opencg_surface.r

        # Create a list of the four planes we need
        left = opencg.XPlane(name=name, boundary=boundary, x0=x0-R)
        right = opencg.XPlane(name=name, boundary=boundary, x0=x0+R)
        bottom = opencg.YPlane(name=name, boundary=boundary, y0=y0-R)
        top = opencg.YPlane(name=name, boundary=boundary, y0=y0+R)
        surfaces = [left, right, bottom, top]

    else:
        msg = 'Unable to create a compatible OpenMC Surface an OpenCG ' \
              'Surface of type "{0}" since it already a compatible ' \
              'Surface type in OpenMC'.format(opencg_surface.type)
        raise ValueError(msg)

    # Add the OpenMC Surface(s) to the global collection of all OpenMC Surfaces
    OPENMC_SURFACES[surface_id] = surfaces

    # Add the OpenCG Surface to the global collection of all OpenCG Surfaces
    OPENCG_SURFACES[surface_id] = opencg_surface

    return surfaces


def get_opencg_cell(openmc_cell):
    """Return an OpenCG cell corresponding to an OpenMC cell.

    Parameters
    ----------
    openmc_cell : openmc.universe.Cell
        OpenMC cell

    Returns
    -------
    opencg_cell : opencg.Cell
        Equivalent OpenCG cell

    """

    cv.check_type('openmc_cell', openmc_cell, openmc.Cell)

    global OPENCG_CELLS
    cell_id = openmc_cell.id

    # If this Cell was already created, use it
    if cell_id in OPENCG_CELLS:
        return OPENCG_CELLS[cell_id]

    # Create an OpenCG Cell to represent this OpenMC Cell
    name = openmc_cell.name
    opencg_cell = opencg.Cell(cell_id, name)

    fill = openmc_cell.fill

    if openmc_cell.fill_type == 'material':
        opencg_cell.fill = get_opencg_material(fill)
    elif openmc_cell.fill_type == 'universe':
        opencg_cell.fill = get_opencg_universe(fill)
    else:
        opencg_cell.fill = get_opencg_lattice(fill)

    if openmc_cell.rotation is not None:
        opencg_cell.rotation = openmc_cell.rotation

    if openmc_cell.translation is not None:
        opencg_cell.translation = openmc_cell.translation

    # Add surfaces to OpenCG cell from OpenMC cell region. Right now this only
    # works if the region is a single half-space or an intersection of
    # half-spaces, i.e., no complex cells.
    region = openmc_cell.region
    if region is not None:
        if isinstance(region, openmc.Halfspace):
            surface = region.surface
            halfspace = -1 if region.side == '-' else 1
            opencg_cell.add_surface(get_opencg_surface(surface), halfspace)
        elif isinstance(region, openmc.Intersection):
            for node in region.nodes:
                if not isinstance(node, openmc.Halfspace):
                    raise NotImplementedError("Complex cells not yet "
                                              "supported in OpenCG.")
                surface = node.surface
                halfspace = -1 if node.side == '-' else 1
                opencg_cell.add_surface(get_opencg_surface(surface), halfspace)
        else:
            raise NotImplementedError("Complex cells not yet supported "
                                       "in OpenCG.")

    # Add the OpenMC Cell to the global collection of all OpenMC Cells
    OPENMC_CELLS[cell_id] = openmc_cell

    # Add the OpenCG Cell to the global collection of all OpenCG Cells
    OPENCG_CELLS[cell_id] = opencg_cell

    return opencg_cell


def get_compatible_opencg_cells(opencg_cell, opencg_surface, halfspace):
    """Generate OpenCG cells that are compatible with OpenMC equivalent to an OpenCG
    cell that is not compatible.

    Parameters
    ----------
    opencg_cell : opencg.Cell
        OpenCG cell
    opencg_surface : opencg.Surface
        OpenCG surface that causes the incompatibility, e.g. an instance of
        XSquarePrism
    halfspace : {-1, 1}
        Which halfspace defined by the surface is contained in the cell

    Returns
    -------
    compatible_cells : list of opencg.Cell
        Collection of cells equivalent to the original one but compatible with
        OpenMC

    """

    cv.check_type('opencg_cell', opencg_cell, opencg.Cell)
    cv.check_type('opencg_surface', opencg_surface, opencg.Surface)
    cv.check_value('halfspace', halfspace, (-1, +1))

    # Initialize an empty list for the new compatible cells
    compatible_cells = []

    # SquarePrism Surfaces
    if opencg_surface.type in ['x-squareprism', 'y-squareprism',
                               'z-squareprism']:

        # Get the compatible Surfaces (XPlanes and YPlanes)
        compatible_surfaces = get_compatible_opencg_surfaces(opencg_surface)

        opencg_cell.remove_surface(opencg_surface)

        # If Cell is inside SquarePrism, add "inside" of Surface halfspaces
        if halfspace == -1:
            opencg_cell.add_surface(compatible_surfaces[0], +1)
            opencg_cell.add_surface(compatible_surfaces[1], -1)
            opencg_cell.add_surface(compatible_surfaces[2], +1)
            opencg_cell.add_surface(compatible_surfaces[3], -1)
            compatible_cells.append(opencg_cell)

        # If Cell is outside SquarePrism, add "outside" of Surface halfspaces
        else:
            # Create 8 Cell clones to represent each of the disjoint planar
            # Surface halfspace intersections
            num_clones = 8

            for clone_id in range(num_clones):
                # Create cloned OpenCG Cell with Surfaces compatible with OpenMC
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
        cell.remove_redundant_surfaces()

    # Return the list of OpenMC compatible OpenCG Cells
    return compatible_cells


def make_opencg_cells_compatible(opencg_universe):
    """Make all cells in an OpenCG universe compatible with OpenMC.

    Parameters
    ----------
    opencg_universe : opencg.Universe
        Universe to check

    """

    cv.check_type('opencg_universe', opencg_universe, opencg.Universe)

    # Check all OpenCG Cells in this Universe for compatibility with OpenMC
    opencg_cells = opencg_universe.cells

    for cell_id, opencg_cell in opencg_cells.items():

        # Check each of the OpenCG Surfaces for OpenMC compatibility
        surfaces = opencg_cell.surfaces

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
                opencg_universe.remove_cell(opencg_cell)

                # Add the compatible OpenCG Cells to the Universe
                opencg_universe.add_cells(cells)

                # Make recursive call to look at the updated state of the
                # OpenCG Universe and return
                return make_opencg_cells_compatible(opencg_universe)

    # If all OpenCG Cells in the OpenCG Universe are compatible, return
    return


def get_openmc_cell(opencg_cell):
    """Return an OpenMC cell corresponding to an OpenCG cell.

    Parameters
    ----------
    opencg_cell : opencg.Cell
        OpenCG cell

    Returns
    -------
    openmc_cell : openmc.universe.Cell
        Equivalent OpenMC cell

    """

    cv.check_type('opencg_cell', opencg_cell, opencg.Cell)

    global OPENMC_CELLS
    cell_id = opencg_cell.id

    # If this Cell was already created, use it
    if cell_id in OPENMC_CELLS:
        return OPENMC_CELLS[cell_id]

    # Create an OpenCG Cell to represent this OpenMC Cell
    name = opencg_cell.name
    openmc_cell = openmc.Cell(cell_id, name)

    fill = opencg_cell.fill

    if opencg_cell.type == 'universe':
        openmc_cell.fill = get_openmc_universe(fill)
    elif opencg_cell.type == 'lattice':
        openmc_cell.fill = get_openmc_lattice(fill)
    else:
        openmc_cell.fill = get_openmc_material(fill)

    if opencg_cell.rotation is not None:
        rotation = np.asarray(opencg_cell.rotation, dtype=np.float64)
        openmc_cell.rotation = rotation

    if opencg_cell.translation is not None:
        translation = np.asarray(opencg_cell.translation, dtype=np.float64)
        openmc_cell.translation = translation

    surfaces = []
    operators = []
    for surface, halfspace in opencg_cell.surfaces.values():
        surfaces.append(get_openmc_surface(surface))
        operators.append(operator.neg if halfspace == -1 else operator.pos)
    openmc_cell.region = openmc.Intersection(
        *[op(s) for op, s in zip(operators, surfaces)])

    # Add the OpenMC Cell to the global collection of all OpenMC Cells
    OPENMC_CELLS[cell_id] = openmc_cell

    # Add the OpenCG Cell to the global collection of all OpenCG Cells
    OPENCG_CELLS[cell_id] = opencg_cell

    return openmc_cell


def get_opencg_universe(openmc_universe):
    """Return an OpenCG universe corresponding to an OpenMC universe.

    Parameters
    ----------
    openmc_universe : openmc.universe.Universe
        OpenMC universe

    Returns
    -------
    opencg_universe : opencg.Universe
        Equivalent OpenCG universe

    """

    cv.check_type('openmc_universe', openmc_universe, openmc.Universe)

    global OPENCG_UNIVERSES
    universe_id = openmc_universe.id

    # If this Universe was already created, use it
    if universe_id in OPENCG_UNIVERSES:
        return OPENCG_UNIVERSES[universe_id]

    # Create an OpenCG Universe to represent this OpenMC Universe
    name = openmc_universe.name
    opencg_universe = opencg.Universe(universe_id, name)

    # Convert all OpenMC Cells in this Universe to OpenCG Cells
    openmc_cells = openmc_universe.cells

    for cell_id, openmc_cell in openmc_cells.items():
        opencg_cell = get_opencg_cell(openmc_cell)
        opencg_universe.add_cell(opencg_cell)

    # Add the OpenMC Universe to the global collection of all OpenMC Universes
    OPENMC_UNIVERSES[universe_id] = openmc_universe

    # Add the OpenCG Universe to the global collection of all OpenCG Universes
    OPENCG_UNIVERSES[universe_id] = opencg_universe

    return opencg_universe


def get_openmc_universe(opencg_universe):
    """Return an OpenMC universe corresponding to an OpenCG universe.

    Parameters
    ----------
    opencg_universe : opencg.Universe
        OpenCG universe

    Returns
    -------
    openmc_universe : openmc.universe.Universe
        Equivalent OpenMC universe

    """

    cv.check_type('opencg_universe', opencg_universe, opencg.Universe)

    global OPENMC_UNIVERSES
    universe_id = opencg_universe.id

    # If this Universe was already created, use it
    if universe_id in OPENMC_UNIVERSES:
        return OPENMC_UNIVERSES[universe_id]

    # Make all OpenCG Cells and Surfaces in this Universe compatible with OpenMC
    make_opencg_cells_compatible(opencg_universe)

    # Create an OpenMC Universe to represent this OpenCSg Universe
    name = opencg_universe.name
    openmc_universe = openmc.Universe(universe_id, name)

    # Convert all OpenCG Cells in this Universe to OpenMC Cells
    opencg_cells = opencg_universe.cells

    for cell_id, opencg_cell in opencg_cells.items():
        openmc_cell = get_openmc_cell(opencg_cell)
        openmc_universe.add_cell(openmc_cell)

    # Add the OpenMC Universe to the global collection of all OpenMC Universes
    OPENMC_UNIVERSES[universe_id] = openmc_universe

    # Add the OpenCG Universe to the global collection of all OpenCG Universes
    OPENCG_UNIVERSES[universe_id] = opencg_universe

    return openmc_universe


def get_opencg_lattice(openmc_lattice):
    """Return an OpenCG lattice corresponding to an OpenMC lattice.

    Parameters
    ----------
    openmc_lattice : openmc.universe.Lattice
        OpenMC lattice

    Returns
    -------
    opencg_lattice : opencg.Lattice
        Equivalent OpenCG lattice

    """

    cv.check_type('openmc_lattice', openmc_lattice, openmc.Lattice)

    global OPENCG_LATTICES
    lattice_id = openmc_lattice.id

    # If this Lattice was already created, use it
    if lattice_id in OPENCG_LATTICES:
        return OPENCG_LATTICES[lattice_id]

    # Create an OpenCG Lattice to represent this OpenMC Lattice
    name = openmc_lattice.name
    dimension = openmc_lattice.dimension
    pitch = openmc_lattice.pitch
    lower_left = openmc_lattice.lower_left
    universes = openmc_lattice.universes
    outer = openmc_lattice.outer

    # Convert 2D dimension to 3D for OpenCG
    if len(dimension) == 2:
        new_dimension = np.ones(3, dtype=np.int)
        new_dimension[:2] = dimension
        dimension = new_dimension

    # Convert 2D pitch to 3D for OpenCG
    if len(pitch) == 2:
        new_pitch = np.ones(3, dtype=np.float64) * np.finfo(np.float64).max
        new_pitch[:2] = pitch
        pitch = new_pitch

    # Convert 2D lower left to 3D for OpenCG
    if len(lower_left) == 2:
        new_lower_left = np.ones(3, dtype=np.float64) * np.finfo(np.float64).min
        new_lower_left[:2] = lower_left
        lower_left = new_lower_left

    # Convert 2D universes array to 3D for OpenCG
    if len(universes.shape) == 2:
        new_universes = universes.copy()
        new_universes.shape = (1,) + universes.shape
        universes = new_universes

    # Initialize an empty array for the OpenCG nested Universes in this Lattice
    universe_array = np.empty(tuple(np.array(dimension)[::-1]),
                              dtype=opencg.Universe)

    # Create OpenCG Universes for each unique nested Universe in this Lattice
    unique_universes = openmc_lattice.get_unique_universes()

    for universe_id, universe in unique_universes.items():
        unique_universes[universe_id] = get_opencg_universe(universe)

    # Build the nested Universe array
    for z in range(dimension[2]):
        for y in range(dimension[1]):
            for x in range(dimension[0]):
                universe_id = universes[z][y][x].id
                universe_array[z][y][x] = unique_universes[universe_id]

    opencg_lattice = opencg.Lattice(lattice_id, name)
    opencg_lattice.dimension = dimension
    opencg_lattice.width = pitch
    opencg_lattice.universes = universe_array
    if outer is not None:
        opencg_lattice.outside = get_opencg_universe(outer)

    offset = np.array(lower_left, dtype=np.float64) - \
             ((np.array(pitch, dtype=np.float64) *
               np.array(dimension, dtype=np.float64))) / -2.0
    opencg_lattice.offset = offset

    # Add the OpenMC Lattice to the global collection of all OpenMC Lattices
    OPENMC_LATTICES[lattice_id] = openmc_lattice

    # Add the OpenCG Lattice to the global collection of all OpenCG Lattices
    OPENCG_LATTICES[lattice_id] = opencg_lattice

    return opencg_lattice


def get_openmc_lattice(opencg_lattice):
    """Return an OpenMC lattice corresponding to an OpenCG lattice.

    Parameters
    ----------
    opencg_lattice : opencg.Lattice
        OpenCG lattice

    Returns
    -------
    openmc_lattice : openmc.universe.Lattice
        Equivalent OpenMC lattice

    """

    cv.check_type('opencg_lattice', opencg_lattice, opencg.Lattice)

    global OPENMC_LATTICES
    lattice_id = opencg_lattice.id

    # If this Lattice was already created, use it
    if lattice_id in OPENMC_LATTICES:
        return OPENMC_LATTICES[lattice_id]

    dimension = opencg_lattice.dimension
    width = opencg_lattice.width
    offset = opencg_lattice.offset
    universes = opencg_lattice.universes
    outer = opencg_lattice.outside

    # Initialize an empty array for the OpenMC nested Universes in this Lattice
    universe_array = np.empty(tuple(np.array(dimension)[::-1]),
                              dtype=openmc.Universe)

    # Create OpenMC Universes for each unique nested Universe in this Lattice
    unique_universes = opencg_lattice.get_unique_universes()

    for universe_id, universe in unique_universes.items():
        unique_universes[universe_id] = get_openmc_universe(universe)

    # Build the nested Universe array
    for z in range(dimension[2]):
        for y in range(dimension[1]):
            for x in range(dimension[0]):
                universe_id = universes[z][y][x].id
                universe_array[z][y][x] = unique_universes[universe_id]

    # Reverse y-dimension in array to match ordering in OpenCG
    universe_array = universe_array[:, ::-1, :]

    lower_left = np.array(offset, dtype=np.float64) + \
                 ((np.array(width, dtype=np.float64) *
                   np.array(dimension, dtype=np.float64))) / -2.0

    openmc_lattice = openmc.RectLattice(lattice_id=lattice_id)
    openmc_lattice.pitch = width
    openmc_lattice.universes = universe_array
    openmc_lattice.lower_left = lower_left
    if outer is not None:
        openmc_lattice.outer = get_openmc_universe(outer)

    # Add the OpenMC Lattice to the global collection of all OpenMC Lattices
    OPENMC_LATTICES[lattice_id] = openmc_lattice

    # Add the OpenCG Lattice to the global collection of all OpenCG Lattices
    OPENCG_LATTICES[lattice_id] = opencg_lattice

    return openmc_lattice


def get_opencg_geometry(openmc_geometry):
    """Return an OpenCG geometry corresponding to an OpenMC geometry.

    Parameters
    ----------
    openmc_geometry : openmc.universe.Geometry
        OpenMC geometry

    Returns
    -------
    opencg_geometry : opencg.Geometry
        Equivalent OpenCG geometry

    """

    cv.check_type('openmc_geometry', openmc_geometry, openmc.Geometry)

    # Clear dictionaries and auto-generated IDs
    OPENMC_SURFACES.clear()
    OPENCG_SURFACES.clear()
    OPENMC_CELLS.clear()
    OPENCG_CELLS.clear()
    OPENMC_UNIVERSES.clear()
    OPENCG_UNIVERSES.clear()
    OPENMC_LATTICES.clear()
    OPENCG_LATTICES.clear()

    openmc_root_universe = openmc_geometry.root_universe
    opencg_root_universe = get_opencg_universe(openmc_root_universe)

    opencg_geometry = opencg.Geometry()
    opencg_geometry.root_universe = opencg_root_universe
    opencg_geometry.initialize_cell_offsets()
    opencg_geometry.assign_auto_ids()

    return opencg_geometry


def get_openmc_geometry(opencg_geometry):
    """Return an OpenMC geometry corresponding to an OpenCG geometry.

    Parameters
    ----------
    opencg_geometry : opencg.Geometry
        OpenCG geometry

    Returns
    -------
    openmc_geometry : openmc.universe.Geometry
        Equivalent OpenMC geometry

    """

    cv.check_type('opencg_geometry', opencg_geometry, opencg.Geometry)

    # Deep copy the goemetry since it may be modified to make all Surfaces
    # compatible with OpenMC's specifications
    opencg_geometry.assign_auto_ids()
    opencg_geometry = copy.deepcopy(opencg_geometry)

    # Update Cell bounding boxes in Geometry
    opencg_geometry.update_bounding_boxes()

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
    universes = opencg_geometry.get_all_universes()
    for universe_id, universe in universes.items():
        if not isinstance(universe, opencg.Lattice):
            make_opencg_cells_compatible(universe)

    opencg_geometry.assign_auto_ids()

    opencg_root_universe = opencg_geometry.root_universe
    openmc_root_universe = get_openmc_universe(opencg_root_universe)

    openmc_geometry = openmc.Geometry()
    openmc_geometry.root_universe = openmc_root_universe

    return openmc_geometry
