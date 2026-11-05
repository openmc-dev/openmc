import openmc


FIELD_LOWER_LEFT = (0.0, 0.0, 0.0)
FIELD_UPPER_RIGHT = (5.0, 5.0, 5.0)
OUTER_LOWER_LEFT = (-5.0, -5.0, -5.0)
OUTER_UPPER_RIGHT = (10.0, 10.0, 10.0)


def _temperature_field_components(material_name=None):
    model = openmc.Model()

    material = openmc.Material(name=material_name)
    material.add_nuclide("U235", 0.2)
    material.add_nuclide("U238", 0.8)
    material.add_element("O", 2.0)
    material.add_element("H", 4.0)
    material.set_density("g/cm3", 5.0)
    material.add_s_alpha_beta("c_H_in_H2O")
    model.materials = [material]

    mesh = openmc.RegularMesh()
    mesh.lower_left = FIELD_LOWER_LEFT
    mesh.upper_right = FIELD_UPPER_RIGHT
    mesh.dimension = (2, 2, 2)
    temperatures = [294.0 + 100.0 * i for i in range(8)]

    temperature_field = openmc.TemperatureField(mesh, temperatures)
    return model, material, mesh, temperature_field


def _finish_model(model, mesh, temperature_field, lower_left, upper_right):
    settings = openmc.Settings()
    settings.batches = 20
    settings.particles = 200
    if temperature_field is not None:
        settings.temperature_field = temperature_field
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Box(lower_left, upper_right),
        constraints={"fissionable": True},
    )
    settings.temperature = {"tolerance": 1000, "multipole": True}
    model.settings = settings

    mesh_tally = openmc.Tally(name="total reaction rate")
    mesh_tally.filters = [openmc.MeshFilter(mesh)]
    mesh_tally.scores = ["total"]
    model.tallies = [mesh_tally]
    return model


def _box(lower_left, upper_right, **kwargs):
    return openmc.model.RectangularParallelepiped(
        lower_left[0],
        upper_right[0],
        lower_left[1],
        upper_right[1],
        lower_left[2],
        upper_right[2],
        **kwargs,
    )


def _temperature_cells(material, box):
    """Explicit cell temperatures, independently ordered with x varying fastest."""
    xmid = openmc.XPlane(x0=2.5)
    ymid = openmc.YPlane(y0=2.5)
    zmid = openmc.ZPlane(z0=2.5)
    cells = []
    for z in range(2):
        for y in range(2):
            for x in range(2):
                region = (-box & (+xmid if x else -xmid)
                          & (+ymid if y else -ymid)
                          & (+zmid if z else -zmid))
                cell = openmc.Cell(fill=material, region=region)
                cell.temperature = 294.0 + 100.0 * (x + 2*y + 4*z)
                cells.append(cell)
    return cells


def csg_single_cube_reflective(reference=False):
    model, material, mesh, temperature_field = _temperature_field_components()
    box = _box(FIELD_LOWER_LEFT, FIELD_UPPER_RIGHT, boundary_type="reflective")
    model.geometry = openmc.Geometry([openmc.Cell(fill=material, region=-box)])
    if reference:
        model.geometry = openmc.Geometry(_temperature_cells(material, box))
        temperature_field = None
    return _finish_model(
        model, mesh, temperature_field, FIELD_LOWER_LEFT, FIELD_UPPER_RIGHT
    )


def csg_nested_cubes_periodic(reference=False):
    model, material, mesh, temperature_field = _temperature_field_components()
    field_box = _box(FIELD_LOWER_LEFT, FIELD_UPPER_RIGHT)
    field_cell = openmc.Cell(fill=material, region=-field_box)

    xmin = openmc.XPlane(x0=OUTER_LOWER_LEFT[0], boundary_type="periodic")
    xmax = openmc.XPlane(x0=OUTER_UPPER_RIGHT[0], boundary_type="periodic")
    ymin = openmc.YPlane(y0=OUTER_LOWER_LEFT[1], boundary_type="periodic")
    ymax = openmc.YPlane(y0=OUTER_UPPER_RIGHT[1], boundary_type="periodic")
    zmin = openmc.ZPlane(z0=OUTER_LOWER_LEFT[2], boundary_type="periodic")
    zmax = openmc.ZPlane(z0=OUTER_UPPER_RIGHT[2], boundary_type="periodic")
    xmin.periodic_surface = xmax
    ymin.periodic_surface = ymax
    zmin.periodic_surface = zmax

    outer_region = +xmin & -xmax & +ymin & -ymax & +zmin & -zmax
    outer_cell = openmc.Cell(fill=material, region=outer_region & +field_box)
    outer_cell.temperature = 250.0
    model.geometry = openmc.Geometry([field_cell, outer_cell])
    if reference:
        model.geometry = openmc.Geometry(
            _temperature_cells(material, field_box) + [outer_cell])
        temperature_field = None
    return _finish_model(
        model, mesh, temperature_field, OUTER_LOWER_LEFT, OUTER_UPPER_RIGHT
    )


def lattice_nested_cubes_vacuum(reference=False):
    model, material, mesh, temperature_field = _temperature_field_components()
    lattice_box = openmc.model.RectangularParallelepiped(
        -1.75, 1.75, -1.75, 1.75, -1.75, 1.75
    )
    universe = openmc.Universe(
        cells=[openmc.Cell(fill=material, region=-lattice_box)]
    )
    lattice = openmc.RectLattice()
    lattice.lower_left = FIELD_LOWER_LEFT
    lattice.pitch = (2.5, 2.5, 2.5)
    lattice.universes = [
        [[universe, universe], [universe, universe]],
        [[universe, universe], [universe, universe]],
    ]
    lattice.outer = openmc.Universe(cells=[openmc.Cell(fill=material)])
    if reference:
        universes = []
        for z in range(2):
            layer = []
            # RectLattice arrays list y rows from high to low.
            for y in (1, 0):
                row = []
                for x in range(2):
                    cell = openmc.Cell(fill=material, region=-lattice_box)
                    cell.temperature = 294.0 + 100.0 * (x + 2*y + 4*z)
                    row.append(openmc.Universe(cells=[cell]))
                layer.append(row)
            universes.append(layer)
        lattice.universes = universes
        outer_cell = openmc.Cell(fill=material)
        outer_cell.temperature = 250.0
        lattice.outer = openmc.Universe(cells=[outer_cell])
        temperature_field = None

    outer_box = _box(OUTER_LOWER_LEFT, OUTER_UPPER_RIGHT, boundary_type="vacuum")
    outer_cell = openmc.Cell(fill=lattice, region=-outer_box)
    if not reference:
        outer_cell.temperature = 250.0
    model.geometry = openmc.Geometry([outer_cell])
    return _finish_model(
        model, mesh, temperature_field, OUTER_LOWER_LEFT, OUTER_UPPER_RIGHT
    )


def dagmc_nested_cubes_reflective(filename, reference=False):
    model, _, mesh, temperature_field = _temperature_field_components(
        material_name="fuel"
    )
    model.geometry = openmc.Geometry(
        openmc.DAGMCUniverse(filename=filename, auto_geom_ids=True)
    )
    if reference:
        temperature_field = None
    return _finish_model(
        model, mesh, temperature_field, OUTER_LOWER_LEFT, OUTER_UPPER_RIGHT
    )
