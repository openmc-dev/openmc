import openmc
import openmc.lib


def test_universe_instances(run_in_tmpdir):
    """Number of instances of cells in universes that are used through
    nested lattices and directly, compared with the Python API."""
    water = openmc.Material()
    water.add_nuclide('H1', 2.0)
    water.add_nuclide('O16', 1.0)
    water.set_density('g/cm3', 1.0)

    # Pin universes
    cyl = openmc.ZCylinder(r=0.4)
    fuel_pin = openmc.Universe(cells=[
        openmc.Cell(fill=water, region=-cyl),
        openmc.Cell(fill=water, region=+cyl)])
    guide_pin = openmc.Universe(cells=[openmc.Cell(fill=water)])

    # Two assembly types, used several times in a core lattice
    def assembly(pins):
        lat = openmc.RectLattice()
        lat.lower_left = (-1.5, -1.5)
        lat.pitch = (1.0, 1.0)
        lat.universes = pins
        return openmc.Universe(cells=[openmc.Cell(fill=lat)])

    f, g = fuel_pin, guide_pin
    assembly_a = assembly([[f, f, f], [f, g, f], [f, f, f]])
    assembly_b = assembly([[f, g, f], [g, g, g], [f, g, f]])
    core_lat = openmc.RectLattice()
    core_lat.lower_left = (-6.0, -3.0)
    core_lat.pitch = (3.0, 3.0)
    core_lat.universes = [[assembly_a, assembly_b, assembly_a, assembly_a],
                          [assembly_b, assembly_a, assembly_b, assembly_a]]

    # The fuel pin is also used directly next to the core
    box = openmc.model.RectangularParallelepiped(
        -6.0, 7.0, -3.0, 3.0, -1.0, 1.0, boundary_type='vacuum')
    x_core = openmc.XPlane(6.0)
    root = openmc.Universe(cells=[
        openmc.Cell(fill=core_lat, region=-box & -x_core),
        openmc.Cell(fill=fuel_pin, region=-box & +x_core)])

    model = openmc.Model()
    model.geometry = openmc.Geometry(root)
    model.materials = openmc.Materials([water])
    model.settings.particles = 10
    model.settings.batches = 1
    model.settings.run_mode = 'fixed source'
    model.export_to_model_xml()

    model.geometry.determine_paths()
    with openmc.lib.run_in_memory():
        for cell in model.geometry.get_all_cells().values():
            assert openmc.lib.cells[cell.id].num_instances == \
                cell.num_instances

        # Five A and three B assemblies, plus the pin used directly
        for cell in fuel_pin.cells.values():
            assert openmc.lib.cells[cell.id].num_instances == 5*8 + 3*4 + 1
        for cell in guide_pin.cells.values():
            assert openmc.lib.cells[cell.id].num_instances == 5*1 + 3*5
