import openmc
import pytest


@pytest.mark.parametrize("direction", [1.0, -1.0])
def test_ww_surface_on_mesh_boundary(run_in_tmpdir, direction):
    """Weight windows applied at a surface that coincides with a weight window
    mesh boundary must use the mesh element the particle is entering,
    regardless of the direction of travel.

    At the surface checkpoint the particle sits exactly on the plane x = 0.
    Structured mesh index lookups resolve a position exactly on an element
    boundary to the element on the lower-coordinate side, so when the window
    was looked up at the particle's position alone, a particle moving in +x
    was given the window of the element it was leaving. It was then not split
    and the +x case failed, while the -x case passed only because there the
    lower-coordinate element happens to be the one being entered. The lookup
    now nudges the position along the direction of travel.
    """

    # Void-like one-group material so the only weight window checks happen at
    # surface crossings
    groups = openmc.mgxs.EnergyGroups([0.0, 20.0e6])
    xsdata = openmc.XSdata('void', groups)
    xsdata.order = 0
    xsdata.set_total([0.0])
    xsdata.set_absorption([0.0])
    xsdata.set_scatter_matrix([[[0.0]]])
    mg_library = openmc.MGXSLibrary(groups)
    mg_library.add_xsdata(xsdata)
    mg_library.export_to_hdf5('mgxs.h5')

    mat = openmc.Material()
    mat.add_macroscopic('void')
    materials = openmc.Materials([mat])
    materials.cross_sections = 'mgxs.h5'

    # Two cells separated by the plane x = 0, which is also a mesh boundary
    x_min = openmc.XPlane(-1.0, boundary_type='vacuum')
    x_mid = openmc.XPlane(0.0)
    x_max = openmc.XPlane(1.0, boundary_type='vacuum')
    yz = openmc.model.RectangularPrism(2.0, 2.0, axis='x',
                                       boundary_type='vacuum')
    left = openmc.Cell(fill=mat, region=+x_min & -x_mid & -yz)
    right = openmc.Cell(fill=mat, region=+x_mid & -x_max & -yz)
    geometry = openmc.Geometry([left, right])

    # Particle is born in the window of its birth element and must be split
    # five ways when it enters the element on the other side of x = 0
    birth_cell, target_cell = (left, right) if direction > 0 else (right, left)
    mesh = openmc.RectilinearMesh()
    mesh.x_grid = [-1.0, 0.0, 1.0]
    mesh.y_grid = [-1.0, 1.0]
    mesh.z_grid = [-1.0, 1.0]
    birth_bounds = (0.5, 1.5)
    target_bounds = (0.1, 0.2)
    if direction > 0:
        lower = [birth_bounds[0], target_bounds[0]]
        upper = [birth_bounds[1], target_bounds[1]]
    else:
        lower = [target_bounds[0], birth_bounds[0]]
        upper = [target_bounds[1], birth_bounds[1]]
    ww = openmc.WeightWindows(mesh, lower_ww_bounds=lower,
                              upper_ww_bounds=upper)

    settings = openmc.Settings()
    settings.energy_mode = 'multi-group'
    settings.run_mode = 'fixed source'
    settings.particles = 10
    settings.batches = 1
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Point((-0.5*direction, 0.0, 0.0)),
        angle=openmc.stats.Monodirectional((direction, 0.0, 0.0)),
    )
    settings.weight_windows = ww
    settings.weight_window_checkpoints = {'collision': False, 'surface': True}

    tally = openmc.Tally()
    tally.filters = [
        openmc.CellFilter([target_cell]),
        openmc.WeightFilter([0.0, 0.5, 2.0])
    ]
    tally.scores = ['flux']

    model = openmc.Model(geometry, materials, settings, openmc.Tallies([tally]))
    model.run(apply_tally_results=True)

    # All flux in the target cell should be carried by split particles with
    # weight 0.2; total flux is conserved (1 cm path length per source particle)
    flux = tally.mean.squeeze()
    assert flux == pytest.approx([1.0, 0.0])
