import openmc


def test_periodic_surface_roundtrip(run_in_tmpdir):
    # Create a simple model with periodic surfaces
    mat = openmc.Material()
    mat.add_nuclide('H1', 1.0)
    mat.set_density('g/cm3', 1.0)
    cyl = openmc.ZCylinder(r=1.0)
    x0 = openmc.XPlane(-5.0, boundary_type='periodic')
    y0 = openmc.YPlane(-5.0, boundary_type='periodic')
    z0 = openmc.ZPlane(-5.0, boundary_type='periodic')
    x1 = openmc.XPlane(5.0, boundary_type='periodic')
    y1 = openmc.YPlane(5.0, boundary_type='periodic')
    z1 = openmc.ZPlane(5.0, boundary_type='periodic')
    x0.periodic_surface = x1
    y0.periodic_surface = y1
    z0.periodic_surface = z1
    cell1 = openmc.Cell(fill=mat, region=-cyl)
    cell2 = openmc.Cell(fill=mat, region=+cyl & +x0 & -x1 & +y0 & -y1 & +z0 & -z1)
    model = openmc.Model()
    model.geometry = openmc.Geometry([cell1, cell2])
    model.settings.particles = 100
    model.settings.batches = 1
    model.settings.run_mode = 'fixed source'
    model.settings.source = openmc.IndependentSource(
        energy=openmc.stats.delta_function(1.0e4)
    )

    # Run model
    model.run()

    # Load summary data and check periodic surfaces
    summary = openmc.Summary('summary.h5')
    surfs = summary.geometry.get_all_surfaces()
    for s in [x0, y0, z0, x1, y1, z1]:
        assert surfs[s.id].boundary_type == 'periodic'
    pairs = [(x0, x1), (y0, y1), (z0, z1)]
    for s0, s1 in pairs:
        assert surfs[s0.id].periodic_surface == surfs[s1.id]
        assert surfs[s1.id].periodic_surface == surfs[s0.id]
