import openmc


def test_get_tally_filter_type():
    """Test various ways of retrieving tallies from a StatePoint object."""

    mat = openmc.Material()
    mat.add_nuclide('H1', 1.0)
    mat.set_density('g/cm3', 10.0)

    sphere = openmc.Sphere(r=10.0, boundary_type='vacuum')
    cell = openmc.Cell(fill=mat, region=-sphere)
    model = openmc.Model()
    model.geometry = openmc.Geometry([cell])
    model.settings.particles = 10
    model.settings.batches = 2
    model.settings.run_mode = 'fixed source'

    reg_mesh = openmc.RegularMesh().from_domain(cell)
    tally1 = openmc.Tally(tally_id=1)
    mesh_filter = openmc.MeshFilter(reg_mesh)
    tally1.filters = [mesh_filter]
    tally1.scores = ['flux']

    tally2 = openmc.Tally(tally_id=2, name='heating tally')
    tally2.filters = [openmc.CellFilter(cell)]
    tally2.scores = ['heating']

    model.tallies = openmc.Tallies([tally1, tally2])

    sp_filename = model.run()

    sp = openmc.StatePoint(sp_filename)

    tally_found = sp.get_tally(filter_type=openmc.MeshFilter)
    assert tally_found.id == 1

    tally_found = sp.get_tally(filter_type=openmc.CellFilter)
    assert tally_found.id == 2

    tally_found = sp.get_tally(filters=[mesh_filter])
    assert tally_found.id == 1

    tally_found = sp.get_tally(scores=['heating'])
    assert tally_found.id == 2

    tally_found = sp.get_tally(name='heating tally')
    assert tally_found.id == 2

    tally_found = sp.get_tally(name=None)
    assert tally_found.id == 1

    tally_found = sp.get_tally(id=1)
    assert tally_found.id == 1

    tally_found = sp.get_tally(id=2)
    assert tally_found.id == 2
