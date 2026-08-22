import openmc


def test_forced_collision(run_in_tmpdir):
    m = openmc.Material()
    m.add_element('Li', 1.0)
    m.set_density('g/cm3', 1.0)

    surf1 = openmc.Sphere(r=99.9)
    surf2 = openmc.Sphere(r=100, boundary_type='vacuum')
    
    cell1 = openmc.Cell(region=-surf1)
    cell2 = openmc.Cell(fill=m, region=-surf2 & +surf1)
    cell2.forced_collision = True
    model = openmc.Model()
    model.geometry = openmc.Geometry([cell1, cell2])
    model.settings.run_mode = 'fixed source'
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point(),
        energy=openmc.stats.Discrete([5.0e6], [1.0]),
        particle='photon',
    )
    model.settings.particles = 10
    model.settings.batches = 1

    tally = openmc.Tally()
    tally.filters = [openmc.CellFilter([cell2])]
    tally.scores = ['heating']
    model.tallies = openmc.Tallies([tally])
    model.run(apply_tally_results=True)

    assert (tally.mean > 0.0).all(), "No heating detected"
