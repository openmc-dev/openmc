import openmc


def test_negative_positron_heating():
    m = openmc.Material()
    m.add_element('Li', 1.0)
    m.set_density('g/cm3', 10.0)

    surf = openmc.Sphere(r=100.0, boundary_type='reflective')
    cell = openmc.Cell(fill=m, region=-surf)
    model = openmc.Model()
    model.geometry = openmc.Geometry([cell])
    model.settings.run_mode = 'fixed source'
    model.settings.source = openmc.Source(
        space=openmc.stats.Point(),
        energy=openmc.stats.Discrete([5.0e6], [1.0]),
        particle='photon',
    )
    model.settings.particles = 7
    model.settings.batches = 1
    model.settings.electron_treatment = 'led'
    model.settings.seed = 513836

    tally = openmc.Tally()
    tally.filters = [openmc.ParticleFilter(['photon', 'electron', 'positron'])]
    tally.scores = ['heating']
    model.tallies = openmc.Tallies([tally])
    model.run(apply_tally_results=True)

    assert (tally.mean >= 0.0).all(), "Negative heating detected"
