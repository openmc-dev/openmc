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
    model.settings.source = openmc.IndependentSource(
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


def test_direct_ttb_heating_scores_as_photon():
    openmc.reset_auto_ids()

    m = openmc.Material()
    m.add_element('Pb', 1.0)
    m.set_density('g/cm3', 11.34)

    surf = openmc.Sphere(r=10.0, boundary_type='reflective')
    cell = openmc.Cell(fill=m, region=-surf)

    model = openmc.Model()
    model.geometry = openmc.Geometry([cell])
    model.settings.run_mode = 'fixed source'
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point(),
        energy=openmc.stats.Discrete([5.0e6], [1.0]),
        particle='photon',
    )
    model.settings.particles = 20
    model.settings.batches = 1
    model.settings.electron_treatment = 'ttb'
    model.settings.seed = 1

    tally = openmc.Tally()
    tally.filters = [openmc.ParticleFilter(['photon', 'electron', 'positron'])]
    tally.scores = ['heating']

    production = openmc.Tally()
    production.filters = [
        openmc.ParticleFilter('photon'),
        openmc.ParticleProductionFilter(['electron', 'positron']),
    ]
    production.scores = ['events']
    production.estimator = 'analog'

    model.tallies = openmc.Tallies([tally, production])
    model.run(apply_tally_results=True)

    heating = tally.mean.ravel()
    assert heating[0] > 0.0
    assert heating[1] == 0.0
    assert heating[2] == 0.0
    assert (production.mean == 0.0).all()
