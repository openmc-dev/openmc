import openmc


def test_zero_energy_pulse_heights(run_in_tmpdir):
    surf = openmc.Sphere(r=10.0, boundary_type='vacuum')
    cell = openmc.Cell(region=-surf)
    model = openmc.Model()
    model.geometry = openmc.Geometry([cell])
    model.settings.run_mode = 'fixed source'
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point(),
        energy=openmc.stats.Discrete([5.0e6], [1.0]),
        particle='photon',
    )
    model.settings.particles = 100
    model.settings.batches = 1
    model.settings.electron_treatment = 'led'

    tally = openmc.Tally()
    tally.scores = ['pulse-height']
    model.tallies = openmc.Tallies([tally])
    model.run(apply_tally_results=True)

    assert (tally.mean == 0.0).all(), "Zero energy pulses should not be counted"
