import openmc


def test_musurface(run_in_tmpdir):
    sphere = openmc.Sphere(r=1.0, boundary_type='vacuum')
    cell = openmc.Cell(region=-sphere)
    model = openmc.Model()
    model.geometry = openmc.Geometry([cell])
    model.settings.particles = 1000
    model.settings.batches = 10
    E = 1.0
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point(),
        angle=openmc.stats.Isotropic(),
        energy=openmc.stats.delta_function(E),
    )
    model.settings.run_mode = "fixed source"

    filter1 = openmc.MuSurfaceFilter(200)
    filter2 = openmc.SurfaceFilter(sphere)
    tally = openmc.Tally()
    tally.filters = [filter1, filter2]
    tally.scores = ['current']
    model.tallies = openmc.Tallies([tally])

    # Run OpenMC
    sp_filename = model.run()

    # Get current binned by mu
    with openmc.StatePoint(sp_filename) as sp:
        current_mu = sp.tallies[tally.id].mean.ravel()

    # All contributions should show up in last bin
    assert current_mu[-1] == 1.0
    for element in current_mu[:-1]:
        assert element == 0.0


