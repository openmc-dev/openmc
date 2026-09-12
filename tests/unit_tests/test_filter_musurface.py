import math

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


def test_musurface_rotated_universe(run_in_tmpdir):
    """MuSurfaceFilter uses the surface normal in the root coordinate frame.

    The plane lives in a universe filled into a cell rotated 45 degrees about
    z, so its normal in the plane's own frame is (1, 0, 0) while the particle
    direction is stored in the root frame. Binning must use the root-frame
    normal, giving mu = cos(45 deg) rather than 1.
    """
    openmc.reset_auto_ids()

    xplane = openmc.XPlane(0.0)
    inner1 = openmc.Cell(region=-xplane)
    inner2 = openmc.Cell(region=+xplane)
    inner_univ = openmc.Universe(cells=[inner1, inner2])

    sph = openmc.Sphere(r=10.0, boundary_type='vacuum')
    root_cell = openmc.Cell(region=-sph, fill=inner_univ)
    root_cell.rotation = (0.0, 0.0, 45.0)

    model = openmc.Model()
    model.geometry = openmc.Geometry([root_cell])

    src = openmc.IndependentSource()
    src.space = openmc.stats.Point((-5.0, 0.0, 0.0))
    src.angle = openmc.stats.Monodirectional((1.0, 0.0, 0.0))

    model.settings.run_mode = 'fixed source'
    model.settings.batches = 1
    model.settings.particles = 100
    model.settings.source = src

    # 20 equal-width bins from -1 to 1; cos(45 deg) = 0.7071 falls in [0.7, 0.8)
    tally = openmc.Tally()
    tally.filters = [
        openmc.MuSurfaceFilter(20),
        openmc.SurfaceFilter([xplane]),
    ]
    tally.scores = ['current']
    model.tallies = [tally]

    model.run(apply_tally_results=True)
    current_mu = tally.mean.ravel()

    expected_bin = int((math.cos(math.radians(45.0)) + 1.0) / 0.1)
    assert current_mu[expected_bin] == 1.0
    assert current_mu.sum() == 1.0
