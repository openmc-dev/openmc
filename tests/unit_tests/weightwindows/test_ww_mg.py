import pytest

import numpy as np
import openmc


def test_weight_windows_mg(request, run_in_tmpdir):
    # import basic random ray model
    model = openmc.examples.random_ray_three_region_cube()

    # create a mesh tally
    mesh = openmc.RegularMesh.from_domain(model.geometry, (3, 3, 3))
    mesh_tally = openmc.Tally()
    mesh_tally.filters = [openmc.MeshFilter(mesh)]
    mesh_tally.scores = ['flux']
    model.tallies = [mesh_tally]

    # replace random ray settings with fixed source settings
    settings = openmc.Settings()
    settings.particles = 5000
    settings.batches = 10
    settings.energy_mode = 'multi-group'
    settings.run_mode = 'fixed source'
    space = openmc.stats.Point((1, 1, 1))
    energy = openmc.stats.delta_function(1e6)
    source = openmc.IndependentSource(space=space)
    settings.source = source

    model.settings = settings
    statepoint = model.run()

    # extract flux from analog simulation
    with openmc.StatePoint(statepoint) as sp:
        tally_out = sp.get_tally(id=mesh_tally.id)
        flux_analog = tally_out.mean

    # load the weight windows for this problem, apply them, and re-run
    ww_lower_bnds = np.loadtxt(request.path.parent / 'ww_mg.txt')
    weight_windows = openmc.WeightWindows(mesh, lower_ww_bounds=ww_lower_bnds, upper_bound_ratio=5.0)
    model.settings.weight_windows = weight_windows
    model.settings.weight_windows_on = True
    settings.weight_windows = weight_windows

    statepoint = model.run()
    with openmc.StatePoint(statepoint) as sp:
        tally_out = sp.get_tally(id=mesh_tally.id)
        flux_ww = tally_out.mean

    # the sum of the fluxes should approach the same value
    analog_sum = flux_analog.sum()
    ww_sum = flux_ww.sum()
    assert np.allclose(analog_sum, ww_sum, rtol=1e-2)

