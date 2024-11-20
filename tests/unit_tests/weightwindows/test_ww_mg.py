import pytest

import numpy as np
import openmc


def test_weight_windows_mg(run_in_tmpdir):
    model = openmc.examples.random_ray_three_region_cube()

    mesh = openmc.RegularMesh.from_domain(model.geometry, (3, 3, 3))
    mesh_tally = openmc.Tally()
    mesh_tally.filters = [openmc.MeshFilter(mesh)]
    mesh_tally.scores = ['flux']
    model.tallies = [mesh_tally]

    settings = openmc.Settings()
    settings.particles = 1000
    settings.batches = 10
    settings.energy_mode = 'multi-group'
    settings.run_mode = 'fixed source'
    space = openmc.stats.Point((1, 1, 1))
    energy = openmc.stats.delta_function(1e6)
    source = openmc.IndependentSource(space=space)
    settings.source = source

    model.settings = settings
    statepoint = model.run()

    with openmc.StatePoint(statepoint) as sp:
        tally_out = sp.get_tally(id=mesh_tally.id)
        flux_analog = tally_out.mean

    # wwg = openmc.WeightWindowGenerator(mesh)
    # model.settings.weight_window_generators = wwg
    # model.settings.weight_window_checkpoints = {'surface': True, 'collision': True}
    # model.settings.survival_biasing = False

    # statepoint = model.run()

    ww_lower_bnds = np.loadtxt('ww_mg.txt')
    weight_windows = openmc.WeightWindows(mesh, lower_ww_bounds=ww_lower_bnds, upper_bound_ratio=5.0)
    model.settings.weight_windows = weight_windows
    model.settings.weight_windows_on = True

    # weight_windows = openmc.hdf5_to_wws('weight_windows.h5')
    # for ww in weight_windows:
    #     ww.lower_ww_bounds *= 0.1
    #     ww.upper_ww_bounds *= 0.1

    # np.savetxt('ww_mg.txt', weight_windows[0].lower_ww_bounds.flatten())

    settings.weight_windows = weight_windows

    statepoint = model.run()

    with openmc.StatePoint(statepoint) as sp:
        tally_out = sp.get_tally(id=mesh_tally.id)
        flux_ww = tally_out.mean

    print(flux_analog.sum())
    print(flux_ww.sum())
    assert np.allclose(flux_analog.sum(), flux_ww.sum(), rtol=1e-2)

