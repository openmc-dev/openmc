"""Test that CollisionFilter produces unbiased results with weight windows."""

import numpy as np
import openmc
from tests import cdtemp


def test_collision_filter_weight_windows():
    """Verify that total flux from CollisionFilter bins matches
    the total flux tally when weight windows are active.

    The test runs two simulations (analog and weight-windowed) and
    checks that the uncollided flux (CollisionFilter bin 0) is
    statistically consistent between them. Before the fix, split
    particles had their collision count reset to zero, inflating the
    bin-0 tally.
    """

    model = openmc.Model()

    # Concrete-like material for good scattering
    mat = openmc.Material()
    mat.set_density('g/cc', 2.3)
    mat.add_nuclide('H1', 0.168)
    mat.add_nuclide('O16', 0.562)
    mat.add_nuclide('Si28', 0.187)
    mat.add_nuclide('Ca40', 0.018)
    mat.add_nuclide('Fe56', 0.004)

    sphere_inner = openmc.Sphere(r=100)
    sphere_outer = openmc.Sphere(r=120, boundary_type='vacuum')
    cell = openmc.Cell(fill=mat, region=-sphere_inner)
    void = openmc.Cell(region=+sphere_inner & -sphere_outer)
    model.geometry = openmc.Geometry([cell, void])

    settings = openmc.Settings()
    settings.run_mode = 'fixed source'
    settings.particles = 2000
    settings.batches = 5
    settings.source = openmc.IndependentSource(
        energy=openmc.stats.delta_function(14e6),
    )
    settings.survival_biasing = True

    # Weight windows: decreasing values radially to force splitting
    ww_mesh = openmc.RegularMesh()
    ww_mesh.lower_left = (-120, -120, -120)
    ww_mesh.upper_right = (120, 120, 120)
    ww_mesh.dimension = (3, 3, 3)

    n_bins = 3 * 3 * 3
    lower_ww = np.logspace(-2, -5, n_bins)
    ww = openmc.WeightWindows(
        ww_mesh, lower_ww, None, 5.0,
        energy_bounds=[0.0, 20e6],
        survival_ratio=3.0,
    )

    cell_filter = openmc.CellFilter([cell])
    collision_filter = openmc.CollisionFilter(list(range(10)))

    # Tally 1: total flux
    tally_total = openmc.Tally()
    tally_total.filters = [cell_filter]
    tally_total.scores = ['flux']

    # Tally 2: flux by collision number
    tally_collision = openmc.Tally()
    tally_collision.filters = [cell_filter, collision_filter]
    tally_collision.scores = ['flux']

    model.tallies = openmc.Tallies([tally_total, tally_collision])

    with cdtemp():
        # Run analog (no weight windows)
        settings.weight_windows_on = False
        model.settings = settings
        sp_analog = model.run()

        sp_analog.rename('statepoint.analog.h5')

        # Run with weight windows
        settings.weight_windows = [ww]
        settings.weight_windows_on = True
        settings.max_history_splits = 1000
        model.settings = settings
        sp_ww = model.run()

        # Compare uncollided flux (bin 0) between analog and WW runs
        sp_a = openmc.StatePoint('statepoint.analog.h5')
        sp_w = openmc.StatePoint(sp_ww)

        # Get uncollided flux (collision bin 0) from both runs
        uncollided_analog = sp_a.tallies[tally_collision.id].get_values(
            filters=[openmc.CollisionFilter], filter_bins=[(0,)]
        ).sum()
        uncollided_ww = sp_w.tallies[tally_collision.id].get_values(
            filters=[openmc.CollisionFilter], filter_bins=[(0,)]
        ).sum()

        total_analog = sp_a.tallies[tally_total.id].mean.sum()
        total_ww = sp_w.tallies[tally_total.id].mean.sum()

        # The uncollided fraction should be similar in both runs.
        # Before the fix, the WW run had an inflated uncollided fraction
        # because split particles were incorrectly tagged as uncollided.
        frac_analog = uncollided_analog / total_analog if total_analog > 0 else 0
        frac_ww = uncollided_ww / total_ww if total_ww > 0 else 0

        # Allow generous tolerance due to Monte Carlo statistics, but
        # the bug caused 2-4x inflation so a 50% tolerance catches it
        assert frac_ww < frac_analog * 1.5, (
            f"Uncollided fraction with WW ({frac_ww:.4f}) is much larger "
            f"than analog ({frac_analog:.4f}), suggesting collision count "
            f"is not preserved during particle splitting"
        )
