"""Test that weight-window splitting does not bias the IFP generation time."""

import numpy as np
import openmc

from tests import cdtemp
from tests.regression_tests import config


def test_ifp_generation_time_weight_windows():
    """Compare the IFP generation time of an analog run and of a run in which
    most fission neutrons are split by weight windows.

    The model is the U-235 sphere of the ifp regression tests. Weight windows
    with an upper bound of 0.1 are defined on every fourth cell of a fine mesh
    and are absent (-1) on the other cells. A fission neutron born where no
    window applies keeps its own weight as reference, so at its first
    collision in a cell with a window it is split into 10 neutrons. Split
    neutrons must keep the lifetime clock of the neutron they were split from.
    When they restarted it at zero, the generation time with weight windows
    was lower than without: by 14.3% with NNDC data, and by 13.4% on average
    over 5 seeds with ENDF/B-VII.1 data.
    """
    model = openmc.Model()

    material = openmc.Material()
    material.add_nuclide('U235', 1.0)
    material.set_density('g/cm3', 16.0)
    radius = 10.0
    sphere = openmc.Sphere(r=radius, boundary_type='vacuum')
    cell = openmc.Cell(region=-sphere, fill=material)
    model.geometry = openmc.Geometry([cell])

    model.settings.particles = 50000
    model.settings.batches = 30
    model.settings.inactive = 10
    model.settings.ifp_n_generation = 1
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Box(*cell.bounding_box),
        constraints={'fissionable': True})

    tally = openmc.Tally()
    tally.scores = ['ifp-time-numerator', 'ifp-denominator']
    model.tallies = [tally]

    mesh = openmc.RegularMesh()
    mesh.lower_left = (-radius, -radius, -radius)
    mesh.upper_right = (radius, radius, radius)
    mesh.dimension = (20, 20, 20)
    i, j, k = np.indices(mesh.dimension)
    has_window = (i + j + k) % 4 == 0
    weight_windows = openmc.WeightWindows(
        mesh,
        lower_ww_bounds=np.where(has_window, 0.05, -1.0),
        upper_ww_bounds=np.where(has_window, 0.1, -1.0),
        energy_bounds=[0.0, 2.0e7],
        particle_type='neutron',
    )

    generation_time = {}
    with cdtemp():
        for weight_windows_on in (False, True):
            if weight_windows_on:
                model.settings.weight_windows = [weight_windows]
                model.settings.weight_windows_on = True
            sp_path = model.run(event_based=config['event'])
            with openmc.StatePoint(sp_path) as sp:
                kinetics = sp.get_kinetics_parameters()
                generation_time[weight_windows_on] = \
                    kinetics.generation_time.nominal_value

    # Across seeds, the standard deviation of this ratio is about 0.5%
    ratio = generation_time[True] / generation_time[False]
    assert abs(ratio - 1.0) <= 0.025, (
        f'IFP generation time with weight windows ({generation_time[True]:.4e}'
        f' s) differs from the analog value ({generation_time[False]:.4e} s)'
    )
