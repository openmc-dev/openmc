import numpy as np
import openmc
import pytest
from openmc.utility_funcs import change_directory

from tests.testing_harness import PyAPITestHarness


@pytest.mark.parametrize("shared_secondary,subdir", [
    (False, "local"),
    (True, "shared"),
])
def test_weightwindows_pulse_height(shared_secondary, subdir):
    with change_directory(subdir):
        openmc.reset_auto_ids()
        model = openmc.Model()

        # Define materials (NaI scintillator)
        NaI = openmc.Material()
        NaI.set_density('g/cc', 3.7)
        NaI.add_element('Na', 1.0)
        NaI.add_element('I', 1.0)

        model.materials = openmc.Materials([NaI])

        # Define geometry: NaI sphere inside vacuum sphere
        s1 = openmc.Sphere(r=1)
        s2 = openmc.Sphere(r=2, boundary_type='vacuum')
        inner_sphere = openmc.Cell(name='inner sphere', fill=NaI, region=-s1)
        outer_sphere = openmc.Cell(name='outer sphere', region=+s1 & -s2)
        model.geometry = openmc.Geometry([inner_sphere, outer_sphere])

        # Define settings
        model.settings.run_mode = 'fixed source'
        model.settings.batches = 5
        model.settings.particles = 100
        model.settings.photon_transport = True
        model.settings.shared_secondary_bank = shared_secondary
        model.settings.max_history_splits = 50
        model.settings.weight_window_checkpoints = {
            'surface': True,
            'collision': True,
        }
        model.settings.source = openmc.IndependentSource(
            energy=openmc.stats.delta_function(1e6),
            particle='photon'
        )

        # Define pulse-height tally
        tally = openmc.Tally(name="pht tally")
        tally.scores = ['pulse-height']
        cell_filter = openmc.CellFilter(inner_sphere)
        energy_filter = openmc.EnergyFilter(np.linspace(0, 1_000_000, 101))
        tally.filters = [cell_filter, energy_filter]
        model.tallies = [tally]

        # Define weight windows on a simple mesh covering the geometry
        ww_mesh = openmc.RegularMesh()
        ww_mesh.lower_left = (-2, -2, -2)
        ww_mesh.upper_right = (2, 2, 2)
        ww_mesh.dimension = (1, 1, 1)

        # Single energy bin for photons
        e_bnds = [0.0, 2e6]

        # Uniform weight window bounds (low enough to trigger some splitting)
        lower_bounds = np.array([0.01])

        ww = openmc.WeightWindows(ww_mesh, lower_bounds, None, 5.0, e_bnds, 'photon')
        model.settings.weight_windows = [ww]

        harness = PyAPITestHarness('statepoint.5.h5', model)
        harness.main()
