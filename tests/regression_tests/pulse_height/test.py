import numpy as np
import openmc
import pytest
from openmc.utility_funcs import change_directory

from tests.testing_harness import PyAPITestHarness


@pytest.mark.parametrize("shared_secondary,particle", [
    (False, "photon"),
    (False, "neutron"),
    (True, "photon"),
    (True, "neutron")
])
def test_pulse_height(shared_secondary, particle):
    subdir = f"shared_{particle}" if shared_secondary else f"local_{particle}"
    with change_directory(subdir):
        openmc.reset_auto_ids()
        model = openmc.Model()

        # Define materials
        NaI = openmc.Material()
        NaI.set_density('g/cm3', 3.7)
        NaI.add_element('Na', 1.0)
        NaI.add_element('I', 1.0)

        # Define geometry: two spheres in each other
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
        model.settings.source = openmc.IndependentSource(
            energy=openmc.stats.delta_function(1e6),
            particle=particle
        )

        # Define tallies
        tally = openmc.Tally(name="pht tally")
        tally.scores = ['pulse-height']
        cell_filter = openmc.CellFilter(inner_sphere)
        energy_filter = openmc.EnergyFilter(np.linspace(0, 1e6, 101))
        tally.filters = [cell_filter, energy_filter]
        model.tallies = [tally]

        harness = PyAPITestHarness('statepoint.5.h5', model)
        harness.main()
