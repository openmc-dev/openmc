import numpy as np
import openmc
import pytest

from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def sphere_model():
    
    model = openmc.model.Model()

    # Define materials
    NaI = openmc.Material()
    NaI.set_density('g/cc', 3.7)
    NaI.add_element('Na', 1.0)
    NaI.add_element('I', 1.0)

    model.materials = openmc.Materials([NaI])

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
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point(),
        energy=openmc.stats.Discrete([1e6], [1]),
        particle='photon'
    )

    # Define tallies
    tally = openmc.Tally(name="geb pht tally")
    tally.gaussian_broadening = (1000., 4., 0.0002)
    tally.scores = ['pulse-height']
    cell_filter = openmc.CellFilter(inner_sphere)
    energy_filter = openmc.EnergyFilter(np.linspace(0, 1_000_000, 101))
    tally.filters = [cell_filter, energy_filter]
    model.tallies = [tally]

    return model



def test_gaussian_broadening(sphere_model):
    harness = PyAPITestHarness('statepoint.5.h5', sphere_model)
    harness.main()
