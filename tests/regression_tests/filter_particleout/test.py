import numpy as np
import openmc
import pytest

from tests.testing_harness import PyAPITestHarness

@pytest.fixture
def model():
    openmc.reset_auto_ids()
    model = openmc.Model()

    # =============================================================================
    # Materials
    # =============================================================================

    u238 = openmc.Material(name='U-238')
    u238.add_nuclide('U238', 1.0)
    u238.set_density('g/cm3', 19.1)
    model.materials = openmc.Materials([u238])
    
    # =============================================================================
    # Geometry
    # =============================================================================

    sphere = openmc.Sphere(r=1e90, boundary_type='vacuum')
    cell = openmc.Cell(fill=u238, region=-sphere)
    universe = openmc.Universe(cells=[cell])
    model.geometry = openmc.Geometry(universe)

    # =============================================================================
    # Settings
    # =============================================================================

    source = openmc.Source()
    source.space = openmc.stats.Point()
    source.energy = openmc.stats.Discrete([2.0e6], [1.0])
    model.settings = openmc.Settings()
    model.settings.run_mode = 'fixed source'
    model.settings.particles = 10_000
    model.settings.batches = 10
    model.settings.source = source
    model.settings.photon_transport = True
    
    # =============================================================================
    # Tallies
    # =============================================================================

    energy_bins = np.logspace(4, 8, 100)  # 10 keV to 100 MeV
    particle_filter = openmc.ParticleFilter(['neutron'])
    particleout_filter = openmc.ParticleoutFilter(['photon'])
    energyout_filter = openmc.EnergyoutFilter(energy_bins)
    cell_filter = openmc.CellFilter([cell])
    material_filter = openmc.MaterialFilter([u238])

    # Define tally
    tally = openmc.Tally(name='photon spectrum')
    tally.filters = [particle_filter, particleout_filter, energyout_filter]
    tally.scores = ['total']
    tally.estimator = 'analog'
    model.tallies = openmc.Tallies([tally])
    return model


def test_filter_particleout(model):
    harness = PyAPITestHarness("statepoint.10.h5", model=model)
    harness.main()
