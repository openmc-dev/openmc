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

    energy_bins = np.array(
      [1.00000000e+04, 1.09749877e+04, 1.20450354e+04, 1.32194115e+04,
       1.45082878e+04, 1.59228279e+04, 1.74752840e+04, 1.91791026e+04,
       2.10490414e+04, 2.31012970e+04, 2.53536449e+04, 2.78255940e+04,
       3.05385551e+04, 3.35160265e+04, 3.67837977e+04, 4.03701726e+04,
       4.43062146e+04, 4.86260158e+04, 5.33669923e+04, 5.85702082e+04,
       6.42807312e+04, 7.05480231e+04, 7.74263683e+04, 8.49753436e+04,
       9.32603347e+04, 1.02353102e+05, 1.12332403e+05, 1.23284674e+05,
       1.35304777e+05, 1.48496826e+05, 1.62975083e+05, 1.78864953e+05,
       1.96304065e+05, 2.15443469e+05, 2.36448941e+05, 2.59502421e+05,
       2.84803587e+05, 3.12571585e+05, 3.43046929e+05, 3.76493581e+05,
       4.13201240e+05, 4.53487851e+05, 4.97702356e+05, 5.46227722e+05,
       5.99484250e+05, 6.57933225e+05, 7.22080902e+05, 7.92482898e+05,
       8.69749003e+05, 9.54548457e+05, 1.04761575e+06, 1.14975700e+06,
       1.26185688e+06, 1.38488637e+06, 1.51991108e+06, 1.66810054e+06,
       1.83073828e+06, 2.00923300e+06, 2.20513074e+06, 2.42012826e+06,
       2.65608778e+06, 2.91505306e+06, 3.19926714e+06, 3.51119173e+06,
       3.85352859e+06, 4.22924287e+06, 4.64158883e+06, 5.09413801e+06,
       5.59081018e+06, 6.13590727e+06, 6.73415066e+06, 7.39072203e+06,
       8.11130831e+06, 8.90215085e+06, 9.77009957e+06, 1.07226722e+07,
       1.17681195e+07, 1.29154967e+07, 1.41747416e+07, 1.55567614e+07,
       1.70735265e+07, 1.87381742e+07, 2.05651231e+07, 2.25701972e+07,
       2.47707636e+07, 2.71858824e+07, 2.98364724e+07, 3.27454916e+07,
       3.59381366e+07, 3.94420606e+07, 4.32876128e+07, 4.75081016e+07,
       5.21400829e+07, 5.72236766e+07, 6.28029144e+07, 6.89261210e+07,
       7.56463328e+07, 8.30217568e+07, 9.11162756e+07, 1.00000000e+08])

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
