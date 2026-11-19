import openmc
import pytest
from openmc.utility_funcs import change_directory

from tests.testing_harness import PyAPITestHarness


@pytest.mark.parametrize("shared_secondary,subdir", [
    (False, "local"),
    (True, "shared"),
])
def test_particle_production_fission(shared_secondary, subdir):
    """Fixed-source model with fissionable material to test that
    ParticleProductionFilter correctly counts fission-born neutrons,
    with both local and shared secondary bank modes."""
    with change_directory(subdir):
        openmc.reset_auto_ids()
        model = openmc.Model()

        mat = openmc.Material()
        mat.set_density('g/cm3', 18.0)
        mat.add_nuclide('U235', 1.0)
        model.materials.append(mat)

        sph = openmc.Sphere(r=5.0, boundary_type='vacuum')
        cell = openmc.Cell(fill=mat, region=-sph)
        model.geometry = openmc.Geometry([cell])

        source = openmc.IndependentSource()
        source.energy = openmc.stats.delta_function(1.0e6)

        model.settings.particles = 100
        model.settings.run_mode = 'fixed source'
        model.settings.batches = 2
        model.settings.source = source
        model.settings.create_fission_neutrons = True
        model.settings.shared_secondary_bank = shared_secondary

        # ParticleProductionFilter tracking fission neutron production
        ppf = openmc.ParticleProductionFilter(['neutron'])
        neutron_filter = openmc.ParticleFilter(['neutron'])
        tally = openmc.Tally()
        tally.filters = [neutron_filter, ppf]
        tally.scores = ['events']
        tally.estimator = 'analog'
        model.tallies = [tally]

        harness = PyAPITestHarness('statepoint.2.h5', model)
        harness.main()
