"""Test the settings.no_reduce feature to ensure tallies are correctly
reduced across MPI processes."""

import openmc
import pytest

from tests.testing_harness import config


@pytest.mark.parametrize('no_reduce', [True, False])
def test_no_reduce(no_reduce, run_in_tmpdir):
    """Test that tally results are correct with and without no_reduce."""

    # Create simple sphere model with vacuum
    model = openmc.Model()
    sphere = openmc.Sphere(r=1.0, boundary_type='vacuum')
    cell = openmc.Cell(region=-sphere)
    model.geometry = openmc.Geometry([cell])
    model.settings.run_mode = 'fixed source'
    model.settings.batches = 10
    model.settings.particles = 100
    model.settings.source = openmc.IndependentSource(space=openmc.stats.Point())
    model.settings.no_reduce = no_reduce

    # Tally: surface current on vacuum boundary
    surf_filter = openmc.SurfaceFilter(sphere)
    tally = openmc.Tally()
    tally.filters = [surf_filter]
    tally.scores = ['current']
    model.tallies = [tally]

    # Run OpenMC with proper MPI arguments if needed
    kwargs = {'apply_tally_results': True, 'openmc_exec': config['exe']}
    if config['mpi']:
        kwargs['mpi_args'] = [config['mpiexec'], '-n', config['mpi_np']]
    model.run(**kwargs)

    # The tally should be ~1.0 (every particle crosses the surface once)
    tally_mean = tally.mean.flatten()[0]
    assert tally_mean == pytest.approx(1.0)
