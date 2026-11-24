import numpy as np
import openmc
import pytest

from tests.mg_photon import (
    ABSORPTION, GROUP_EDGES, TOTAL, check_exact_solution, make_model,
)
from tests.regression_tests import config


def test_mg_photon_survival_biasing(run_in_tmpdir):
    """Survival biasing must leave the exact answers unchanged"""
    model = make_model()
    model.settings.survival_biasing = True
    statepoint = model.run(
        openmc_exec=config['exe'], event_based=config['event'])
    check_exact_solution(statepoint)


def test_mg_photon_nuclide_scatter(run_in_tmpdir):
    """Check that nuclide tallies use the derived multiplicity."""
    model = make_model('derived')
    tally = openmc.Tally(name='nuclide scatter')
    tally.filters = [openmc.EnergyFilter(GROUP_EDGES)]
    tally.nuclides = ['photon', 'total']
    tally.scores = ['scatter']
    model.tallies.append(tally)
    model.run(openmc_exec=config['exe'], event_based=config['event'],
              apply_tally_results=True)
    np.testing.assert_allclose(
        tally.get_values(nuclides=['photon']),
        tally.get_values(nuclides=['total']),
        rtol=1e-5, atol=1e-8, equal_nan=False)


@pytest.mark.parametrize(
    ('particle', 'energy', 'error'),
    [
        ('photon', 3.0e6, 'Source energy above range'),
        ('photon', 500.0, 'outside the MGXS group structure'),
        ('neutron', 1.0e6, 'does not match'),
    ],
)
def test_mg_photon_invalid_source(run_in_tmpdir, particle, energy, error):
    model = make_model()
    source = model.settings.source[0]
    source.particle = particle
    source.energy = openmc.stats.delta_function(energy)
    with pytest.raises(RuntimeError, match=error):
        model.run(openmc_exec=config['exe'])


def test_mg_photon_pulse_height(run_in_tmpdir):
    model = make_model()
    cells = list(model.geometry.get_all_cells().values())
    tally = openmc.Tally()
    tally.filters = [openmc.CellFilter(cells),
                     openmc.EnergyFilter(GROUP_EDGES)]
    tally.scores = ['pulse-height']
    model.tallies.append(tally)
    with pytest.raises(
        RuntimeError, match='Pulse-height tallies are not supported'
    ):
        model.run(openmc_exec=config['exe'])


def test_mg_photon_production_without_scattering(run_in_tmpdir):
    absorption = ABSORPTION.copy()
    absorption[-1] = TOTAL[-1]
    model = make_model('derived', absorption=absorption)
    with pytest.raises(
        RuntimeError, match='produces photons but has no scattering'
    ):
        model.run(openmc_exec=config['exe'])


def test_mg_photon_scattering_without_production(run_in_tmpdir):
    model = make_model('derived')
    library = openmc.MGXSLibrary.from_hdf5('mgxs.h5')
    xsdata = library.xsdatas[0]
    scatter_matrix = xsdata.scatter_matrix[0].copy()
    scatter_matrix[0] = 0.0
    xsdata.set_scatter_matrix(scatter_matrix)
    library.export_to_hdf5('mgxs.h5')
    with pytest.raises(RuntimeError, match='but no photon production'):
        model.run(openmc_exec=config['exe'])


def test_mg_photon_random_ray(run_in_tmpdir):
    model = make_model()
    model.settings.source[0].space = openmc.stats.Point()
    model.settings.random_ray = {
        'distance_inactive': 20.0,
        'distance_active': 100.0,
        'ray_source': openmc.IndependentSource(
            space=openmc.stats.Box((-5.0,) * 3, (5.0,) * 3)),
    }
    model.tallies = []
    with pytest.raises(RuntimeError, match='not supported by the random ray'):
        model.run(openmc_exec=config['exe'])
