import numpy as np
import openmc
import pytest

from tests.mg_photon import (
    ABSORPTION, GROUP_EDGES, TOTAL, check_exact_solution, make_model,
)
from tests.regression_tests import config


def test_mg_photon_survival_biasing(run_in_tmpdir):
    """Survival biasing must leave the exact answers unchanged"""
    model = make_model(survival_biasing=True)
    statepoint = model.run(
        openmc_exec=config['exe'], event_based=config['event'])
    check_exact_solution(statepoint)


def test_mg_photon_nuclide_scatter(run_in_tmpdir):
    """Check that nuclide and material tallies use the same derived multiplicity."""
    model = make_model('derived')
    tally = openmc.Tally(name='nuclide scatter')
    tally.filters = [openmc.EnergyFilter(GROUP_EDGES)]
    tally.nuclides = ['photon', 'total']
    tally.scores = ['scatter']
    model.tallies.append(tally)
    statepoint = model.run(
        openmc_exec=config['exe'], event_based=config['event'])
    with openmc.StatePoint(statepoint) as sp:
        tally = sp.get_tally(name='nuclide scatter')
        library = tally.get_values(nuclides=['photon']).ravel()
        material = tally.get_values(nuclides=['total']).ravel()
    assert np.allclose(library, material), \
        f'Nuclide scatter {library} differs from material {material}'


@pytest.mark.parametrize(
    ('invalid_input', 'error'),
    [
        ('energy above library', 'Source energy above range'),
        ('energy below library', 'outside the MGXS group structure'),
        ('source mismatch', 'does not match'),
        ('pulse-height tally', 'Pulse-height tallies are not supported'),
        ('no scattering', 'produces photons but has no scattering'),
    ],
)
def test_mg_photon_validation(run_in_tmpdir, invalid_input, error):
    if invalid_input == 'no scattering':
        # A group that only absorbs leaves nothing to carry its photon
        # production when the multiplicity has to be derived
        absorption = ABSORPTION.copy()
        absorption[-1] = TOTAL[-1]
        model = make_model('derived', absorption=absorption)
    else:
        model = make_model()

    source = model.settings.source[0]
    if invalid_input == 'energy above library':
        source.energy = openmc.stats.Discrete([3.0e6], [1.0])
    elif invalid_input == 'energy below library':
        source.energy = openmc.stats.Discrete([500.0], [1.0])
    elif invalid_input == 'source mismatch':
        source.particle = 'neutron'
    elif invalid_input == 'pulse-height tally':
        cells = list(model.geometry.get_all_cells().values())
        tally = openmc.Tally()
        tally.filters = [
            openmc.CellFilter(cells),
            openmc.EnergyFilter(GROUP_EDGES),
        ]
        tally.scores = ['pulse-height']
        model.tallies.append(tally)

    with pytest.raises(RuntimeError, match=error):
        model.run(openmc_exec=config['exe'])
