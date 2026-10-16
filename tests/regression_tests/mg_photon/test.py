"""Multigroup photon transport with the Monte Carlo solver.

A reflective cube with a spatially uniform, isotropic source behaves as an
infinite medium, so the group track lengths (volume-integrated fluxes) per
source photon follow from a group balance and can be checked exactly. The
library uses P1 scattering and a top group that produces more photons than it
scatters (a row-constant multiplicity), which exercises angular sampling and
the weight multiplicity used to represent secondary photon production.

"""
from pathlib import Path

import numpy as np
import openmc
import pytest

from openmc.utility_funcs import change_directory
from tests.regression_tests import config
from tests.testing_harness import PyAPITestHarness

GROUP_EDGES = [1.0e3, 1.0e5, 5.0e5, 2.0e6]

# Macroscopic cross sections in 1/cm; index 0 is the highest-energy group
TOTAL = np.array([0.25, 0.40, 1.20])
ABSORPTION = np.array([0.05, 0.15, 1.00])

# P0 production (nu-scatter) matrix indexed [g_in, g_out]. The top group
# produces 0.25/cm but only scatters TOTAL - ABSORPTION = 0.20/cm; the extra
# photons stand in for fluorescence and annihilation photons.
PRODUCTION = np.array([
    [0.10, 0.08, 0.07],
    [0.00, 0.15, 0.10],
    [0.00, 0.00, 0.20],
])
MU_BAR = 0.2  # average scattering cosine, P1/P0

SPEED_OF_LIGHT = 2.99792458e10  # cm/s


def _multiplicity():
    """Row-constant multiplicity so each collision yields PRODUCTION/TOTAL"""
    ratio = PRODUCTION.sum(axis=1) / (TOTAL - ABSORPTION)
    return np.repeat(ratio[:, np.newaxis], len(TOTAL), axis=1)


def _exact_track_lengths():
    """Group track length per source photon in an infinite medium.

    Each group balances collisions against source and in-production:
    TOTAL[g] * L[g] = S[g] + sum_g' PRODUCTION[g', g] * L[g']

    """
    source = np.array([1.0, 0.0, 0.0])  # 1 MeV photons are born in group 1
    return np.linalg.solve(np.diag(TOTAL) - PRODUCTION.T, source)


def _make_model():
    groups = openmc.mgxs.EnergyGroups(group_edges=GROUP_EDGES)
    photon = openmc.XSdata('photon', groups)
    photon.order = 1
    photon.set_total(TOTAL)
    photon.set_absorption(ABSORPTION)
    photon.set_scatter_matrix(
        np.stack([PRODUCTION, MU_BAR * PRODUCTION], axis=-1))
    photon.set_multiplicity_matrix(_multiplicity())
    # No inverse velocity is given, so photons must default to 1/c

    library = openmc.MGXSLibrary(groups, particle_type='photon')
    library.add_xsdata(photon)
    library.export_to_hdf5('mgxs.h5')

    material = openmc.Material()
    material.set_density('macro', 1.0)
    material.add_macroscopic('photon')

    box = openmc.model.RectangularParallelepiped(
        -5.0, 5.0, -5.0, 5.0, -5.0, 5.0, boundary_type='reflective')
    cell = openmc.Cell(fill=material, region=-box)

    model = openmc.Model()
    model.geometry = openmc.Geometry([cell])
    model.materials = openmc.Materials([material])
    model.materials.cross_sections = 'mgxs.h5'

    tally = openmc.Tally(name='photon')
    tally.filters = [
        openmc.EnergyFilter(GROUP_EDGES),
        openmc.ParticleFilter('photon'),
    ]
    tally.scores = ['flux', 'total', 'absorption', 'scatter', 'nu-scatter',
                    'inverse-velocity']
    model.tallies.append(tally)

    model.settings.energy_mode = 'multi-group'
    model.settings.run_mode = 'fixed source'
    model.settings.photon_transport = True
    model.settings.batches = 20
    model.settings.particles = 1000
    model.settings.source = openmc.IndependentSource(
        particle='photon',
        space=openmc.stats.Box((-5.0, -5.0, -5.0), (5.0, 5.0, 5.0)),
        energy=openmc.stats.Discrete([1.0e6], [1.0]))

    return model


@pytest.fixture
def model():
    return _make_model()


class MGPhotonTestHarness(PyAPITestHarness):
    """Check exact infinite-medium answers before the regression comparison"""

    def _get_results(self, hash_output=False):
        exact = _exact_track_lengths()
        expected = {
            'flux': exact,
            'total': TOTAL * exact,
            'absorption': ABSORPTION * exact,
            'scatter': (TOTAL - ABSORPTION) * exact,
            'nu-scatter': PRODUCTION.sum(axis=1) * exact,
        }

        with openmc.StatePoint(self._sp_name) as statepoint:
            tally = statepoint.get_tally(name='photon')

            # Energy filter bins run from low to high energy, so reverse them
            # to match the group ordering of the library
            for score, reference in expected.items():
                mean = tally.get_values(scores=[score]).ravel()[::-1]
                std_dev = tally.get_values(
                    scores=[score], value='std_dev').ravel()[::-1]
                deviation = np.abs(mean - reference) / std_dev
                if np.any(deviation > 5.0):
                    raise AssertionError(
                        f"Multigroup photon '{score}' disagrees with the "
                        f"infinite-medium solution: {mean} vs. {reference} "
                        f"({deviation} standard deviations)")

            # Photons travel at the speed of light, so the inverse-velocity
            # score must equal flux / c in every group
            flux = tally.get_values(scores=['flux']).ravel()
            inverse_velocity = tally.get_values(
                scores=['inverse-velocity']).ravel()
            if not np.allclose(SPEED_OF_LIGHT * inverse_velocity, flux,
                               rtol=1e-10):
                raise AssertionError(
                    'Multigroup photons are not moving at the speed of light: '
                    f'flux / inverse-velocity = {flux / inverse_velocity}')

        return super()._get_results(hash_output)

    def _cleanup(self):
        super()._cleanup()
        Path('mgxs.h5').unlink(missing_ok=True)


def test_mg_photon(model):
    harness = MGPhotonTestHarness('statepoint.20.h5', model)
    harness.main()


@pytest.mark.parametrize(
    ('invalid_input', 'error'),
    [
        ('energy above library', 'Source energy above range'),
        ('energy below library', 'outside the MGXS group structure'),
        ('source mismatch', 'does not match'),
        ('pulse-height tally', 'Pulse-height tallies are not supported'),
    ],
)
def test_mg_photon_validation(tmp_path, invalid_input, error):
    with change_directory(tmp_path):
        model = _make_model()
        source = model.settings.source[0]
        if invalid_input == 'energy above library':
            source.energy = openmc.stats.Discrete([3.0e6], [1.0])
        elif invalid_input == 'energy below library':
            source.energy = openmc.stats.Discrete([500.0], [1.0])
        elif invalid_input == 'source mismatch':
            source.particle = 'neutron'
        else:
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
