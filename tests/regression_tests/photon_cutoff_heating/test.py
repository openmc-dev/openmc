from pathlib import Path

import openmc
import pytest
from openmc.utility_funcs import change_directory

from tests.testing_harness import PyAPITestHarness
from tests.unit_tests.test_energy_cutoff import inf_medium_model


SOURCE_ENERGY = 1.0e6
PHOTON_CUTOFF = 0.25e6


class PhotonCutoffHeatingTestHarness(PyAPITestHarness):
    def _get_results(self):
        with openmc.StatePoint(self._sp_name) as sp:
            heating = sp.get_tally(name='heating').mean.ravel()
            # The reflective medium absorbs all source energy. Both total
            # and nuclide-specific heating include cutoff deposition.
            assert heating == pytest.approx([SOURCE_ENERGY, SOURCE_ENERGY])

            pulse_height = sp.get_tally(name='pulse-height').mean.ravel()
            assert pulse_height == pytest.approx([0.0, 1.0])

            outgoing = sp.get_tally(name='outgoing').mean.ravel()
            # Retain the physical energy of photons killed by the cutoff
            # for outgoing-energy filters; this bin excludes zero energy.
            assert outgoing[1] > 0.0

        return super()._get_results()


@pytest.mark.parametrize('electron_treatment', ['led', 'ttb'])
@pytest.mark.parametrize('shared_secondary', [False, True])
def test_photon_cutoff_heating(electron_treatment, shared_secondary):
    bank_type = 'shared' if shared_secondary else 'local'
    subdir = Path(f'{electron_treatment}_{bank_type}')
    subdir.mkdir(exist_ok=True)
    with change_directory(subdir):
        openmc.reset_auto_ids()
        # Reuse the existing photon-cutoff model. A source below the pair
        # production threshold isolates the scattered-photon cutoff.
        model = inf_medium_model(PHOTON_CUTOFF, SOURCE_ENERGY)
        model.settings.electron_treatment = electron_treatment
        model.settings.shared_secondary_bank = shared_secondary

        heating = model.tallies[1]
        heating.nuclides = ['total', 'Zr90']

        pulse_height = openmc.Tally(name='pulse-height')
        pulse_height.filters = [
            openmc.CellFilter(list(model.geometry.get_all_cells().values())),
            openmc.EnergyFilter(
                [0.0, SOURCE_ENERGY - 1.0, SOURCE_ENERGY + 1.0]),
        ]
        pulse_height.scores = ['pulse-height']

        outgoing = openmc.Tally(name='outgoing')
        outgoing.filters = [
            openmc.EnergyoutFilter([0.0, 1.0, PHOTON_CUTOFF, SOURCE_ENERGY]),
            openmc.ParticleFilter('photon'),
        ]
        outgoing.scores = ['events']
        outgoing.estimator = 'analog'
        model.tallies.extend([pulse_height, outgoing])

        harness = PhotonCutoffHeatingTestHarness('statepoint.10.h5', model)
        harness.main()
