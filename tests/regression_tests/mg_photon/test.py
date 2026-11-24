from pathlib import Path

import openmc
import pytest

from tests.mg_photon import check_exact_solution, make_model
from tests.testing_harness import PyAPITestHarness


class MGPhotonTestHarness(PyAPITestHarness):
    """Check exact infinite-medium answers before the regression comparison"""

    def _get_results(self, hash_output=False):
        check_exact_solution(self._sp_name)
        return super()._get_results(hash_output)

    def _cleanup(self):
        super()._cleanup()
        Path('mgxs.h5').unlink(missing_ok=True)


# The derived multiplicity equals the explicit row-constant one, and without
# survival biasing the multiplicity only scales particle weights, so both
# libraries must reproduce the same inputs_true.dat and results_true.dat.
@pytest.mark.parametrize('multiplicity', ['explicit', 'derived'])
def test_mg_photon(multiplicity):
    openmc.reset_auto_ids()
    harness = MGPhotonTestHarness(
        'statepoint.20.h5', make_model(multiplicity))
    harness.main()
