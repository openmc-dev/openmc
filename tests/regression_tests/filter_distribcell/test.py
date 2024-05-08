import glob
import os
from pathlib import Path

import pytest

from tests.testing_harness import *


class DistribcellTestHarness(TestHarness):
    def __init__(self, statepoint_name, dir, tallies_out, hash_out):
        self.dir = dir
        self.tallies_out = tallies_out
        self.hash_out = hash_out
        super().__init__(statepoint_name)

    def execute_test(self):
        """Run OpenMC with the appropriate arguments and check the outputs."""
        try:
            dirs = ('case-1', '../case-2', '../case-3', '../case-4')
            sps = ('statepoint.1.*', 'statepoint.1.*', 'statepoint.3.*',
                   'statepoint.1.*')
            tallies_out_present = (True, True, False, True)
            hash_out = (False, False, True, False)

            os.chdir(self.dir)
            self._run_openmc()
            self._test_output_created(self.tallies_out)
            results = self._get_results(self.hash_out)
            self._write_results(results)
            self._compare_results()
        finally:
            os.chdir(self.dir)
            self._cleanup()

    def update_results(self):
        """Update the results_true using the current version of OpenMC."""
        try:
            dirs = ('case-1', '../case-2', '../case-3', '../case-4', '../case-5')
            sps = ('statepoint.1.h5', 'statepoint.1.h5', 'statepoint.3.h5',
                   'statepoint.1.h5', 'statepoint.1.h5')
            tallies_out_present = (True, True, False, True, True)
            hash_out = (False, False, True, False, False)

            os.chdir(self.dir)
            self._run_openmc()
            self._test_output_created(self.tallies_out)
            results = self._get_results(self.hash_out)
            self._write_results(results)
            self._overwrite_results()
        finally:
            os.chdir(self.dir)
            self._cleanup()

    def _test_output_created(self, tallies_out_present):
        """Make sure statepoint.* and tallies.out have been created."""
        statepoint = glob.glob(os.path.join(os.getcwd(), self._sp_name))
        assert len(statepoint) == 1, 'Either multiple or no statepoint files ' \
             'exist.'
        assert statepoint[0].endswith('h5'), \
             'Statepoint file is not a HDF5 file.'
        if tallies_out_present:
            assert os.path.exists(os.path.join(os.getcwd(), 'tallies.out')), \
                 'Tally output file does not exist.'


cases = [
    ('case-1', 'statepoint.1.h5', True, False),
    ('case-2', 'statepoint.1.h5', True, False),
    ('case-3', 'statepoint.3.h5', False, True),
    ('case-4', 'statepoint.1.h5', True, False),
    ('case-5', 'statepoint.1.h5', True, False)
]
@pytest.mark.parametrize(','.join(['dir', 'statepoint', 'tallies_out', 'hash_out']), cases,
                         ids=[cases[0] for cases in cases])
def test_filter_distribcell(request, dir, statepoint, tallies_out, hash_out):
    dir = request.fspath.dirname / Path(dir)
    harness = DistribcellTestHarness(statepoint, dir.resolve(), tallies_out, hash_out)
    harness.main()
