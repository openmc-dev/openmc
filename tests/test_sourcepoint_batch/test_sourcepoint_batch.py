#!/usr/bin/env python

import glob
import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import TestHarness
from openmc import StatePoint


class SourcepointTestHarness(TestHarness):
    def _test_output_created(self):
        """Make sure statepoint.* files have been created."""
        statepoint = glob.glob(os.path.join(os.getcwd(), 'statepoint.*'))
        assert len(statepoint) == 5, '5 statepoint files must exist.'
        assert statepoint[0].endswith('h5'), \
             'Statepoint file is not a HDF5 file.'

    def _get_results(self):
        """Digest info in the statepoint and return as a string."""
        # Get the eigenvalue information.
        outstr = TestHarness._get_results(self)

        # Read the statepoint file.
        statepoint = glob.glob(os.path.join(os.getcwd(), self._sp_name))[0]
        with StatePoint(statepoint) as sp:
            # Add the source information.
            xyz = sp.source[0]['xyz']
            outstr += ' '.join(['{0:12.6E}'.format(x) for x in xyz])
            outstr += "\n"

        return outstr


if __name__ == '__main__':
    harness = SourcepointTestHarness('statepoint.08.h5')
    harness.main()
