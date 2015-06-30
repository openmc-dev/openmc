#!/usr/bin/env python

import sys
sys.path.insert(0, '..')
from testing_harness import *


class SourcepointTestHarness(TestHarness):
    def _test_output_created(self):
        """Make sure statepoint.* files have been created."""
        statepoint = glob.glob(os.path.join(os.getcwd(), 'statepoint.*'))
        assert len(statepoint) == 5, '5 statepoint files must exist.' 
        assert statepoint[0].endswith('binary') \
             or statepoint[0].endswith('h5'), \
             'Statepoint file is not a binary or hdf5 file.'

    def _get_results(self):
        """Digest info in the statepoint and return as a string."""
        # Read the statepoint file.
        statepoint = glob.glob(os.path.join(os.getcwd(), self._sp_name))[0]
        sp = StatePoint(statepoint)
        sp.read_results()
        sp.read_source()

        # Get the eigenvalue information.
        outstr = TestHarness._get_results(self)

        # Add the source information.
        xyz = sp._source[0]._xyz
        outstr += ' '.join(['{0:12.6E}'.format(x) for x in xyz])
        outstr += "\n"

        return outstr


if __name__ == '__main__':
    harness = SourcepointTestHarness('statepoint.08.*')
    harness.main()
