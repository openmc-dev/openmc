#!/usr/bin/env python

import os
import sys
import glob
import hashlib
sys.path.insert(0, os.pardir)
from testing_harness import PyAPITestHarness
import openmc


class TallyAggregationTestHarness(PyAPITestHarness):
    def _build_inputs(self):

        # The summary.h5 file needs to be created to read in the tallies
        self._input_set.settings.output = {'summary': True}

        # Initialize the nuclides
        u235 = openmc.Nuclide('U-235')
        u238 = openmc.Nuclide('U-238')
        pu239 = openmc.Nuclide('Pu-239')

        # Initialize the filters
        energy_filter = openmc.Filter(type='energy', bins=[0.0, 0.253e-6,
                                                           1.0e-3, 1.0, 20.0])
        distrib_filter   = openmc.Filter(type='distribcell', bins=[60])

        # Initialized the tallies
        tally = openmc.Tally(name='distribcell tally')
        tally.filters = [energy_filter, distrib_filter]
        tally.scores = ['nu-fission', 'total']
        tally.nuclides = [u235, u238, pu239]
        tallies_file = openmc.Tallies([tally])

        # Export tallies to file
        self._input_set.tallies = tallies_file
        super(TallyAggregationTestHarness, self)._build_inputs()

    def _get_results(self, hash_output=True):
        """Digest info in the statepoint and return as a string."""

        # Read the statepoint file.
        statepoint = glob.glob(os.path.join(os.getcwd(), self._sp_name))[0]
        sp = openmc.StatePoint(statepoint)

        # Extract the tally of interest
        tally = sp.get_tally(name='distribcell tally')

        # Perform tally aggregations across filter bins, nuclides and scores
        outstr = ''

        # Sum across all energy filter bins
        tally_sum = tally.summation(filter_type='energy')
        outstr += ', '.join(map(str, tally_sum.mean))
        outstr += ', '.join(map(str, tally_sum.std_dev))

        # Sum across all distribcell filter bins
        tally_sum = tally.summation(filter_type='distribcell')
        outstr += ', '.join(map(str, tally_sum.mean))
        outstr += ', '.join(map(str, tally_sum.std_dev))

        # Sum across all nuclides
        tally_sum = tally.summation(nuclides=['U-235', 'U-238', 'Pu-239'])
        outstr += ', '.join(map(str, tally_sum.mean))
        outstr += ', '.join(map(str, tally_sum.std_dev))

        # Sum across all scores
        tally_sum = tally.summation(scores=['nu-fission', 'total'])
        outstr += ', '.join(map(str, tally_sum.mean))
        outstr += ', '.join(map(str, tally_sum.std_dev))

        # Hash the results if necessary
        if hash_output:
            sha512 = hashlib.sha512()
            sha512.update(outstr.encode('utf-8'))
            outstr = sha512.hexdigest()

        return outstr

    def _cleanup(self):
        super(TallyAggregationTestHarness, self)._cleanup()
        f = os.path.join(os.getcwd(), 'tallies.xml')
        if os.path.exists(f): os.remove(f)

if __name__ == '__main__':
    harness = TallyAggregationTestHarness('statepoint.10.h5', True)
    harness.main()
