#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import TestHarness, PyAPITestHarness
import openmc


class FilterDelayedgroupTestHarness(PyAPITestHarness):
    def _build_inputs(self):
        filt = openmc.Filter(type='delayedgroup',
                             bins=(1, 2, 3, 4, 5, 6))
        tally = openmc.Tally(tally_id=1)
        tally.add_filter(filt)
        tally.add_score('delayed-nu-fission')
        self._input_set.tallies = openmc.TalliesFile()
        self._input_set.tallies.add_tally(tally)

        super(FilterDelayedgroupTestHarness, self)._build_inputs()

    def _cleanup(self):
        super(FilterDelayedgroupTestHarness, self)._cleanup()
        f = os.path.join(os.getcwd(), 'tallies.xml')
        if os.path.exists(f): os.remove(f)


if __name__ == '__main__':
    harness = FilterDelayedgroupTestHarness('statepoint.10.*', True)
    harness.main()
