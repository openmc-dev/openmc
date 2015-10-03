#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import TestHarness, PyAPITestHarness
import openmc


class ScoreTotalYNTestHarness(PyAPITestHarness):
    def _build_inputs(self):
        filt = openmc.Filter(type='cell', bins=(10, 21, 22, 23))
        tallies = [openmc.Tally(tally_id=i) for i in range(1, 5)]
        [t.add_filter(filt) for t in tallies]
        tallies[0].add_score('total')
        [t.add_score('total-y4') for t in tallies[1:]]
        [t.add_nuclide('U-235') for t in tallies[1:]]
        [t.add_nuclide('total') for t in tallies[1:]]
        tallies[1].estimator = 'tracklength'
        tallies[2].estimator = 'analog'
        tallies[3].estimator = 'collision'
        self._input_set.tallies = openmc.TalliesFile()
        [self._input_set.tallies.add_tally(t) for t in tallies]

        super(ScoreTotalYNTestHarness, self)._build_inputs()

    def _cleanup(self):
        super(ScoreTotalYNTestHarness, self)._cleanup()
        f = os.path.join(os.getcwd(), 'tallies.xml')
        if os.path.exists(f): os.remove(f)


if __name__ == '__main__':
    harness = ScoreTotalYNTestHarness('statepoint.10.*', True)
    harness.main()
