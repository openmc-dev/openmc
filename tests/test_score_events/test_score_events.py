#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import TestHarness, PyAPITestHarness
import openmc


class ScoreEventsTestHarness(PyAPITestHarness):
    def _build_inputs(self):
        filt = openmc.Filter(type='cell', bins=(21, 27))
        tallies = [openmc.Tally(tally_id=i) for i in range(1, 3)]
        [t.add_filter(filt) for t in tallies]
        [t.add_score('events') for t in tallies]
        tallies[0].estimator = 'tracklength'
        tallies[1].estimator = 'analog'
        self._input_set.tallies = openmc.TalliesFile()
        [self._input_set.tallies.add_tally(t) for t in tallies]

        super(ScoreEventsTestHarness, self)._build_inputs()

    def _cleanup(self):
        super(ScoreEventsTestHarness, self)._cleanup()
        f = os.path.join(os.getcwd(), 'tallies.xml')
        if os.path.exists(f): os.remove(f)


if __name__ == '__main__':
    harness = ScoreEventsTestHarness('statepoint.10.*', True)
    harness.main()
