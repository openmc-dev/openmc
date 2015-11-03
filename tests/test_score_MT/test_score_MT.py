#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import TestHarness, PyAPITestHarness
import openmc


class ScoreMTTestHarness(PyAPITestHarness):
    def _build_inputs(self):
        filt = openmc.Filter(type='cell', bins=(10, 21, 22, 23))
        tallies = [openmc.Tally(tally_id=i) for i in range(1, 4)]
        [t.add_filter(filt) for t in tallies]
        [t.add_score('n2n') for t in tallies]
        [t.add_score('16') for t in tallies]
        [t.add_score('51') for t in tallies]
        [t.add_score('102') for t in tallies]
        tallies[0].estimator = 'tracklength'
        tallies[1].estimator = 'analog'
        tallies[2].estimator = 'collision'
        self._input_set.tallies = openmc.TalliesFile()
        [self._input_set.tallies.add_tally(t) for t in tallies]

        super(ScoreMTTestHarness, self)._build_inputs()

    def _cleanup(self):
        super(ScoreMTTestHarness, self)._cleanup()
        f = os.path.join(os.getcwd(), 'tallies.xml')
        if os.path.exists(f): os.remove(f)


if __name__ == '__main__':
    harness = ScoreMTTestHarness('statepoint.10.*', True)
    harness.main()
