#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, os.pardir))
from testing_harness import TestHarness, PyAPITestHarness
import openmc


class ScoreScatterYNTestHarness(PyAPITestHarness):
    def _build_inputs(self):
        filt = openmc.Filter(type='cell', bins=(21, ))
        t1 = openmc.Tally(tally_id=1)
        t1.add_filter(filt)
        t1.add_score('scatter-0')
        t1.estimator = 'analog'
        t2 = openmc.Tally(tally_id=2)
        t2.add_filter(filt)
        t2.add_score('scatter-y4')
        self._input_set.tallies = openmc.TalliesFile()
        self._input_set.tallies.add_tally(t1)
        self._input_set.tallies.add_tally(t2)

        PyAPITestHarness._build_inputs(self)

    def _cleanup(self):
        PyAPITestHarness._cleanup(self)
        f = os.path.join(os.getcwd(), 'tallies.xml')
        if os.path.exists(f): os.remove(f)


if __name__ == '__main__':
    harness = ScoreScatterYNTestHarness('statepoint.10.*', True)
    harness.main()
