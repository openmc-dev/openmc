#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, os.pardir))
from testing_harness import TestHarness, PyAPITestHarness
import openmc


class ScoreNuScatterTestHarness(PyAPITestHarness):
    def _build_inputs(self):
        filt = openmc.Filter(type='cell', bins=(10, 21, 22, 23))
        t = openmc.Tally(tally_id=1)
        t.add_filter(filt)
        t.add_score('nu-scatter')
        self._input_set.tallies = openmc.TalliesFile()
        self._input_set.tallies.add_tally(t)

        PyAPITestHarness._build_inputs(self)

    def _cleanup(self):
        PyAPITestHarness._cleanup(self)
        f = os.path.join(os.getcwd(), 'tallies.xml')
        if os.path.exists(f): os.remove(f)


if __name__ == '__main__':
    harness = ScoreNuScatterTestHarness('statepoint.10.*', True)
    harness.main()
