#!/usr/bin/env python

import sys
sys.path.insert(0, '..')
from testing_harness import TestHarness, PyAPITestHarness
import openmc
import os


class ScoreNuScatterNTestHarness(PyAPITestHarness):
    def _build_inputs(self):
        filt = openmc.Filter(type='cell', bins=(21, 22, 23))
        t = openmc.Tally(tally_id=1)
        t.add_filter(filt)
        t.add_score('nu-scatter')
        t.add_score('nu-scatter-1')
        t.add_score('nu-scatter-2')
        t.add_score('nu-scatter-3')
        t.add_score('nu-scatter-4')
        self._input_set.tallies = openmc.TalliesFile()
        self._input_set.tallies.add_tally(t)

        PyAPITestHarness._build_inputs(self)

    def _cleanup(self):
        PyAPITestHarness._cleanup(self)
        f = os.path.join(os.getcwd(), 'tallies.xml')
        if os.path.exists(f): os.remove(f)


if __name__ == '__main__':
    harness = ScoreNuScatterNTestHarness('statepoint.10.*', True)
    harness.main()
