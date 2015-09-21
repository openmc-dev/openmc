#!/usr/bin/env python

import sys
sys.path.insert(0, '..')
from testing_harness import TestHarness, PyAPITestHarness
import openmc
import os


class FilterCellbornTestHarness(PyAPITestHarness):
    def _build_inputs(self):
        filt = openmc.Filter(type='cellborn', bins=(10, 21, 22, 23))
        tally = openmc.Tally(tally_id=1)
        tally.add_filter(filt)
        tally.add_score('total')
        self._input_set.tallies = openmc.TalliesFile()
        self._input_set.tallies.add_tally(tally)

        PyAPITestHarness._build_inputs(self)

    def _cleanup(self):
        PyAPITestHarness._cleanup(self)
        f = os.path.join(os.getcwd(), 'tallies.xml')
        if os.path.exists(f): os.remove(f)


if __name__ == '__main__':
    harness = FilterCellbornTestHarness('statepoint.10.*', True)
    harness.main()
