#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import TestHarness, PyAPITestHarness
import openmc


class FilterCellbornTestHarness(PyAPITestHarness):
    def _build_inputs(self):
        filt = openmc.Filter(type='cellborn', bins=(10, 21, 22, 23))
        tally = openmc.Tally(tally_id=1)
        tally.add_filter(filt)
        tally.add_score('total')
        self._input_set.tallies = openmc.TalliesFile()
        self._input_set.tallies.add_tally(tally)

        super(FilterCellbornTestHarness, self)._build_inputs()

    def _cleanup(self):
        super(FilterCellbornTestHarness, self)._cleanup()
        f = os.path.join(os.getcwd(), 'tallies.xml')
        if os.path.exists(f): os.remove(f)


if __name__ == '__main__':
    harness = FilterCellbornTestHarness('statepoint.10.*', True)
    harness.main()
