#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import TestHarness, PyAPITestHarness
import openmc


class FilterEnergyoutTestHarness(PyAPITestHarness):
    def _build_inputs(self):
        filt = openmc.Filter(type='energyout',
                             bins=(0.0, 0.253e-6, 1.0e-3, 1.0, 20.0))
        tally = openmc.Tally(tally_id=1)
        tally.add_filter(filt)
        tally.add_score('scatter')
        self._input_set.tallies = openmc.TalliesFile()
        self._input_set.tallies.add_tally(tally)

        super(FilterEnergyoutTestHarness, self)._build_inputs()

    def _cleanup(self):
        super(FilterEnergyoutTestHarness, self)._cleanup()
        f = os.path.join(os.getcwd(), 'tallies.xml')
        if os.path.exists(f): os.remove(f)


if __name__ == '__main__':
    harness = FilterEnergyoutTestHarness('statepoint.10.*', True)
    harness.main()
