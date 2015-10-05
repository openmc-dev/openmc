#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import TestHarness, PyAPITestHarness
import openmc


class FilterGroupTransferTestHarness(PyAPITestHarness):
    def _build_inputs(self):
        filt1 = openmc.Filter(type='energy',
                              bins=(0.0, 0.253e-6, 1.0e-3, 1.0, 20.0))
        filt2 = openmc.Filter(type='energyout',
                              bins=(0.0, 0.253e-6, 1.0e-3, 1.0, 20.0))
        tally = openmc.Tally(tally_id=1)
        tally.add_filter(filt1)
        tally.add_filter(filt2)
        tally.add_score('scatter')
        tally.add_score('nu-fission')
        self._input_set.tallies = openmc.TalliesFile()
        self._input_set.tallies.add_tally(tally)

        super(FilterGroupTransferTestHarness, self)._build_inputs()

    def _cleanup(self):
        super(FilterGroupTransferTestHarness, self)._cleanup()
        f = os.path.join(os.getcwd(), 'tallies.xml')
        if os.path.exists(f): os.remove(f)


if __name__ == '__main__':
    harness = FilterGroupTransferTestHarness('statepoint.10.*', True)
    harness.main()
