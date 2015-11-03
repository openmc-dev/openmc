#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
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

        self._input_set.build_default_materials_and_geometry()
        self._input_set.build_default_settings()

        self._input_set.settings.inactive = 0

        self._input_set.export()

    def _cleanup(self):
        super(ScoreNuScatterTestHarness, self)._cleanup()
        f = os.path.join(os.getcwd(), 'tallies.xml')
        if os.path.exists(f): os.remove(f)


if __name__ == '__main__':
    harness = ScoreNuScatterTestHarness('statepoint.10.*', True)
    harness.main()
