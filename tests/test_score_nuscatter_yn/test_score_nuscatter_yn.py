#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import TestHarness, PyAPITestHarness
import openmc


class ScoreNuScatterYNTestHarness(PyAPITestHarness):
    def _build_inputs(self):
        filt = openmc.Filter(type='cell', bins=(21, ))
        t1 = openmc.Tally(tally_id=1)
        t1.add_filter(filt)
        t1.add_score('nu-scatter-0')
        t2 = openmc.Tally(tally_id=2)
        t2.add_filter(filt)
        t2.add_score('nu-scatter-y3')
        self._input_set.tallies = openmc.TalliesFile()
        self._input_set.tallies.add_tally(t1)
        self._input_set.tallies.add_tally(t2)

        self._input_set.build_default_materials_and_geometry()
        self._input_set.build_default_settings()

        self._input_set.settings.inactive = 0

        self._input_set.export()

    def _cleanup(self):
        super(ScoreNuScatterYNTestHarness, self)._cleanup()
        f = os.path.join(os.getcwd(), 'tallies.xml')
        if os.path.exists(f): os.remove(f)


if __name__ == '__main__':
    harness = ScoreNuScatterYNTestHarness('statepoint.10.*', True)
    harness.main()
