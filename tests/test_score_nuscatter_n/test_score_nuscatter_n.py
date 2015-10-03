#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import TestHarness, PyAPITestHarness
import openmc


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

        self._input_set.build_default_materials_and_geometry()
        self._input_set.build_default_settings()

        self._input_set.settings.inactive = 0

        self._input_set.export()

    def _cleanup(self):
        super(ScoreNuScatterNTestHarness, self)._cleanup()
        f = os.path.join(os.getcwd(), 'tallies.xml')
        if os.path.exists(f): os.remove(f)


if __name__ == '__main__':
    harness = ScoreNuScatterNTestHarness('statepoint.10.*', True)
    harness.main()
