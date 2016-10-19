#!/usr/bin/env python

import os
import sys

sys.path.insert(0, os.pardir)
from testing_harness import PyAPITestHarness
from input_set import MGInputSet


class MGMaxOrderTestHarness(PyAPITestHarness):
    def __init__(self, statepoint_name, tallies_present, mg=False):
        PyAPITestHarness.__init__(self, statepoint_name, tallies_present)
        self._input_set = MGInputSet()

    def _build_inputs(self):
        """Write input XML files."""
        reps = ['iso']
        self._input_set.build_default_materials_and_geometry(reps=reps)
        self._input_set.build_default_settings()
        # Set P1 scattering
        self._input_set.settings.max_order = 1
        self._input_set.export()


if __name__ == '__main__':
    harness = MGMaxOrderTestHarness('statepoint.10.*', False, mg=True)
    harness.main()
