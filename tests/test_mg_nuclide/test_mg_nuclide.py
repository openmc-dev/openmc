#!/usr/bin/env python

import os
import sys

sys.path.insert(0, os.pardir)
from testing_harness import PyAPITestHarness
from input_set import MGInputSet


class MGNuclideTestHarness(PyAPITestHarness):
    def __init__(self, statepoint_name, tallies_present, mg=False):
        PyAPITestHarness.__init__(self, statepoint_name, tallies_present)
        self._input_set = MGInputSet()

    def _build_inputs(self):
        """Write input XML files."""
        self._input_set.build_default_materials_and_geometry(as_macro=False)
        self._input_set.build_default_settings()
        self._input_set.export()


if __name__ == '__main__':
    harness = MGNuclideTestHarness('statepoint.10.*', False, mg=True)
    harness.main()
