#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import TestHarness, PyAPITestHarness
import openmc


class MGBasicTestHarness(PyAPITestHarness):
    def _build_inputs(self):
        """Write input XML files."""
        self._input_set.build_default_materials_and_geometry()
        self._input_set.build_default_settings()
        self._input_set.settings.survival_biasing = True
        self._input_set.export()


if __name__ == '__main__':
    harness = MGBasicTestHarness('statepoint.10.h5', False, mg=True)
    harness.main()
