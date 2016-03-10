#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import TestHarness, PyAPITestHarness
import openmc


class MGBasicTestHarness(PyAPITestHarness):
    def _build_inputs(self):
        super(MGBasicTestHarness, self)._build_inputs()


if __name__ == '__main__':
    harness = MGBasicTestHarness('statepoint.10.*', False, mg=True)
    harness.main()
