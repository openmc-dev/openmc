#!/usr/bin/env python

import sys
sys.path.insert(0, '..')
from testing_harness import PlotTestHarness


if __name__ == '__main__':
    harness = PlotTestHarness('statepoint.10.*', True)
    harness.execute_test()
