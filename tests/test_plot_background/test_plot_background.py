#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import PlotTestHarness


if __name__ == '__main__':
    harness = PlotTestHarness(('1_plot.ppm', ))
    harness.main()
