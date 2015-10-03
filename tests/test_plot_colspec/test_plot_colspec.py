#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, os.pardir))
from testing_harness import PlotTestHarness


if __name__ == '__main__':
    harness = PlotTestHarness(('1_plot.ppm', ))
    harness.main()
