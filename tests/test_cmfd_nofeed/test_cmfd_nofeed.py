#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, os.pardir))
from testing_harness import CMFDTestHarness


if __name__ == '__main__':
    harness = CMFDTestHarness('statepoint.20.*', True)
    harness.main()
