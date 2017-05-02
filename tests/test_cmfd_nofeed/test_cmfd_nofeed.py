#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import CMFDTestHarness


if __name__ == '__main__':
    harness = CMFDTestHarness('statepoint.20.h5')
    harness.main()
