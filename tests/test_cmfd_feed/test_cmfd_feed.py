#!/usr/bin/env python

import sys
sys.path.insert(0, '..')
from testing_harness import CMFDTestHarness


if __name__ == '__main__':
    harness = CMFDTestHarness('statepoint.20.*', True)
    harness.main()
