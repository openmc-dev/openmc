#!/usr/bin/env python

import sys
sys.path.insert(0, '..')
from testing_harness import TestHarness


if __name__ == '__main__':
    harness = TestHarness('statepoint.5.*', True)
    harness.main()
