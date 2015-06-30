#!/usr/bin/env python

import sys
sys.path.insert(0, '..')
from testing_harness import HashedTestHarness


if __name__ == '__main__':
    harness = HashedTestHarness('statepoint.10.*', True)
    harness.main()
