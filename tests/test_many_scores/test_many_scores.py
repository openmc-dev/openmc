#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, os.pardir))
from testing_harness import TestHarness


if __name__ == '__main__':
    harness = TestHarness('statepoint.5.*', True)
    harness.main()
