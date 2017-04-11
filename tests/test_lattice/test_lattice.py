#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import TestHarness


if __name__ == '__main__':
    harness = TestHarness('statepoint.10.h5')
    harness.main()
