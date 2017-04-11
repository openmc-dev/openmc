#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import PyAPITestHarness
from openmc.examples import slab_mg


if __name__ == '__main__':
    model = slab_mg()
    harness = PyAPITestHarness('statepoint.10.h5', model)
    harness.main()
