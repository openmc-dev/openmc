#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import PyAPITestHarness


if __name__ == '__main__':
    # Force iso-in-lab scattering.
    harness = PyAPITestHarness('statepoint.10.h5')
    harness._model.materials.make_isotropic_in_lab()
    harness.main()
