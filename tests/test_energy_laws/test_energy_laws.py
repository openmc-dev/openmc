#!/usr/bin/env python

"""The purpose of this test is to provide coverage of Watt, N-body phase space,
and evaporation energy distributions. The only nuclide that uses a Watt fission
spectrum in ENDF/B-VII.1 is U-233. The only nuclide that has a reaction using
the N-body phase space distribution is H-2(n.2n). Several nuclides have
reactions with evaporation spectra. In this test, the material is composed of
U-233, H-2, and Na-23 (which has several reactions with evaporation spectra.

"""

import glob
import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import TestHarness


if __name__ == '__main__':
    harness = TestHarness('statepoint.10.*')
    harness.main()
