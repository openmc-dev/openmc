#!/usr/bin/env python

import os
import sys
import glob
import hashlib
sys.path.insert(0, os.pardir)
from testing_harness import PyAPITestHarness
import openmc
import openmc.mgxs


class IsoInLabTestHarness(PyAPITestHarness):

    def _build_inputs(self):
        """Write input XML files with iso-in-lab scattering."""

        self._input_set.build_default_materials_and_geometry()
        self._input_set.build_default_settings()
        self._input_set.materials.make_isotropic_in_lab()
        self._input_set.export()


if __name__ == '__main__':
    harness = IsoInLabTestHarness('statepoint.10.*')
    harness.main()
