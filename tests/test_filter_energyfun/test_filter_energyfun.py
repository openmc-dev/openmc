#!/usr/bin/env python

import os
import sys
import glob
import hashlib
sys.path.insert(0, os.pardir)
from testing_harness import PyAPITestHarness
import openmc


class FilterEnergyFunHarness(PyAPITestHarness):
    def _build_inputs(self):
        # Build the default material, geometry, settings inputs.
        self._input_set.build_default_materials_and_geometry()
        self._input_set.build_default_settings()

        # Add Am241 to the fuel.
        self._input_set.materials[1].add_nuclide('Am241', 1e-7)

        # Make an EnergyFunctionFilter for the Am242m / Am242 branching ratio.
        x = [1e-5, 3.69e-1, 1e3, 1e5, 6e5, 1e6, 2e6, 4e6, 3e7]
        y = [0.1, 0.1, 0.1333, 0.158, 0.18467, 0.25618, 0.4297, 0.48, 0.48]
        filt = openmc.EnergyFunctionFilter(x, y)

        # Make tallies.
        tallies = [openmc.Tally() for i in range(2)]
        for t in tallies:
            t.scores = ['102']
            t.nuclides = ['Am241']
        tallies[1].filters = [filt]
        self._input_set.tallies = openmc.Tallies(tallies)

        # Export inputs to xml.
        self._input_set.export()

    def _get_results(self):
        # Read the statepoint file.
        statepoint = glob.glob(os.path.join(os.getcwd(), self._sp_name))[0]
        sp = openmc.StatePoint(statepoint)

        # Use tally arithmetic to compute the branching ratio.
        br_tally = sp.tallies[10001] / sp.tallies[10000]

        # Output the tally in a Pandas DataFrame.
        return br_tally.get_pandas_dataframe().to_string() + '\n'

    def _cleanup(self):
        super(FilterEnergyFunHarness, self)._cleanup()
        f = os.path.join(os.getcwd(), 'tallies.xml')
        if os.path.exists(f): os.remove(f)


if __name__ == '__main__':
    harness = FilterEnergyFunHarness('statepoint.10.*', True)
    harness.main()
