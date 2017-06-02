#!/usr/bin/env python

import os
import sys
import glob
sys.path.insert(0, os.pardir)
from testing_harness import PyAPITestHarness
import openmc


class FilterEnergyFunHarness(PyAPITestHarness):
    def __init__(self, *args, **kwargs):
        super(FilterEnergyFunHarness, self).__init__(*args, **kwargs)

        # Add Am241 to the fuel.
        self._model.materials[1].add_nuclide('Am241', 1e-7)

        # Define Am242m / Am242 branching ratio from ENDF/B-VII.1 data.
        x = [1e-5, 3.69e-1, 1e3, 1e5, 6e5, 1e6, 2e6, 4e6, 3e7]
        y = [0.1, 0.1, 0.1333, 0.158, 0.18467, 0.25618, 0.4297, 0.48, 0.48]

        # Make an EnergyFunctionFilter directly from the x and y lists.
        filt1 = openmc.EnergyFunctionFilter(x, y)

        # Also make a filter with the .from_tabulated1d constructor.  Make sure
        # the filters are identical.
        tab1d = openmc.data.Tabulated1D(x, y)
        filt2 = openmc.EnergyFunctionFilter.from_tabulated1d(tab1d)
        assert filt1 == filt2, 'Error with the .from_tabulated1d constructor'

        # Make tallies.
        tallies = [openmc.Tally(1), openmc.Tally(2)]
        for t in tallies:
            t.scores = ['(n,gamma)']
            t.nuclides = ['Am241']
        tallies[1].filters = [filt1]
        self._model.tallies = tallies

    def _get_results(self):
        # Read the statepoint file.
        statepoint = glob.glob(os.path.join(os.getcwd(), self._sp_name))[0]
        sp = openmc.StatePoint(statepoint)

        # Use tally arithmetic to compute the branching ratio.
        br_tally = sp.tallies[2] / sp.tallies[1]

        # Output the tally in a Pandas DataFrame.
        return br_tally.get_pandas_dataframe().to_string() + '\n'


if __name__ == '__main__':
    harness = FilterEnergyFunHarness('statepoint.10.h5')
    harness.main()
