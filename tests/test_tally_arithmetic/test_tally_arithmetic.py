#!/usr/bin/env python

import os
import sys
import glob
import hashlib
sys.path.insert(0, os.pardir)
from testing_harness import TestHarness
import openmc


class TallyArithmeticTestHarness(TestHarness):
    def _build_inputs(self):

        u235 = openmc.Nuclide('U-235')
        u238 = openmc.Nuclide('U-238')
        pu239 = openmc.Nuclide('Pu-239')

        # Instantiate energy filter
        energy_filter = openmc.Filter(type='energy', bins=[0., 20.])

        # Create tallies
        tally_1 = openmc.Tally(name='tally 1')
        tally_1.add_filter(openmc.Filter(type='cell', bins=[1]))
        tally_1.add_filter(energy_filter)
        tally_1.add_score('fission')
        tally_1.add_score('nu-fission')
        tally_1.add_nuclide(u235)
        tally_1.add_nuclide(u238)

        tally_2 = openmc.Tally(name='tally 2')
        tally_2.add_filter(openmc.Filter(type='material', bins=[1]))
        tally_2.add_filter(energy_filter)
        tally_2.add_score('fission')
        tally_2.add_score('absorption')
        tally_2.add_nuclide(u235)
        tally_2.add_nuclide(pu239)

        # Export tallies to file
        tallies_file = openmc.TalliesFile()
        tallies_file.add_tally(tally_1)
        tallies_file.add_tally(tally_2)
        tallies_file.export_to_xml()

    def _get_results(self, hash_output=False):
        """Digest info in the statepoint and return as a string."""

        # Read the statepoint file.
        statepoint = glob.glob(os.path.join(os.getcwd(), self._sp_name))[0]
        sp = openmc.StatePoint(statepoint)

        # Read the summary file.
        summary = glob.glob(os.path.join(os.getcwd(), 'summary.h5'))[0]
        su = openmc.Summary(summary)
        sp.link_with_summary(su)

        # Load the tallies
        tally_1 = sp.get_tally(name='tally 1')
        tally_2 = sp.get_tally(name='tally 2')

        # Perform all the tally arithmetic operations and output results
        outstr = ''
        tally_3 = tally_1 + tally_2
        outstr += tally_3.__repr__()
        outstr += str(tally_3.mean)

        tally_3 = tally_1 - tally_2
        outstr += tally_3.__repr__()
        outstr += str(tally_3.mean)

        tally_3 = tally_1 * tally_2
        outstr += tally_3.__repr__()
        outstr += str(tally_3.mean)

        tally_3 = tally_1 / tally_2
        outstr += tally_3.__repr__()
        outstr += str(tally_3.mean)

        print(outstr)

        # Hash the results if necessary
        if hash_output:
            sha512 = hashlib.sha512()
            sha512.update(outstr.encode('utf-8'))
            outstr = sha512.hexdigest()

        return outstr

if __name__ == '__main__':
    harness = TallyArithmeticTestHarness('statepoint.10.*', True)
    harness.main()
