#!/usr/bin/env python

import os
import sys
import glob
import hashlib
import itertools
sys.path.insert(0, os.pardir)
from testing_harness import PyAPITestHarness
import openmc


class TallySliceMergeTestHarness(PyAPITestHarness):
    def _build_inputs(self):

        # The summary.h5 file needs to be created to read in the tallies
        self._input_set.settings.output = {'summary': True}

        # Initialize the tallies file
        tallies_file = openmc.TalliesFile()

        # Define nuclides and scores to add to both tallies
        self.nuclides = ['U-235', 'U-238']
        self.scores = ['fission', 'nu-fission']

        # Define filters for energy and spatial domain

        low_energy = openmc.Filter(type='energy', bins=[0., 0.625e-6])
        high_energy = openmc.Filter(type='energy', bins=[0.625e-6, 20.])
        merged_energies = low_energy.merge(high_energy)

        cell_21 = openmc.Filter(type='cell', bins=[21])
        cell_27 = openmc.Filter(type='cell', bins=[27])
        distribcell_filter = openmc.Filter(type='distribcell', bins=[21])

        self.cell_filters = [cell_21, cell_27]
        self.energy_filters = [low_energy, high_energy]

        # Initialize cell tallies with filters, nuclides and scores
        tallies = []
        for cell_filter in self.energy_filters:
            for energy_filter in self.cell_filters:
                for nuclide in self.nuclides:
                    for score in self.scores:
                        tally = openmc.Tally()
                        tally.estimator = 'tracklength'
                        tally.add_score(score)
                        tally.add_nuclide(nuclide)
                        tally.add_filter(cell_filter)
                        tally.add_filter(energy_filter)
                        tallies.append(tally)

        # Merge all cell tallies together
        while len(tallies) != 1:
            halfway = int(len(tallies) / 2)
            zip_split = zip(tallies[:halfway], tallies[halfway:])
            tallies = list(map(lambda xy: xy[0].merge(xy[1]), zip_split))

        # Specify a name for the tally
        tallies[0].name = 'cell tally'

        # Initialize a distribcell tally
        distribcell_tally = openmc.Tally(name='distribcell tally')
        distribcell_tally.estimator = 'tracklength'
        distribcell_tally.add_filter(distribcell_filter)
        distribcell_tally.add_filter(merged_energies)
        for score in self.scores:
            distribcell_tally.add_score(score)
        for nuclide in self.nuclides:
            distribcell_tally.add_nuclide(nuclide)

        # Add tallies to a TalliesFile
        tallies_file = openmc.TalliesFile()
        tallies_file.add_tally(tallies[0])
        tallies_file.add_tally(distribcell_tally)

        # Export tallies to file
        self._input_set.tallies = tallies_file
        super(TallySliceMergeTestHarness, self)._build_inputs()

    def _get_results(self, hash_output=False):
        """Digest info in the statepoint and return as a string."""

        # Read the statepoint file.
        statepoint = glob.glob(os.path.join(os.getcwd(), self._sp_name))[0]
        sp = openmc.StatePoint(statepoint)

        # Read the summary file.
        summary = glob.glob(os.path.join(os.getcwd(), 'summary.h5'))[0]
        su = openmc.Summary(summary)
        sp.link_with_summary(su)

        # Extract the cell tally
        tallies = [sp.get_tally(name='cell tally')]

        # Slice the tallies by cell filter bins
        cell_filter_prod = itertools.product(tallies, self.cell_filters)
        tallies = map(lambda tf: tf[0].get_slice(filters=[tf[1].type], 
                      filter_bins=[tf[1].get_bin(0)]), cell_filter_prod)

        # Slice the tallies by energy filter bins
        energy_filter_prod = itertools.product(tallies, self.energy_filters)
        tallies = map(lambda tf: tf[0].get_slice(filters=[tf[1].type],
                      filter_bins=[(tf[1].get_bin(0),)]), energy_filter_prod)

        # Slice the tallies by nuclide
        nuclide_prod = itertools.product(tallies, self.nuclides)
        tallies = map(lambda tn: tn[0].get_slice(nuclides=[tn[1]]), nuclide_prod)

        # Slice the tallies by score
        score_prod = itertools.product(tallies, self.scores)
        tallies = map(lambda ts: ts[0].get_slice(scores=[ts[1]]), score_prod)
        tallies = list(tallies)

        # Initialize an output string
        outstr = ''

        # Append sliced Tally Pandas DataFrames to output string
        for tally in tallies:
            df = tally.get_pandas_dataframe()
            outstr += df.to_string()

        # Merge all tallies together
        while len(tallies) != 1:
            halfway = int(len(tallies) / 2)
            zip_split = zip(tallies[:halfway], tallies[halfway:])
            tallies = list(map(lambda xy: xy[0].merge(xy[1]), zip_split))

        # Append merged Tally Pandas DataFrame to output string
        df = tallies[0].get_pandas_dataframe()
        outstr += df.to_string()

        # Extract the distribcell tally
        distribcell_tally = sp.get_tally(name='distribcell tally')

        # Sum up a few subdomains from the distribcell tally 
        sum1 = distribcell_tally.summation(filter_type='distribcell', 
                                           filter_bins=[0,100,2000,30000])
        # Sum up a few subdomains from the distribcell tally
        sum2 = distribcell_tally.summation(filter_type='distribcell', 
                                           filter_bins=[500,5000,50000])

        # Merge the distribcell tally slices
        merge_tally = sum1.merge(sum2)

        # Append merged Tally Pandas DataFrame to output string
        df = merge_tally.get_pandas_dataframe()
        outstr += df.to_string()

        # Hash the results if necessary
        if hash_output:
            sha512 = hashlib.sha512()
            sha512.update(outstr.encode('utf-8'))
            outstr = sha512.hexdigest()

        return outstr

    def _cleanup(self):
        super(TallySliceMergeTestHarness, self)._cleanup()
        f = os.path.join(os.getcwd(), 'tallies.xml')
        if os.path.exists(f): os.remove(f)

if __name__ == '__main__':
    harness = TallySliceMergeTestHarness('statepoint.10.h5', True)
    harness.main()
