import hashlib
import itertools

import openmc

from tests.testing_harness import PyAPITestHarness


class TallySliceMergeTestHarness(PyAPITestHarness):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Define nuclides and scores to add to both tallies
        self.nuclides = ['U235', 'U238']
        self.scores = ['fission', 'nu-fission']

        # Define filters for energy and spatial domain

        low_energy = openmc.EnergyFilter([0., 0.625])
        high_energy = openmc.EnergyFilter([0.625, 20.e6])
        merged_energies = low_energy.merge(high_energy)

        cell_21 = openmc.CellFilter(21)
        cell_27 = openmc.CellFilter(27)
        distribcell_filter = openmc.DistribcellFilter(21)

        mesh = openmc.RegularMesh(name='mesh')
        mesh.dimension = [2, 2]
        mesh.lower_left = [-50., -50.]
        mesh.upper_right = [+50., +50.]
        mesh_filter = openmc.MeshFilter(mesh)

        self.cell_filters = [cell_21, cell_27]
        self.energy_filters = [low_energy, high_energy]

        # Initialize cell tallies with filters, nuclides and scores
        tallies = []
        for energy_filter in self.energy_filters:
            for cell_filter in self.cell_filters:
                for nuclide in self.nuclides:
                    for score in self.scores:
                        tally = openmc.Tally()
                        tally.estimator = 'tracklength'
                        tally.scores.append(score)
                        tally.nuclides.append(nuclide)
                        tally.filters.append(cell_filter)
                        tally.filters.append(energy_filter)
                        tallies.append(tally)

        # Merge all cell tallies together
        while len(tallies) != 1:
            halfway = len(tallies) // 2
            zip_split = zip(tallies[:halfway], tallies[halfway:])
            tallies = list(map(lambda xy: xy[0].merge(xy[1]), zip_split))

        # Specify a name for the tally
        tallies[0].name = 'cell tally'

        # Initialize a distribcell tally
        distribcell_tally = openmc.Tally(name='distribcell tally')
        distribcell_tally.estimator = 'tracklength'
        distribcell_tally.filters = [distribcell_filter, merged_energies]
        for score in self.scores:
            distribcell_tally.scores.append(score)
        for nuclide in self.nuclides:
            distribcell_tally.nuclides.append(nuclide)

        mesh_tally = openmc.Tally(name='mesh tally')
        mesh_tally.estimator = 'tracklength'
        mesh_tally.filters = [mesh_filter, merged_energies]
        mesh_tally.scores = self.scores
        mesh_tally.nuclides = self.nuclides

        # Add tallies to a Tallies object
        self._model.tallies = [tallies[0], distribcell_tally, mesh_tally]

    def _get_results(self, hash_output=False):
        """Digest info in the statepoint and return as a string."""

        # Read the statepoint file.
        sp = openmc.StatePoint(self._sp_name)

        # Extract the cell tally
        tallies = [sp.get_tally(name='cell tally')]

        # Slice the tallies by cell filter bins
        cell_filter_prod = itertools.product(tallies, self.cell_filters)
        tallies = map(lambda tf: tf[0].get_slice(filters=[type(tf[1])],
                                                 filter_bins=[(tf[1].bins[0],)]),
                      cell_filter_prod)

        # Slice the tallies by energy filter bins
        energy_filter_prod = itertools.product(tallies, self.energy_filters)
        tallies = map(lambda tf: tf[0].get_slice(filters=[type(tf[1])],
                                                 filter_bins=[(tf[1].bins[0],)]),
                      energy_filter_prod)

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
        outstr += df.to_string() + '\n'

        # Extract the distribcell tally
        distribcell_tally = sp.get_tally(name='distribcell tally')

        # Sum up a few subdomains from the distribcell tally
        sum1 = distribcell_tally.summation(filter_type=openmc.DistribcellFilter,
                                           filter_bins=[0, 100, 2000, 30000])
        # Sum up a few subdomains from the distribcell tally
        sum2 = distribcell_tally.summation(filter_type=openmc.DistribcellFilter,
                                           filter_bins=[500, 5000, 50000])

        # Merge the distribcell tally slices
        merge_tally = sum1.merge(sum2)

        # Append merged Tally Pandas DataFrame to output string
        df = merge_tally.get_pandas_dataframe()
        outstr += df.to_string() + '\n'

        # Extract the mesh tally
        mesh_tally = sp.get_tally(name='mesh tally')

        # Sum up a few subdomains from the mesh tally
        sum1 = mesh_tally.summation(filter_type=openmc.MeshFilter,
                                    filter_bins=[(1, 1), (1, 2)])
        # Sum up a few subdomains from the mesh tally
        sum2 = mesh_tally.summation(filter_type=openmc.MeshFilter,
                                    filter_bins=[(2, 1), (2, 2)])

        # Merge the mesh tally slices
        merge_tally = sum1.merge(sum2)

        # Append merged Tally Pandas DataFrame to output string
        df = merge_tally.get_pandas_dataframe()
        outstr += df.to_string() + '\n'

        # Hash the results if necessary
        if hash_output:
            sha512 = hashlib.sha512()
            sha512.update(outstr.encode('utf-8'))
            outstr = sha512.hexdigest()

        return outstr


def test_tally_slice_merge():
    harness = TallySliceMergeTestHarness('statepoint.10.h5')
    harness.main()
