import os
import glob
import hashlib

import openmc

from tests.testing_harness import PyAPITestHarness


class AsymmetricLatticeTestHarness(PyAPITestHarness):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Extract universes encapsulating fuel and water assemblies
        geometry = self._model.geometry
        water = geometry.get_universes_by_name('water assembly (hot)')[0]
        fuel = geometry.get_universes_by_name('fuel assembly (hot)')[0]

        # Construct a 3x3 lattice of fuel assemblies
        core_lat = openmc.RectLattice(name='3x3 Core Lattice', lattice_id=202)
        core_lat.lower_left = (-32.13, -32.13)
        core_lat.pitch = (21.42, 21.42)
        core_lat.universes = [[fuel,  water, water],
                              [fuel,  fuel,  fuel],
                              [water, water, water]]

        # Create bounding surfaces
        min_x = openmc.XPlane(-32.13, boundary_type='reflective')
        max_x = openmc.XPlane(+32.13, boundary_type='reflective')
        min_y = openmc.YPlane(-32.13, boundary_type='reflective')
        max_y = openmc.YPlane(+32.13, boundary_type='reflective')
        min_z = openmc.ZPlane(0, boundary_type='reflective')
        max_z = openmc.ZPlane(+32.13, boundary_type='reflective')

        # Define root universe
        root_univ = openmc.Universe(universe_id=0, name='root universe')
        root_cell = openmc.Cell(cell_id=1)
        root_cell.region = +min_x & -max_x & +min_y & -max_y & +min_z & -max_z
        root_cell.fill = core_lat
        root_univ.add_cell(root_cell)

        # Over-ride geometry in the input set with this 3x3 lattice
        self._model.geometry.root_universe = root_univ

        # Initialize a "distribcell" filter for the fuel pin cell
        distrib_filter = openmc.DistribcellFilter(27)

        # Initialize the tallies
        tally = openmc.Tally(name='distribcell tally', tally_id=27)
        tally.filters.append(distrib_filter)
        tally.scores.append('nu-fission')

        # Assign the tallies file to the input set
        self._model.tallies.append(tally)

        # Specify summary output and correct source sampling box
        self._model.settings.source = openmc.Source(space=openmc.stats.Box(
            [-32, -32, 0], [32, 32, 32], only_fissionable = True))

    def _get_results(self, hash_output=True):
        """Digest info in statepoint and summary and return as a string."""

        # Read the statepoint file
        statepoint = glob.glob(os.path.join(os.getcwd(), self._sp_name))[0]
        sp = openmc.StatePoint(statepoint)

        # Extract the tally of interest
        tally = sp.get_tally(name='distribcell tally')

        # Create a string of all mean, std. dev. values for both tallies
        outstr = ''
        outstr += '\n'.join(map('{:.8e}'.format, tally.mean.flatten())) + '\n'
        outstr += '\n'.join(map('{:.8e}'.format, tally.std_dev.flatten())) + '\n'

        # Extract fuel assembly lattices from the summary
        cells = sp.summary.geometry.get_all_cells()
        fuel_cell = cells[27]

        # Append a string of lattice distribcell offsets to the string
        outstr += '\n'.join(fuel_cell.paths) + '\n'

        # Hash the results if necessary
        if hash_output:
            sha512 = hashlib.sha512()
            sha512.update(outstr.encode('utf-8'))
            outstr = sha512.hexdigest()

        return outstr


def test_asymmetric_lattice():
    harness = AsymmetricLatticeTestHarness('statepoint.10.h5')
    harness.main()
