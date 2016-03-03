#!/usr/bin/env python

import os
import sys
import glob
import hashlib
sys.path.insert(0, os.pardir)
from testing_harness import PyAPITestHarness
import openmc
from openmc.source import Source
from openmc.stats import Box


class AsymmetricLatticeTestHarness(PyAPITestHarness):

    def _build_inputs(self):
        """Build an axis-asymmetric lattice of fuel assemblies"""

        # Build full core geometry from underlying input set
        self._input_set.build_default_materials_and_geometry()

        # Extract universes encapsulating fuel and water assemblies
        geometry = self._input_set.geometry.geometry
        water = geometry.get_universes_by_name('water assembly (hot)')[0]
        fuel = geometry.get_universes_by_name('fuel assembly (hot)')[0]

        # Construct a 3x3 lattice of fuel assemblies
        core_lat = openmc.RectLattice(name='3x3 Core Lattice', lattice_id=202)
        core_lat.dimension = (3, 3)
        core_lat.lower_left = (-32.13, -32.13)
        core_lat.pitch = (21.42, 21.42)
        core_lat.universes = [[fuel,  water, water],
                              [fuel,  fuel,  fuel],
                              [water, water, water]]

        # Create bounding surfaces
        min_x = openmc.XPlane(x0=-32.13, boundary_type='reflective')
        max_x = openmc.XPlane(x0=+32.13, boundary_type='reflective')
        min_y = openmc.YPlane(y0=-32.13, boundary_type='reflective')
        max_y = openmc.YPlane(y0=+32.13, boundary_type='reflective')
        min_z = openmc.ZPlane(z0=0, boundary_type='reflective')
        max_z = openmc.ZPlane(z0=+32.13, boundary_type='reflective')

        # Define root universe
        root_univ = openmc.Universe(universe_id=0, name='root universe')
        root_cell = openmc.Cell(cell_id=1)
        root_cell.region = +min_x & -max_x & +min_y & -max_y & +min_z & -max_z
        root_cell.fill = core_lat
        root_univ.add_cell(root_cell)

        # Over-ride geometry in the input set with this 3x3 lattice
        self._input_set.geometry.geometry.root_universe = root_univ

        # Initialize a "distribcell" filter for the fuel pin cell
        distrib_filter = openmc.Filter(type='distribcell', bins=[27])

        # Initialize the tallies
        tally = openmc.Tally(name='distribcell tally', tally_id=27)
        tally.add_filter(distrib_filter)
        tally.add_score('nu-fission')

        # Initialize the tallies file
        tallies_file = openmc.TalliesFile()
        tallies_file.add_tally(tally)

        # Assign the tallies file to the input set
        self._input_set.tallies = tallies_file

        # Build default settings
        self._input_set.build_default_settings()

        # Specify summary output and correct source sampling box
        source = Source(space=Box([-32, -32, 0], [32, 32, 32]))
        source.space.only_fissionable = True
        self._input_set.settings.source = source
        self._input_set.settings.output = {'summary': True}

        # Write input XML files
        self._input_set.export()

    def _get_results(self, hash_output=True):
        """Digest info in statepoint and summary and return as a string."""

        # Read the statepoint file
        statepoint = glob.glob(os.path.join(os.getcwd(), self._sp_name))[0]
        sp = openmc.StatePoint(statepoint)

        # Read the summary file
        summary = glob.glob(os.path.join(os.getcwd(), 'summary.h5'))[0]
        su = openmc.Summary(summary)
        sp.link_with_summary(su)

        # Extract the tally of interest
        tally = sp.get_tally(name='distribcell tally')

        # Create a string of all mean, std. dev. values for both tallies
        outstr = ''
        outstr += ', '.join(map(str, tally.mean.flatten())) + '\n'
        outstr += ', '.join(map(str, tally.std_dev.flatten())) + '\n'

        # Extract fuel assembly lattices from the summary
        core = su.get_cell_by_id(1)
        fuel = su.get_cell_by_id(80)
        fuel = fuel.fill
        core = core.fill

        # Append a string of lattice distribcell offsets to the string
        outstr += ', '.join(map(str, fuel.offsets.flatten())) + '\n'
        outstr += ', '.join(map(str, core.offsets.flatten())) + '\n'

        # Hash the results if necessary
        if hash_output:
            sha512 = hashlib.sha512()
            sha512.update(outstr.encode('utf-8'))
            outstr = sha512.hexdigest()

        return outstr

    def _cleanup(self):
        super(AsymmetricLatticeTestHarness, self)._cleanup()
        f = os.path.join(os.getcwd(), 'tallies.xml')
        if os.path.exists(f): os.remove(f)


if __name__ == '__main__':
    harness = AsymmetricLatticeTestHarness('statepoint.10.h5', True)
    harness.main()
