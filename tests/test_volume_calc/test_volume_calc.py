#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import PyAPITestHarness
import openmc


class VolumeTest(PyAPITestHarness):
    def _build_inputs(self):
        # Define materials
        water = openmc.Material(1)
        water.add_nuclide('H1', 2.0)
        water.add_nuclide('O16', 1.0)
        water.add_nuclide('B10', 0.0001)
        water.add_s_alpha_beta('c_H_in_H2O', '71t')
        water.set_density('g/cc', 1.0)

        fuel = openmc.Material(2)
        fuel.add_nuclide('U235', 1.0)
        fuel.add_nuclide('Mo99', 0.1)
        fuel.set_density('g/cc', 4.5)

        materials = openmc.Materials((water, fuel))
        materials.default_xs = '71c'
        materials.export_to_xml()

        cyl = openmc.ZCylinder(1, R=1.0, boundary_type='vacuum')
        top_sphere = openmc.Sphere(2, z0=5., R=1., boundary_type='vacuum')
        top_plane = openmc.ZPlane(3, z0=5.)
        bottom_sphere = openmc.Sphere(4, z0=-5., R=1., boundary_type='vacuum')
        bottom_plane = openmc.ZPlane(5, z0=-5.)

        # Define geometry
        inside_cyl = openmc.Cell(1, fill=fuel, region=-cyl & -top_plane & +bottom_plane)
        top_hemisphere = openmc.Cell(2, fill=water, region=-top_sphere & +top_plane)
        bottom_hemisphere = openmc.Cell(3, fill=water, region=-bottom_sphere & -top_plane)
        root = openmc.Universe(0, cells=(inside_cyl, top_hemisphere, bottom_hemisphere))

        geometry = openmc.Geometry()
        geometry.root_universe = root
        geometry.export_to_xml()

        # Set up stochastic volume calculation
        vol_calc = openmc.VolumeCalculation(
            [inside_cyl, top_hemisphere, bottom_hemisphere],
            100000)

        # Define settings
        settings = openmc.Settings()
        settings.particles = 1000
        settings.batches = 4
        settings.inactive = 0
        settings.source = openmc.Source(space=openmc.stats.Box(
            [-1., -1., -5.], [1., 1., 5.]))
        settings.volume_calculations = vol_calc
        settings.export_to_xml()

    def _get_results(self):
        # Read the statepoint file.
        statepoint = os.path.join(os.getcwd(), self._sp_name)
        sp = openmc.StatePoint(statepoint)

        # Write out k-combined.
        outstr = 'k-combined: {:12.6e} {:12.6e}\n'.format(*sp.k_combined)

        # Read volume calculation results
        vol = openmc.VolumeCalculation.from_hdf5(
            os.path.join(os.getcwd(), 'volume_1.h5'))

        # Write cell volumes and total # of atoms for each nuclide
        for cell_id, results in sorted(vol.results.items()):
            outstr += 'Cell {0}: {1[0]:.4f} +/- {1[1]:.4f} cm^3\n'.format(
                cell_id, results['volume'])
        outstr += str(vol.atoms_dataframe) + '\n'

        return outstr

if __name__ == '__main__':
    harness = VolumeTest('statepoint.4.h5')
    harness.main()
