import os
import glob
import sys

import openmc

from tests.testing_harness import PyAPITestHarness


class VolumeTest(PyAPITestHarness):
    def _build_inputs(self):
        # Define materials
        water = openmc.Material(1)
        water.add_nuclide('H1', 2.0)
        water.add_nuclide('O16', 1.0)
        water.add_nuclide('B10', 0.0001)
        water.add_s_alpha_beta('c_H_in_H2O')
        water.set_density('g/cc', 1.0)

        fuel = openmc.Material(2)
        fuel.add_nuclide('U235', 1.0)
        fuel.add_nuclide('Mo99', 0.1)
        fuel.set_density('g/cc', 4.5)

        materials = openmc.Materials((water, fuel))
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

        geometry = openmc.Geometry(root)
        geometry.export_to_xml()

        # Set up stochastic volume calculation
        ll, ur = root.bounding_box
        vol_calcs = [
            openmc.VolumeCalculation(list(root.cells.values()), 100000),
            openmc.VolumeCalculation([water, fuel], 100000, ll, ur),
            openmc.VolumeCalculation([root], 100000, ll, ur)
        ]

        # Define settings
        settings = openmc.Settings()
        settings.run_mode = 'volume'
        settings.volume_calculations = vol_calcs
        settings.export_to_xml()

    def _get_results(self):
        outstr = ''
        for i, filename in enumerate(sorted(glob.glob('volume_*.h5'))):
            outstr += 'Volume calculation {}\n'.format(i)

            # Read volume calculation results
            volume_calc = openmc.VolumeCalculation.from_hdf5(filename)

            # Write cell volumes and total # of atoms for each nuclide
            for uid, volume in sorted(volume_calc.volumes.items()):
                outstr += 'Domain {}: {} cm^3\n'.format(uid, volume)
            outstr += str(volume_calc.atoms_dataframe) + '\n'

        return outstr

    def _test_output_created(self):
        pass

def test_volume_calc():
    harness = VolumeTest('')
    harness.main()
