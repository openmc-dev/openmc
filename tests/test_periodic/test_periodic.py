#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import PyAPITestHarness
import openmc


class PeriodicTest(PyAPITestHarness):
    def _build_inputs(self):
        # Define materials
        water = openmc.Material(1)
        water.add_nuclide('H-1', 2.0)
        water.add_nuclide('O-16', 1.0)
        water.add_s_alpha_beta('HH2O', '71t')
        water.set_density('g/cc', 1.0)

        fuel = openmc.Material(2)
        fuel.add_nuclide('U-235', 1.0)
        fuel.set_density('g/cc', 4.5)

        materials = openmc.Materials((water, fuel))
        materials.default_xs = '71c'
        materials.export_to_xml()

        # Define geometry
        x_min = openmc.XPlane(1, x0=-5., boundary_type='periodic')
        x_max = openmc.XPlane(2, x0=5., boundary_type='periodic')
        x_max.periodic_surface = x_min

        y_min = openmc.YPlane(3, y0=-5., boundary_type='periodic')
        y_max = openmc.YPlane(4, y0=5., boundary_type='periodic')

        z_min = openmc.ZPlane(5, z0=-5., boundary_type='reflective')
        z_max = openmc.ZPlane(6, z0=5., boundary_type='reflective')
        z_cyl = openmc.ZCylinder(7, x0=-2.5, y0=2.5, R=2.0)

        outside_cyl = openmc.Cell(1, fill=water, region=(
            +x_min & -x_max & +y_min & -y_max & +z_min & -z_max & +z_cyl))
        inside_cyl = openmc.Cell(2, fill=fuel, region=+z_min & -z_max & -z_cyl)
        root_universe = openmc.Universe(0, cells=(outside_cyl, inside_cyl))

        geometry = openmc.Geometry()
        geometry.root_universe = root_universe
        geometry.export_to_xml()

        # Define settings
        settings = openmc.Settings()
        settings.particles = 1000
        settings.batches = 4
        settings.inactive = 0
        settings.source = openmc.Source(space=openmc.stats.Box(
            *outside_cyl.region.bounding_box))
        settings.export_to_xml()


if __name__ == '__main__':
    harness = PeriodicTest('statepoint.4.h5')
    harness.main()
