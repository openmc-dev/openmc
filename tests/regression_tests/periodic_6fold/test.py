import openmc
import numpy as np

from tests.testing_harness import PyAPITestHarness


class Periodic6FoldTest(PyAPITestHarness):
    def _build_inputs(self):
        # Define materials
        water = openmc.Material(1)
        water.add_nuclide('H1', 2.0)
        water.add_nuclide('O16', 1.0)
        water.add_s_alpha_beta('c_H_in_H2O')
        water.set_density('g/cc', 1.0)

        fuel = openmc.Material(2)
        fuel.add_nuclide('U235', 1.0)
        fuel.set_density('g/cc', 4.5)

        materials = openmc.Materials((water, fuel))
        materials.default_temperature = '294K'
        materials.export_to_xml()

        # Define the geometry.  Note that this geometry is somewhat non-sensical
        # (it essentially defines a circle of half-cylinders), but it is
        # designed so that periodic and reflective BCs will give different
        # answers.
        theta1 = (-1/6 + 1/2) * np.pi
        theta2 = (1/6 - 1/2) * np.pi
        plane1 = openmc.Plane(a=np.cos(theta1), b=np.sin(theta1),
                              boundary_type='periodic')
        plane2 = openmc.Plane(a=np.cos(theta2), b=np.sin(theta2),
                              boundary_type='periodic')

        x_max = openmc.XPlane(x0=5., boundary_type='reflective')

        z_cyl = openmc.ZCylinder(x0=3*np.cos(np.pi/6), y0=3*np.sin(np.pi/6),
                                 r=2.0)

        outside_cyl = openmc.Cell(1, fill=water, region=(
            +plane1 & +plane2 & -x_max & +z_cyl))
        inside_cyl = openmc.Cell(2, fill=fuel, region=(
            +plane1 & +plane2 & -z_cyl))
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
            (0, 0, 0), (5, 5, 0)))
        settings.export_to_xml()


def test_periodic():
    harness = Periodic6FoldTest('statepoint.4.h5')
    harness.main()
