#!/usr/bin/env python

from math import pi
import os
import sys

import numpy as np

sys.path.insert(0, os.pardir)
from testing_harness import PyAPITestHarness
import openmc


class SourceTestHarness(PyAPITestHarness):
    def _build_inputs(self):
        mat1 = openmc.Material(material_id=1)
        mat1.set_density('g/cm3', 4.5)
        mat1.add_nuclide(openmc.Nuclide('U-235', '71c'), 1.0)
        materials = openmc.Materials([mat1])
        materials.export_to_xml()

        sphere = openmc.Sphere(surface_id=1, R=10.0, boundary_type='vacuum')
        inside_sphere = openmc.Cell(cell_id=1)
        inside_sphere.region = -sphere
        inside_sphere.fill = mat1

        root = openmc.Universe(universe_id=0)
        root.add_cell(inside_sphere)
        geometry = openmc.Geometry()
        geometry.root_universe = root
        geometry.export_to_xml()

        # Create an array of different sources
        x_dist = openmc.stats.Uniform(-3., 3.)
        y_dist = openmc.stats.Discrete([-4., -1., 3.], [0.2, 0.3, 0.5])
        z_dist = openmc.stats.Tabular([-2., 0., 2.], [0.2, 0.3, 0.2])
        spatial1 = openmc.stats.CartesianIndependent(x_dist, y_dist, z_dist)
        spatial2 = openmc.stats.Box([-4., -4., -4.], [4., 4., 4.])
        spatial3 = openmc.stats.Point([1.2, -2.3, 0.781])

        mu_dist = openmc.stats.Discrete([-1., 0., 1.], [0.5, 0.25, 0.25])
        phi_dist = openmc.stats.Uniform(0., 6.28318530718)
        angle1 = openmc.stats.PolarAzimuthal(mu_dist, phi_dist)
        angle2 = openmc.stats.Monodirectional(reference_uvw=[0., 1., 0.])
        angle3 = openmc.stats.Isotropic()

        E = np.logspace(-6, 1)
        p = np.sin(np.linspace(0., pi))
        p /= sum(np.diff(E)*p[:-1])
        energy1 = openmc.stats.Maxwell(1.2895)
        energy2 = openmc.stats.Watt(0.988, 2.249)
        energy3 = openmc.stats.Tabular(E, p, interpolation='histogram')

        source1 = openmc.Source(spatial1, angle1, energy1, strength=0.5)
        source2 = openmc.Source(spatial2, angle2, energy2, strength=0.3)
        source3 = openmc.Source(spatial3, angle3, energy3, strength=0.2)

        settings = openmc.Settings()
        settings.batches = 10
        settings.inactive = 5
        settings.particles = 1000
        settings.source = [source1, source2, source3]
        settings.export_to_xml()


if __name__ == '__main__':
    harness = SourceTestHarness('statepoint.10.h5')
    harness.main()
