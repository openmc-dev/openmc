from math import pi, cos

import numpy as np
import openmc

from tests.testing_harness import PyAPITestHarness


class SourceTestHarness(PyAPITestHarness):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        mat1 = openmc.Material(material_id=1, temperature=294)
        mat1.set_density("g/cm3", 4.5)
        mat1.add_nuclide(openmc.Nuclide("U235"), 1.0)
        self._model.materials = openmc.Materials([mat1])

        sphere = openmc.Sphere(surface_id=1, r=10.0, boundary_type="vacuum")
        inside_sphere = openmc.Cell(cell_id=1)
        inside_sphere.region = -sphere
        inside_sphere.fill = mat1

        root = openmc.Universe(universe_id=0)
        root.add_cell(inside_sphere)
        self._model.geometry = openmc.Geometry(root)

        # Create an array of different sources
        x_dist = openmc.stats.Uniform(-3.0, 3.0)
        y_dist = openmc.stats.Discrete([-4.0, -1.0, 3.0], [0.2, 0.3, 0.5])
        z_dist = openmc.stats.Tabular([-2.0, 0.0, 2.0], [0.2, 0.3, 0.2])
        r_dist = openmc.stats.Uniform(2.0, 3.0)
        r_dist1 = openmc.stats.PowerLaw(2.0, 3.0, 1.0)
        r_dist2 = openmc.stats.PowerLaw(2.0, 3.0, 2.0)
        cos_theta_dist = openmc.stats.Discrete(
            [cos(pi / 4), 0.0, cos(3 * pi / 4)], [0.3, 0.4, 0.3]
        )
        phi_dist = openmc.stats.Uniform(0.0, 2 * pi)
        spatial1 = openmc.stats.CartesianIndependent(x_dist, y_dist, z_dist)
        spatial2 = openmc.stats.Box([-4.0, -4.0, -4.0], [4.0, 4.0, 4.0])
        spatial3 = openmc.stats.Point([1.2, -2.3, 0.781])
        spatial4 = openmc.stats.SphericalIndependent(
            r_dist, cos_theta_dist, phi_dist, origin=(1.0, 1.0, 0.0)
        )
        spatial5 = openmc.stats.CylindricalIndependent(
            r_dist, phi_dist, z_dist, origin=(1.0, 1.0, 0.0)
        )
        spatial6 = openmc.stats.SphericalIndependent(
            r_dist2, cos_theta_dist, phi_dist, origin=(1.0, 1.0, 0.0)
        )
        spatial7 = openmc.stats.CylindricalIndependent(
            r_dist1, phi_dist, z_dist, origin=(1.0, 1.0, 0.0)
        )

        mu_dist = openmc.stats.Discrete([-1.0, 0.0, 1.0], [0.5, 0.25, 0.25])
        phi_dist = openmc.stats.Uniform(0.0, 6.28318530718)
        angle1 = openmc.stats.PolarAzimuthal(mu_dist, phi_dist)
        angle2 = openmc.stats.Monodirectional(reference_uvw=[0.0, 1.0, 0.0])
        angle3 = openmc.stats.Isotropic()

        # Note that the definition for E is equivalent to logspace(0, 7) but we
        # manually take powers because of last-digit differences that may cause
        # test failures with different versions of numpy
        E = np.array([10**x for x in np.linspace(0, 7)])
        p = np.sin(np.linspace(0.0, pi))
        p /= sum(np.diff(E) * p[:-1])
        energy1 = openmc.stats.Maxwell(1.2895e6)
        energy2 = openmc.stats.Watt(0.988e6, 2.249e-6)
        energy3 = openmc.stats.Tabular(E, p, interpolation="histogram")
        energy4 = openmc.stats.Mixture([1, 2, 3], [energy1, energy2, energy3])

        time1 = openmc.stats.Uniform(2, 5)

        source1 = openmc.IndependentSource(spatial1, angle1, energy1, strength=0.3)
        source2 = openmc.IndependentSource(spatial2, angle2, energy2, strength=0.1)
        source3 = openmc.IndependentSource(spatial3, angle3, energy3, strength=0.1)
        source4 = openmc.IndependentSource(spatial4, angle3, energy3, strength=0.1)
        source5 = openmc.IndependentSource(spatial5, angle3, energy3, strength=0.1)
        source6 = openmc.IndependentSource(spatial5, angle3, energy4, strength=0.1)
        source7 = openmc.IndependentSource(
            spatial6, angle3, energy4, time1, strength=0.1
        )
        source8 = openmc.IndependentSource(
            spatial7, angle3, energy4, time1, strength=0.1
        )

        settings = openmc.Settings()
        settings.batches = 10
        settings.inactive = 5
        settings.particles = 1000
        settings.source = [
            source1,
            source2,
            source3,
            source4,
            source5,
            source6,
            source7,
            source8,
        ]
        self._model.settings = settings


def test_source():
    harness = SourceTestHarness("statepoint.10.h5", model=openmc.Model())
    harness.main()
