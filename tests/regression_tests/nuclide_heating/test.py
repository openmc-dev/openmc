from math import pi

import numpy as np
import openmc

from tests.testing_harness import PyAPITestHarness


class NuclideHeatingTestHarness(PyAPITestHarness):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        mat = openmc.Material()
        mat.add_nuclide("Cu63", 0.5, "ao")
        mat.add_nuclide("Cu65", 0.5, "ao")
        mat.set_density("g/cm3", 1.0)
        self._model.materials = openmc.Materials([mat])

        sphere = openmc.Sphere(r=20, boundary_type="reflective")
        inside_sphere = openmc.Cell()
        inside_sphere.region = -sphere
        inside_sphere.fill = mat
        self._model.geometry = openmc.Geometry([inside_sphere])

        source = openmc.IndependentSource()
        source.space = openmc.stats.Point((0, 0, 0))
        source.angle = openmc.stats.Isotropic()
        source.energy = openmc.stats.Discrete([10.0e6], [1.0])
        source.particle = "photon"

        settings = openmc.Settings()
        settings.particles = 10000
        settings.batches = 1
        settings.photon_transport = True
        settings.electron_treatment = "ttb"
        settings.cutoff = {"energy_photon": 1000}
        settings.run_mode = "fixed source"
        settings.source = source
        self._model.settings = settings

        tally1 = openmc.Tally()
        tally1.scores = ["heating"]
        tally1.nuclides = ["Cu63", "Cu65"]

        tally2 = openmc.Tally()
        tally2.scores = ["heating"]

        tallies = openmc.Tallies([tally1, tally2])
        self._model.tallies = tallies


def test_nuclide_heating():
    harness = NuclideHeatingTestHarness("statepoint.1.h5", model=openmc.Model())
    harness.main()
