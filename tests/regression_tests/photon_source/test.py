from math import pi

import numpy as np
import openmc

from tests.testing_harness import PyAPITestHarness


class SourceTestHarness(PyAPITestHarness):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        mat = openmc.Material()
        mat.set_density("g/cm3", 0.998207)
        mat.add_element("H", 0.111894)
        mat.add_element("O", 0.888106)
        self._model.materials = openmc.Materials([mat])

        sphere = openmc.Sphere(r=1.0e9, boundary_type="reflective")
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
        settings.cutoff = {"energy_photon": 1000.0}
        settings.run_mode = "fixed source"
        settings.source = source
        self._model.settings = settings

        particle_filter = openmc.ParticleFilter("photon")
        tally = openmc.Tally()
        tally.filters = [particle_filter]
        tally.scores = ["flux", "(n,gamma)"]
        tallies = openmc.Tallies([tally])
        self._model.tallies = tallies


def test_photon_source():
    harness = SourceTestHarness("statepoint.1.h5", model=openmc.Model())
    harness.main()
