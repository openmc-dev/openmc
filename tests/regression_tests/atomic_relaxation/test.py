from math import pi

import numpy as np
import openmc

from tests.testing_harness import PyAPITestHarness


class SourceTestHarness(PyAPITestHarness):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        mat = openmc.Material()
        mat.add_element('Pb', 1.0)
        mat.set_density('g/cm3', 11.35)
        self._model.materials = openmc.Materials([mat])

        sphere = openmc.Sphere(r=1.0e9, boundary_type='reflective')
        inside_sphere = openmc.Cell(fill=mat, region=-sphere)
        self._model.geometry = openmc.Geometry([inside_sphere])

        # Isotropic point source of 1 MeV photons at the origin
        source = openmc.IndependentSource()
        source.space = openmc.stats.Point((0, 0, 0))
        source.angle = openmc.stats.Isotropic()
        source.energy = openmc.stats.delta_function(1.0e6)
        source.particle = 'photon'

        # Fixed-source photon transport with atomic relaxation disabled
        settings = openmc.Settings()
        settings.particles = 10000
        settings.batches = 1
        settings.photon_transport = True
        settings.electron_treatment = 'led'
        settings.atomic_relaxation = False
        settings.run_mode = 'fixed source'
        settings.source = source
        self._model.settings = settings

        particle_filter = openmc.ParticleFilter(['photon', 'electron'])
        tally = openmc.Tally()
        tally.filters = [particle_filter]
        tally.scores = ['flux', 'heating']
        self._model.tallies = openmc.Tallies([tally])


def test_atomic_relaxation():
    harness = SourceTestHarness('statepoint.1.h5', model=openmc.Model())
    harness.main()
