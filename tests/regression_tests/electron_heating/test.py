import os
import glob

import numpy as np
import pytest

import openmc

from tests.testing_harness import PyAPITestHarness


class ElectronHeatingTest(PyAPITestHarness):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Define materials
        water = openmc.Material(1)
        water.add_nuclide("H1", 2.0)
        water.add_nuclide("O16", 1.0)
        water.set_density("g/cc", 1.0)

        materials = openmc.Materials((water,))
        self._model.materials = materials

        sphere = openmc.Sphere(surface_id=1, r=1.0, boundary_type="reflective")

        # Define geometry
        sph = openmc.Cell(1, fill=water, region=-sphere)
        root = openmc.Universe(0, cells=(sph,))

        self._model.geometry = openmc.Geometry(root)

        # Define source

        source = openmc.IndependentSource()
        source.space = openmc.stats.Point((0, 0, 0))
        source.angle = openmc.stats.Isotropic()
        source.energy = openmc.stats.Discrete([10.0e6], [1.0])
        source.particle = "electron"

        # Define settings
        settings = openmc.Settings()
        settings.particles = 10000
        settings.batches = 1
        settings.photon_transport = True
        settings.electron_treatment = "ttb"
        settings.cutoff = {"energy_photon": 1000.0}
        settings.run_mode = "fixed source"
        settings.source = source
        self._model.settings = settings

        # Define tallies

        tally = openmc.Tally()
        tally.scores = ["heating"]
        tallies = openmc.Tallies([tally])
        self._model.tallies = tallies


def test_electron_heating_calc():
    harness = ElectronHeatingTest("statepoint.1.h5", model=openmc.Model())
    harness.main()
