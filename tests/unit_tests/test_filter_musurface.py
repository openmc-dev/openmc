import math

import numpy as np
import pytest
from uncertainties import unumpy

import openmc


def test_spherical_mesh_estimators(run_in_tmpdir):
    #"""Test that collision/tracklength estimators agree for SphericalMesh"""

    sphere = openmc.Sphere(r=10.0, boundary_type='vacuum')
    cell = openmc.Cell(fill=None, region=-sphere)
    model = openmc.Model()
    model.geometry = openmc.Geometry([cell])
    model.settings.particles = 1000
    model.settings.inactive = 10
    model.settings.batches = 20

    sph_mesh = openmc.MuSurfaceFilter(np.linspace(-1.0,1.0,201))

    tally = openmc.Tally()
    tally.filters = [openmc.MuSurfaceFilter(sphere)]
    tally.scores = ['current']

    model.tallies = openmc.Tallies([tally])

    # Run OpenMC
    sp_filename = model.run()

    # Get radial flux distribution
    with openmc.StatePoint(sp_filename) as sp:
        current_mu = sp.tallies[tally.id].mean.ravel()

    assert (current_mu[-1] ==  current_mu.max())
    assert (current_mu[-1] ==  1.0)
    assert (current_mu[0:-2] == 0)


