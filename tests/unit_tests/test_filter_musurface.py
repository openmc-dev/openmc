import math

import numpy as np
import pytest
from uncertainties import unumpy

import openmc


def test_musurface(run_in_tmpdir):
    sphere = openmc.Sphere(r=1.0, boundary_type='vacuum')
    cell = openmc.Cell(fill=None, region=-sphere)
    model = openmc.Model()
    model.geometry = openmc.Geometry([cell])
    model.settings.particles = 10000
    model.settings.batches = 20
    E = 1.0
    source=openmc.IndependentSource()
    source.space = openmc.stats.Point()
    source.angle= openmc.stats.Isotropic()
    source.energy = openmc.stats.Discrete([E,], [1.0,]) 
    model.settings.source = source
    model.settings.run_mode = "fixed source"

    filter1 = openmc.MuSurfaceFilter(np.linspace(-1.0,1.0,201))
    filter2 = openmc.SurfaceFilter(sphere)
    tally = openmc.Tally()
    tally.filters = [filter1, filter2]
    tally.scores = ['current']

    model.tallies = openmc.Tallies([tally])

    # Run OpenMC
    sp_filename = model.run()

    # Get radial flux distribution
    with openmc.StatePoint(sp_filename) as sp:
        current_mu = sp.tallies[tally.id].mean.ravel()
    
    assert (current_mu[-1] ==  1.0)
    for element in current_mu[:-2]:
        assert (element == 0)


