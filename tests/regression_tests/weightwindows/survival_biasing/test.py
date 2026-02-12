import pytest
import numpy as np

import openmc
from openmc.stats import Discrete, Point

from tests.testing_harness import HashedPyAPITestHarness


@pytest.fixture
def model():
    # Material
    w = openmc.Material(name='Tungsten')
    w.add_element('W', 1.0)
    w.set_density('g/cm3', 19.25)

    materials = openmc.Materials([w])

    # Geometry surfaces
    x0 = openmc.XPlane(x0=0.0, boundary_type='reflective')
    x1 = openmc.XPlane(x0=160.0, boundary_type='vacuum')
    y0 = openmc.YPlane(y0=0.0, boundary_type='reflective')
    y1 = openmc.YPlane(y0=160.0, boundary_type='reflective')
    z0 = openmc.ZPlane(z0=0.0, boundary_type='reflective')
    z1 = openmc.ZPlane(z0=160.0, boundary_type='reflective')

    region = +x0 & -x1 & +y0 & -y1 & +z0 & -z1
    cell = openmc.Cell(region=region, fill=w)
    root = openmc.Universe(cells=[cell])
    geometry = openmc.Geometry(root)

    # Source: planar on x=0, mono-directional along +x, 14.1 MeV neutrons
    space = openmc.stats.CartesianIndependent(
        openmc.stats.Discrete([0.01], [1.0]),
        openmc.stats.Uniform(0.0, 160.0),
        openmc.stats.Uniform(0.0, 160.0),
    )
    angle = openmc.stats.Monodirectional((1.0, 0.0, 0.0))
    energy = openmc.stats.Discrete([14.1e6], [1.0])

    source = openmc.Source(space=space, angle=angle, energy=energy)

    settings = openmc.Settings()
    settings.run_mode = 'fixed source'
    settings.batches = 5
    settings.particles = 50
    settings.source = source

    model = openmc.Model(geometry=geometry, materials=materials, settings=settings)

    # Mesh tally: 1 cm voxels, flux only
    mesh = openmc.RegularMesh()
    mesh.dimension = (20, 20, 1)
    mesh.lower_left = (0.0, 0.0, 0.0)
    mesh.upper_right = (160.0, 160.0, 160.0)

    mesh_filter = openmc.MeshFilter(mesh)
    flux_tally = openmc.Tally(name='flux')
    flux_tally.filters = [mesh_filter]
    flux_tally.scores = ['flux']
    tallies = openmc.Tallies([flux_tally])
    model.tallies = tallies

    lower_ww_bounds = np.loadtxt('ww_n.txt')

    weight_windows = openmc.WeightWindows(mesh,
                                          lower_ww_bounds,
                                          upper_bound_ratio=5.0,
                                          particle_type='neutron')

    model.settings.weight_windows = weight_windows
    model.settings.weight_window_checkpoints = {'surface': True,
                                                'collision': True}
    model.settings.survival_biasing = True

    return model


def test_weight_windows_with_survival_biasing(model):
    harness = HashedPyAPITestHarness('statepoint.5.h5', model)
    harness.main()
