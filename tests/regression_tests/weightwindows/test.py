from copy import deepcopy
import pytest
import numpy as np

import openmc
from openmc.stats import Discrete, Point

from tests.testing_harness import HashedPyAPITestHarness


@pytest.fixture
def model():
    model = openmc.Model()

    # materials (M4 steel alloy)
    m4 = openmc.Material()
    m4.set_density('g/cc', 2.3)
    m4.add_nuclide('H1', 0.168018676)
    m4.add_nuclide("H2", 1.93244e-05)
    m4.add_nuclide("O16", 0.561814465)
    m4.add_nuclide("O17", 0.00021401)
    m4.add_nuclide("Na23", 0.021365)
    m4.add_nuclide("Al27", 0.021343)
    m4.add_nuclide("Si28", 0.187439342)
    m4.add_nuclide("Si29", 0.009517714)
    m4.add_nuclide("Si30", 0.006273944)
    m4.add_nuclide("Ca40", 0.018026179)
    m4.add_nuclide("Ca42", 0.00012031)
    m4.add_nuclide("Ca43", 2.51033e-05)
    m4.add_nuclide("Ca44", 0.000387892)
    m4.add_nuclide("Ca46", 7.438e-07)
    m4.add_nuclide("Ca48", 3.47727e-05)
    m4.add_nuclide("Fe54", 0.000248179)
    m4.add_nuclide("Fe56", 0.003895875)
    m4.add_nuclide("Fe57", 8.99727e-05)
    m4.add_nuclide("Fe58", 1.19737e-05)

    s0 = openmc.Sphere(r=240)
    s1 = openmc.Sphere(r=250, boundary_type='vacuum')

    c0 = openmc.Cell(fill=m4, region=-s0)
    c1 = openmc.Cell(region=+s0 & -s1)

    model.geometry = openmc.Geometry([c0, c1])

    # settings
    settings = model.settings
    settings.run_mode = 'fixed source'
    settings.particles = 200
    settings.batches = 2
    settings.max_splits = 200
    settings.photon_transport = True
    space = Point((0.001, 0.001, 0.001))
    energy = Discrete([14E6], [1.0])

    settings.source = openmc.Source(space=space, energy=energy)

    # tally
    mesh = openmc.RegularMesh()
    mesh.lower_left = (-240, -240, -240)
    mesh.upper_right = (240, 240, 240)
    mesh.dimension = (5, 10, 15)

    mesh_filter = openmc.MeshFilter(mesh)

    e_bnds = [0.0, 0.5, 2E7]
    energy_filter = openmc.EnergyFilter(e_bnds)

    particle_filter = openmc.ParticleFilter(['neutron', 'photon'])

    tally = openmc.Tally()
    tally.filters = [mesh_filter, energy_filter, particle_filter]
    tally.scores = ['flux']

    model.tallies.append(tally)

    # weight windows

    # load pre-generated weight windows
    # (created using the same tally as above)
    ww_n_lower_bnds = np.loadtxt('ww_n.txt')
    ww_p_lower_bnds = np.loadtxt('ww_p.txt')

    # create a mesh matching the one used
    # to generate the weight windows
    ww_mesh = openmc.RegularMesh()
    ww_mesh.lower_left = (-240, -240, -240)
    ww_mesh.upper_right = (240, 240, 240)
    ww_mesh.dimension = (5, 6, 7)

    ww_n = openmc.WeightWindows(ww_mesh,
                                ww_n_lower_bnds,
                                None,
                                10.0,
                                e_bnds,
                                max_lower_bound_ratio=1.5)

    ww_p = openmc.WeightWindows(ww_mesh,
                                ww_p_lower_bnds,
                                None,
                                10.0,
                                e_bnds,
                                max_lower_bound_ratio=1.5)

    model.settings.weight_windows = [ww_n, ww_p]

    return model


def test_weightwindows(model):
    test = HashedPyAPITestHarness('statepoint.2.h5', model)
    test.main()
