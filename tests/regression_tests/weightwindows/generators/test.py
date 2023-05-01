import os

import numpy as np

import openmc

import pytest

def test_ww_generator(run_in_tmpdir):
    # create a simple spherical shield model
    model = openmc.Model()

    water = openmc.Material()
    water.set_density('g/cc', 1.0)
    water.add_nuclide('H1', 0.66)
    water.add_nuclide('O16', 0.34)

    model.materials = openmc.Materials([water])

    s = openmc.Sphere(r=50, boundary_type='vacuum')
    c = openmc.Cell(fill=water, region=-s)

    model.geometry = openmc.Geometry([c])

    model.settings.particles = 500
    model.settings.batches = 5
    model.settings.run_mode = 'fixed source'
    model.settings.max_splits = 100

    mesh = openmc.RegularMesh.from_domain(model.geometry.root_universe)
    energy_bounds = np.linspace(0, 1E6, 70)
    particle = 'neutron'

    wwg = openmc.WeightWindowGenerator(mesh, energy_bounds, particle)
    wwg.max_realizations = 1
    wwg.on_the_fly = True
    wwg.update_params = {'ratio' : 5.0}

    model.settings.weight_window_generators = wwg
    model.export_to_xml()

    model_in = openmc.Model.from_xml()
    assert len(model_in.settings.weight_window_generators) == 1
    wwg_in = model_in.settings.weight_window_generators[0]
    assert wwg_in.max_realizations == 1
    assert wwg_in.on_the_fly == True
    assert wwg_in.update_interval == 1
    assert wwg_in.update_params['ratio'] == 5.0

    openmc.run()
    assert os.path.exists('weight_windows.h5')