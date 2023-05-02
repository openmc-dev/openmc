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
    meshes = {mesh.id : mesh}

    wwg = openmc.WeightWindowGenerator(mesh, energy_bounds, particle)
    wwg.max_realizations = 1
    wwg.on_the_fly = True
    wwg.update_params = {'ratio' : 5.0,
                         'threshold': 0.8,
                         'value' : 'mean'}

    model.settings.weight_window_generators = wwg
    model.export_to_xml()

    openmc.run()
    # we test the effectiveness of the update ethod elsewhere, so
    # just test that the generation happens successfully here
    assert os.path.exists('weight_windows.h5')

    wws_mean = openmc.hdf5_to_wws(meshes=meshes)
    assert len(wws_mean) == 1

    # check that generation using the relative error works too
    wwg.update_params['value'] = 'rel_err'
    model.run()

    wws_rel_err = openmc.hdf5_to_wws(meshes=meshes)
    assert len(wws_rel_err) == 1

    # we should not get the same set of weight windows when switching to use of
    # rel. err.
    assert (wws_mean[0].lower_ww_bounds != wws_rel_err[0].lower_ww_bounds).any()
