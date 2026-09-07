import openmc
import numpy as np
import pytest

@pytest.fixture
def one_group_lib():
    groups = openmc.mgxs.EnergyGroups([0.0, 20.0e6])
    xsdata = openmc.XSdata('slab_mat', groups)
    xsdata.order = 0
    xsdata.set_total([0.0])               
    xsdata.set_absorption([0.0])
    xsdata.set_scatter_matrix([[[0.0]]])
    
    mg_library = openmc.MGXSLibrary(groups)
    mg_library.add_xsdata(xsdata)
    name = 'mgxs.h5'
    mg_library.export_to_hdf5(name)
    yield name

@pytest.fixture
def slab_model(one_group_lib):
    model = openmc.Model()
    mat = openmc.Material(name='slab_material')
    mat.set_density('macro', 1.0)
    mat.add_macroscopic('slab_mat')

    model.materials = openmc.Materials([mat])
    model.materials.cross_sections = one_group_lib
    
    x_min = openmc.XPlane(x0=0.0, boundary_type='vacuum')
    x_max = openmc.XPlane(x0=10.0, boundary_type='vacuum')

    y_min = openmc.YPlane(y0=-10.0, boundary_type='vacuum')
    y_max = openmc.YPlane(y0=10.0, boundary_type='vacuum')
    z_min = openmc.ZPlane(z0=-10.0, boundary_type='vacuum')
    z_max = openmc.ZPlane(z0=19.0, boundary_type='vacuum')

    cell = openmc.Cell(fill=mat, region=+z_min & -x_max & +y_min & -y_max & +z_min & -z_max)
    model.geometry = openmc.Geometry([cell])

    model.settings = openmc.Settings()
    model.settings.energy_mode = 'multi-group'
    model.settings.run_mode = 'fixed source'
    model.settings.batches = 3
    model.settings.particles = 10

    source = openmc.IndependentSource()
    source.space = openmc.stats.Point((5.0, 0.0, 0.0))
    model.settings.source = source
    return model

def test_inverse_velocity(run_in_tmpdir, slab_model):
    tally = openmc.Tally()
    tally.scores = ['flux','inverse-velocity']
    slab_model.tallies = [tally]
    slab_model.run(apply_tally_results=True)
    inverse_velocity = tally.mean.squeeze()[1]/tally.mean.squeeze()[0]

    assert inverse_velocity == pytest.approx(1.6144e-5, rel=1e-4)
