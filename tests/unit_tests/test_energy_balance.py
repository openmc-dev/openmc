import openmc
import pytest
import numpy as np


@pytest.fixture
def two_cell_model():
    """Simple two-cell slab model with a fixed source.
       Cell1 occupies x in [-10, 0], cell2 x in [0, 10]. 
    """
    openmc.reset_auto_ids()
    model = openmc.Model()

    m1 = openmc.Material()
    m1.add_element('Li', 1.0)
    m1.set_density('g/cm3', 1.0)

    m2 = openmc.Material()
    m2.add_element('Al', 1.0)
    m2.set_density('g/cm3', 1.0)
    
    xmin = openmc.XPlane(-10.0, boundary_type="reflective")
    xmid = openmc.XPlane(0.0)
    xmax = openmc.XPlane(10.0, boundary_type="reflective")
    ymin = openmc.YPlane(-10.0, boundary_type="reflective")
    ymax = openmc.YPlane(10.0, boundary_type="reflective")
    zmin = openmc.ZPlane(-10.0, boundary_type="reflective")
    zmax = openmc.ZPlane(10.0, boundary_type="reflective")
    cell1 = openmc.Cell(fill=m1, region=+xmin & -xmid & +ymin & -ymax & +zmin & -zmax)
    cell2 = openmc.Cell(fill=m2, region=+xmid & -xmax & +ymin & -ymax & +zmin & -zmax)
    model.geometry = openmc.Geometry([cell1, cell2])

    src = openmc.IndependentSource()
    src.space = openmc.stats.Point((-5.0, 0.0, 0.0))
    src.particle = 'photon'
    src.energy = openmc.stats.Discrete([1e6],[1.0])
    
    model.settings.run_mode = 'fixed source'
    model.settings.batches = 1
    model.settings.particles = 100
    model.settings.source = src

    return model, xmid, cell1, cell2, m1, m2

@pytest.fixture
def two_sphere_model():
    openmc.reset_auto_ids()
    model = openmc.Model()

    mat = openmc.Material()
    mat.add_nuclide('H2', 1.0)
    mat.set_density('g/cm3', 1.0)

    surf1 = openmc.Sphere(r=0.01)
    surf2 = openmc.Sphere(r=1000, boundary_type='vacuum')
    cell1 = openmc.Cell(region=-surf1)
    cell2 = openmc.Cell(fill=mat, region=+surf1 & -surf2)

    model.geometry = openmc.Geometry([cell1, cell2])

    src = openmc.IndependentSource()
    src.energy = openmc.stats.Discrete([1e6],[1.0])

    model.settings.run_mode = 'fixed source'
    model.settings.batches = 1
    model.settings.particles = 1000
    model.settings.source = src

    return model, cell1, cell2, mat

def test_energy_balance_simple(two_cell_model, run_in_tmpdir):
    model, xmid, cell1, cell2, *_  = two_cell_model

    tally1 = openmc.Tally()
    tally1.scores = ['heating']
    tally1.filters = [openmc.CellFilter([cell1, cell2])]

    tally2 = openmc.Tally()
    tally2.scores = ['current']
    tally2.filters = [openmc.EnergyFunctionFilter([0.0,20e6],[0.0,20e6]), openmc.SurfaceFilter([xmid])]

    
    model.tallies = [tally1, tally2]
    model.run(apply_tally_results=True)

    assert tally1.mean.sum() == pytest.approx(1e6)

    assert tally2.mean[0] == pytest.approx(tally1.mean[1])

def test_current_conservation(two_cell_model, run_in_tmpdir):
    model, xmid, cell1, cell2, m1, m2 = two_cell_model

    energyfunc_filter = openmc.EnergyFunctionFilter([0.0,20e6],[0.0,20e6])

    tally1 = openmc.Tally()
    tally1.scores = ['current']
    tally1.filters = [energyfunc_filter, openmc.SurfaceFilter([xmid])]

    tally2 = openmc.Tally()
    tally2.scores = ['current']
    tally2.filters = [energyfunc_filter,
                      openmc.CellFromFilter([cell1]),
                      openmc.CellFilter([cell2])]

    tally3 = openmc.Tally()
    tally3.scores = ['current']
    tally3.filters = [energyfunc_filter,
                      openmc.CellFromFilter([cell2]),
                      openmc.CellFilter([cell1])]

    tally4 = openmc.Tally()
    tally4.scores = ['current']
    tally4.filters = [energyfunc_filter,
                      openmc.MaterialFromFilter([m1]),
                      openmc.MaterialFilter([m2])]

    tally5 = openmc.Tally()
    tally5.scores = ['current']
    tally5.filters = [energyfunc_filter,
                      openmc.MaterialFromFilter([m2]),
                      openmc.MaterialFilter([m1])]
    
    
    model.tallies = [tally1, tally2, tally3, tally4, tally5]
    model.run(apply_tally_results=True)

    assert pytest.approx(tally1.mean.sum()) == tally2.mean[0] + tally3.mean[0]

    assert tally2.mean[0] == pytest.approx(tally4.mean[0])

    assert tally3.mean[0] == pytest.approx(tally5.mean[0])

def test_cellfrom_heating(run_in_tmpdir, two_sphere_model):

    model, cell1, cell2, mat = two_sphere_model

    tally1 = openmc.Tally()
    tally1.scores = ['heating']
    tally1.filters = [openmc.CellFromFilter([cell1, cell2]), openmc.CellFilter([cell2])]

    tally2 = openmc.Tally()
    tally2.scores = ['heating']
    tally2.filters = [openmc.MaterialFromFilter([mat]), openmc.MaterialFilter([mat])]

    model.tallies = [tally1, tally2]

    model.run(apply_tally_results=True)

    assert np.all(tally1.mean > 0)

    assert tally1.mean[1] == tally2.mean[0] 
    
