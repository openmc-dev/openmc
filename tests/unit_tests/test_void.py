import numpy as np
import openmc
import pytest


@pytest.fixture
def model():
    openmc.reset_auto_ids()
    model = openmc.model.Model()

    zn = openmc.Material()
    zn.set_density('g/cm3', 7.14)
    zn.add_nuclide('Zn64', 1.0)
    model.materials.append(zn)

    radii = np.linspace(1.0, 100.0)
    surfs = [openmc.Sphere(r=r) for r in radii]
    surfs[-1].boundary_type = 'vacuum'
    cells = [openmc.Cell(fill=(zn if i % 2 == 0 else None), region=region)
             for i, region in enumerate(openmc.model.subdivide(surfs))]
    model.geometry = openmc.Geometry(cells)

    model.settings.run_mode = 'fixed source'
    model.settings.batches = 3
    model.settings.particles = 1000
    model.settings.source = openmc.IndependentSource(space=openmc.stats.Point())

    cell_filter = openmc.CellFilter(cells[1::2])
    material_filter = openmc.MaterialFilter([None])
    
    tally1 = openmc.Tally()
    tally1.filters = [cell_filter]
    tally1.scores = ['flux']
    
    tally2 = openmc.Tally()
    tally2.filters = [material_filter]
    tally2.scores = ['flux']
    
    model.tallies.append(tally1)
    model.tallies.append(tally2)

    return model

@pytest.fixture
def empty_sphere():
    openmc.reset_auto_ids()
    model = openmc.Model()
    surf = openmc.Sphere(r=10, boundary_type='vacuum')
    cell = openmc.Cell(region=-surf)
    model.geometry = openmc.Geometry([cell])

    model.settings.run_mode = 'fixed source'
    model.settings.batches = 3
    model.settings.particles = 1000

    tally = openmc.Tally()
    tally.scores = ['total', 'elastic']
    tally.nuclides = ['U235']
    tally.multiply_density = False
    model.tallies.append(tally)

    return model


def test_equivalent_microxs(empty_sphere, run_in_tmpdir):
    sp_file = empty_sphere.run()
    with openmc.StatePoint(sp_file) as sp:
        tally1 = sp.tallies[1]

    mat = openmc.Material()
    mat.add_nuclide('H1', 1e-16)

    empty_sphere.geometry.get_all_cells()[1].fill = mat

    sp_file = empty_sphere.run()
    with openmc.StatePoint(sp_file) as sp:
        tally2 = sp.tallies[1]

    assert np.isclose(tally1.mean.sum(), tally2.mean.sum(), rtol=1e-10, atol=0)
    assert tally1.mean.sum() > 0


def test_equivalent_void_specification(model, run_in_tmpdir):
    sp_file = model.run()
    with openmc.StatePoint(sp_file) as sp:
        tally1 = sp.tallies[1]
        tally2 = sp.tallies[2]
        
    assert np.isclose(tally1.mean.sum(),tally2.mean.sum(), rtol=1e-10, atol=0)
