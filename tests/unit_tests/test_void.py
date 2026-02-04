import numpy as np
import openmc
import pytest


@pytest.fixture
def model():
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

def test_equivalt_void_specification(model, run_in_tmpdir):
    sp_file = model.run()
    with openmc.StatePoint(sp_file) as sp:
        tally1 = sp.tallies[1]
        tally2 = sp.tallies[2]
        
    assert np.isclose(tally1.mean.sum(),tally2.mean.sum(), rtol=1e-10, atol=0)
