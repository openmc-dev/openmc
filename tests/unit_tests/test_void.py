import numpy as np
import openmc
import pytest


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
