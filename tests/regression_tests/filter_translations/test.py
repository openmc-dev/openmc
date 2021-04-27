from itertools import permutations
import numpy as np

import openmc
import pytest

from tests.testing_harness import HashedPyAPITestHarness


translations = list(permutations((10, 5, 0))) + list(permutations((-10, -5, 0)))
@pytest.fixture(params=translations)
def model(request):

    translation = np.array(request.param)

    model = openmc.model.Model()

    fuel = openmc.Material()
    fuel.set_density('g/cm3', 10.0)
    fuel.add_nuclide('U235', 1.0)
    zr = openmc.Material()
    zr.set_density('g/cm3', 1.0)
    zr.add_nuclide('Zr90', 1.0)
    model.materials.extend([fuel, zr])

    box1 = openmc.model.rectangular_prism(10.0, 10.0)
    box2 = openmc.model.rectangular_prism(20.0, 20.0, boundary_type='reflective')
    top = openmc.ZPlane(z0=10.0, boundary_type='vacuum')
    bottom = openmc.ZPlane(z0=-10.0, boundary_type='vacuum')
    cell1 = openmc.Cell(fill=fuel, region=box1 & +bottom & -top)
    cell2 = openmc.Cell(fill=zr, region=~box1 & box2 & +bottom & -top)
    model.geometry = openmc.Geometry([cell1, cell2])

    model.settings.batches = 5
    model.settings.inactive = 0
    model.settings.particles = 1000

    llc = np.array([-9, -9, -9]) - translation
    urc = np.array([9, 9, 9]) - translation

    reg_mesh = openmc.RegularMesh()
    reg_mesh.dimension = [9, 9, 9]
    reg_mesh.lower_left = llc
    reg_mesh.upper_right = urc

    recti_mesh = openmc.RectilinearMesh()
    recti_mesh.x_grid = np.linspace(llc[0], urc[0], 17)
    recti_mesh.y_grid = np.linspace(llc[1], urc[1], 17)
    recti_mesh.z_grid = np.linspace(llc[2], urc[2], 17)

    # Create filters
    filters = [openmc.MeshFilter(reg_mesh),
               openmc.MeshFilter(recti_mesh)]
    for f in filters:
        f.translation = translation

    # Create tallies
    for f in filters:
        tally = openmc.Tally()
        tally.filters = [f]
        tally.scores = ['total']
        model.tallies.append(tally)

    input_name = 'input_dat_' + '_'.join(map(str, translation)) + '.dat'
    return model, input_name


def test_filter_mesh_translations(model):
    harness = HashedPyAPITestHarness('statepoint.5.h5', *model)
    harness.main()


