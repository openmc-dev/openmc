import filecmp
from itertools import product
from pathlib import Path

import numpy as np

import openmc
import openmc.lib
import pytest

pytest.importorskip('vtk')


def ids_func(param):
    return f"{param['library']}_{param['elem_type']}"

test_params = (['libmesh', 'moab'],
               ['tets', 'hexes'])

test_cases = []
for library, elem_type in product(*test_params):
    test_case = {'library' : library,
                 'elem_type' : elem_type}
    test_cases.append(test_case)

@pytest.mark.parametrize("test_opts", test_cases, ids=ids_func)
def test_unstructured_mesh_to_vtk(run_in_tmpdir, request, test_opts):

    if test_opts['library'] == 'moab' and test_opts['elem_type'] == 'hexes':
        pytest.skip('Hexes are not supported with MOAB')

    if test_opts['library'] == 'libmesh' and not openmc.lib._libmesh_enabled():
        pytest.skip('LibMesh is not enabled in this build.')

    if test_opts['library'] == 'moab' and not openmc.lib._dagmc_enabled():
        pytest.skip('DAGMC (and MOAB) mesh not enabled in this build.')

    # pull in a simple model -- just need to create the statepoint file
    openmc.reset_auto_ids()
    model = openmc.examples.pwr_pin_cell()

    if test_opts['elem_type'] == 'tets':
        filename = Path('tets.exo')
    else:
        filename = Path('hexes.exo')

    # create a basic tally using the unstructured mesh
    umesh = openmc.UnstructuredMesh(request.node.path.parent / filename,
                                    test_opts['library'])
    umesh.output = False
    mesh_filter = openmc.MeshFilter(umesh)
    tally = openmc.Tally()
    tally.filters = [mesh_filter]
    tally.scores = ['flux']
    tally.estimator = 'collision'
    model.tallies = openmc.Tallies([tally])
    sp_file = model.run()

    # check VTK output after reading mesh from statepoint file
    with openmc.StatePoint(sp_file) as sp:
        umesh = sp.meshes[umesh.id]

    test_data = {'ids' : np.arange(umesh.n_elements)}
    umesh.write_data_to_vtk('umesh.vtk',
                            datasets=test_data,
                            volume_normalization=False)

    # compare file content with reference file
    ref_file = Path(f"{test_opts['library']}_{test_opts['elem_type']}_ref.vtk")
    assert filecmp.cmp('umesh.vtk', request.node.path.parent / ref_file)