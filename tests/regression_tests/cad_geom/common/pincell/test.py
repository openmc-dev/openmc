from pathlib import Path

import openmc
import openmc.lib

import h5py
import numpy as np
import pytest

from tests.testing_harness import PyAPITestHarness, config

pytestmark = pytest.mark.skipif(
    not openmc.lib.feature_enabled('dagmc'),
    reason="DAGMC CAD geometry is not enabled.")


test_libs = []
if openmc.lib.feature_enabled('dagmc'):
    test_libs.append('dagmc')
if openmc.lib.feature_enabled('xdg'):
    test_libs.append('xdg')

@pytest.mark.parametrize('interface', test_libs)
def test_missing_material_id(cad_pincell_model, interface):
    model = cad_pincell_model(interface)
    # remove the last material, which is identified by ID in the DAGMC file
    model.materials = model.materials[:-1]
    with pytest.raises(RuntimeError) as exec_info:
        model.run()
    exp_error_msg = "Material with name/ID '41' not found for volume (cell) 3"
    assert exp_error_msg in str(exec_info.value)


@pytest.mark.parametrize('interface', test_libs)
def test_missing_material_name(cad_pincell_model, interface):
    model = cad_pincell_model(interface)
    # remove the first material, which is identified by name in the DAGMC file
    model.materials = model.materials[1:]
    with pytest.raises(RuntimeError) as exec_info:
        model.run()
    exp_error_msg = "Material with name/ID 'no-void fuel' not found for volume (cell) 1"
    assert exp_error_msg in str(exec_info.value)

@pytest.mark.parametrize('interface', test_libs)
def test_surf_source(cad_pincell_model, interface):
    model = cad_pincell_model()
    # create a surface source read on this model to ensure
    # particles are being generated correctly
    n = 100
    model.settings.surf_source_write = {'surface_ids': [1], 'max_particles': n}

    # If running in MPI mode, setup proper keyword arguments for run()
    kwargs = {'openmc_exec': config['exe']}
    if config['mpi']:
        kwargs['mpi_args'] = [config['mpiexec'], '-n', config['mpi_np']]
    model.run(**kwargs)

    with h5py.File('surface_source.h5') as fh:
        assert fh.attrs['filetype'] == b'source'
        arr = fh['source_bank'][...]
    expected_size = n * int(config['mpi_np']) if config['mpi'] else n
    assert arr.size == expected_size

    # check that all particles are on surface 1 (radius = 7)
    xs = arr[:]['r']['x']
    ys = arr[:]['r']['y']
    rad = np.sqrt(xs**2 + ys**2)
    assert np.allclose(rad, 7.0)


@pytest.mark.parametrize("interface", test_libs)
def test_transport(cad_pincell_model, interface):
    input_file = f'inputs_{interface}.dat'
    model = cad_pincell_model(interface)
    harness = PyAPITestHarness('statepoint.5.h5', model, input_file)
    harness.main()