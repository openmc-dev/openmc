import h5py
import numpy as np
import pytest

import openmc
import openmc.lib


def test_ww_roundtrip(request, run_in_tmpdir, monkeypatch):
    # Load weight windows from a wwinp file
    wwinp_file = request.path.with_name('wwinp_n')
    wws = openmc.WeightWindowsList.from_wwinp(wwinp_file)

    def fail_xml_export(*args, **kwargs):
        pytest.fail('Weight windows should not be serialized to XML')

    monkeypatch.setattr(openmc.WeightWindows, 'to_xml_element', fail_xml_export)

    # Roundtrip them, writing to HDF5 and reading back in
    wws.export_to_hdf5('ww.h5')
    wws_new = openmc.WeightWindowsList.from_hdf5('ww.h5')

    # Check that the new weight windows are the same as the original
    assert len(wws) == len(wws_new)
    for ww, ww_new in zip(wws, wws_new):
        assert ww.particle_type == ww_new.particle_type
        assert (ww.lower_ww_bounds == ww_new.lower_ww_bounds).all()
        assert (ww.upper_ww_bounds == ww_new.upper_ww_bounds).all()
        assert ww.survival_ratio == ww_new.survival_ratio
        assert ww.num_energy_bins == ww_new.num_energy_bins
        assert ww.max_split == ww_new.max_split
        assert ww.weight_cutoff == ww_new.weight_cutoff
        assert ww.mesh.id == ww_new.mesh.id


def test_export_hdf5_format(request, run_in_tmpdir):
    # openmc_weight_windows_import expects this on-disk layout.
    wws = openmc.WeightWindowsList.from_wwinp(request.path.with_name('wwinp_n'))
    wws.export_to_hdf5('ww.h5')

    with h5py.File('ww.h5') as f:
        assert f.attrs['filetype'] == b'weight_windows'
        assert list(f.attrs['version']) == [1, 0]
        wws_group = f['weight_windows']
        assert int(wws_group.attrs['n_weight_windows']) == len(wws)
        for ww in wws:
            group = wws_group[f'weight_windows_{ww.id}']
            assert group['lower_ww_bounds'].ndim == 2
            assert group['lower_ww_bounds'].shape[0] == ww.num_energy_bins
            assert 'max_lower_bound_ratio' in group


def test_export_periodic_mesh_metadata(run_in_tmpdir):
    cylindrical = openmc.CylindricalMesh(
        r_grid=[0.0, 1.0, 2.0], phi_grid=[0.0, np.pi],
        z_grid=[-1.0, 1.0], origin=(1.0, 2.0, 3.0), mesh_id=10,
        name='cylindrical')
    spherical = openmc.SphericalMesh(
        r_grid=[0.0, 2.0], theta_grid=[0.0, np.pi],
        phi_grid=[0.0, 2.0 * np.pi], origin=(-1.0, -2.0, -3.0), mesh_id=11,
        name='spherical')
    windows = openmc.WeightWindowsList([
        openmc.WeightWindows(cylindrical, [1.0, 2.0],
                             upper_bound_ratio=5.0, id=10),
        openmc.WeightWindows(spherical, [3.0],
                             upper_bound_ratio=5.0, id=11),
    ])

    windows.export_to_hdf5('ww.h5')
    roundtrip = openmc.WeightWindowsList.from_hdf5('ww.h5')
    meshes = {ww.mesh.id: ww.mesh for ww in roundtrip}

    assert meshes[10].name == cylindrical.name
    assert meshes[11].name == spherical.name
    assert np.allclose(meshes[10].origin, cylindrical.origin)
    assert np.allclose(meshes[11].origin, spherical.origin)
    for ww in roundtrip:
        assert ww.energy_bounds[0] == 0.0
        assert ww.energy_bounds[-1] == np.finfo(float).max


def test_mesh_to_lib_object(run_in_tmpdir):
    mesh = openmc.RegularMesh(mesh_id=17, name='runtime mesh')
    mesh.dimension = (2, 3)
    mesh.lower_left = (0.0, 1.0)
    mesh.upper_right = (2.0, 4.0)

    with pytest.raises(RuntimeError, match='must be initialized'):
        mesh.to_lib_object()

    model = openmc.Model()
    sphere = openmc.Sphere(boundary_type='vacuum')
    model.geometry = openmc.Geometry([openmc.Cell(region=-sphere)])
    model.settings.particles = 1
    model.settings.batches = 1

    with openmc.lib.TemporarySession(model):
        lib_mesh = mesh.to_lib_object()
        assert isinstance(lib_mesh, openmc.lib.RegularMesh)
        assert lib_mesh.id == mesh.id
        assert lib_mesh.name == mesh.name
        assert tuple(lib_mesh.dimension) == mesh.dimension
        assert np.allclose(lib_mesh.lower_left, mesh.lower_left)
        assert np.allclose(lib_mesh.upper_right, mesh.upper_right)


@pytest.mark.parametrize('library', ('libmesh', 'moab'))
def test_export_hdf5_unstructured_mesh(request, run_in_tmpdir, library):
    if library == 'libmesh' and not openmc.lib._libmesh_enabled():
        pytest.skip('LibMesh not enabled in this build.')
    if library == 'moab' and not openmc.lib._dagmc_enabled():
        pytest.skip('DAGMC (and MOAB) not enabled in this build.')

    mesh = openmc.UnstructuredMesh(
        request.path.with_name('test_mesh_tets.exo'), library, mesh_id=20,
        name='unstructured', length_multiplier=2.0)
    ww = openmc.WeightWindows(mesh, np.ones(12_000), upper_bound_ratio=5.0)
    openmc.WeightWindowsList([ww]).export_to_hdf5('ww.h5')

    with h5py.File('ww.h5') as f:
        mesh_group = f['meshes'][f'mesh {mesh.id}']
        assert mesh_group['type'][()] == b'unstructured'
        assert mesh_group['name'][()] == b'unstructured'
        assert mesh_group['length_multiplier'][()] == 2.0
