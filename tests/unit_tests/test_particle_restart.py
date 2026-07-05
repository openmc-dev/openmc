import h5py
import numpy as np
import pytest

import openmc


def write_particle_restart_file(filename, *, current_batch=10,
                                current_generation=1, generations_per_batch=1,
                                n_particles=1000, run_mode='fixed source',
                                particle_id=542, particle_type=2112,
                                weight=1.0, energy=1.0e6,
                                xyz=(0., 0., 0.), uvw=(1., 0., 0.),
                                time=0.0):
    """Create a file matching Particle::write_restart in src/particle.cpp."""
    with h5py.File(filename, 'w') as f:
        f.attrs['filetype'] = np.bytes_('particle restart')
        f.attrs['version'] = np.array([2, 1])
        f.create_dataset('current_batch', data=current_batch)
        f.create_dataset('generations_per_batch', data=generations_per_batch)
        f.create_dataset('current_generation', data=current_generation)
        f.create_dataset('n_particles', data=n_particles)
        f.create_dataset('run_mode', data=np.bytes_(run_mode))
        f.create_dataset('id', data=particle_id)
        f.create_dataset('type', data=particle_type)
        f.create_dataset('weight', data=weight)
        f.create_dataset('energy', data=energy)
        f.create_dataset('xyz', data=np.asarray(xyz))
        f.create_dataset('uvw', data=np.asarray(uvw))
        f.create_dataset('time', data=time)


@pytest.fixture
def particle_files(tmp_path):
    params = [
        {'current_batch': 1, 'current_generation': 1, 'particle_id': 285880,
         'particle_type': 2112, 'weight': 1.0, 'energy': 2.0e6,
         'xyz': (1., 2., 3.), 'uvw': (1., 0., 0.)},
        {'current_batch': 1, 'current_generation': 2, 'particle_id': 300145,
         'particle_type': 22, 'weight': 0.5, 'energy': 5.0e5,
         'xyz': (-4., 5., -6.), 'uvw': (0., 1., 0.)},
        {'current_batch': 3, 'current_generation': 1, 'particle_id': 12,
         'particle_type': 2112, 'weight': 2.0, 'energy': 14.1e6,
         'xyz': (7., -8., 9.), 'uvw': (0., 0., -1.)},
    ]
    files = []
    for i, kwargs in enumerate(params):
        path = tmp_path / f'particle_{kwargs["current_batch"]}_{i}.h5'
        write_particle_restart_file(path, **kwargs)
        files.append(path)
    return files, params


def test_read_particle_restart(particle_files):
    files, params = particle_files
    p = openmc.Particle(files[0])
    assert p.current_batch == params[0]['current_batch']
    assert p.current_generation == params[0]['current_generation']
    assert p.id == params[0]['particle_id']
    assert p.type == openmc.ParticleType.NEUTRON
    assert p.weight == params[0]['weight']
    assert p.energy == params[0]['energy']
    assert p.run_mode == 'fixed source'
    np.testing.assert_array_equal(p.xyz, params[0]['xyz'])
    np.testing.assert_array_equal(p.uvw, params[0]['uvw'])


def test_write_vtkhdf(particle_files, tmp_path):
    files, params = particle_files
    out = tmp_path / 'lost_particles.vtkhdf'
    openmc.Particle.write_vtkhdf(files, out)

    with h5py.File(out, 'r') as f:
        root = f['VTKHDF']
        assert root.attrs['Type'] == b'PolyData'
        assert root.attrs['Version'][0] >= 2

        n = len(files)
        assert root['NumberOfPoints'][()] == [n]
        assert root['Points'].shape == (n, 3)
        np.testing.assert_array_equal(
            root['Points'][()], [p['xyz'] for p in params])

        # Each point is a vertex cell
        vertices = root['Vertices']
        assert vertices['NumberOfCells'][()] == [n]
        assert vertices['NumberOfConnectivityIds'][()] == [n]
        np.testing.assert_array_equal(
            vertices['Connectivity'][()], np.arange(n))
        np.testing.assert_array_equal(
            vertices['Offsets'][()], np.arange(n + 1))

        # Other PolyData topologies present but empty
        for name in ('Lines', 'Polygons', 'Strips'):
            group = root[name]
            assert group['NumberOfCells'][()] == [0]
            assert group['NumberOfConnectivityIds'][()] == [0]
            assert group['Connectivity'].shape == (0,)
            np.testing.assert_array_equal(group['Offsets'][()], [0])

        point_data = root['PointData']
        np.testing.assert_array_equal(
            point_data['energy'][()], [p['energy'] for p in params])
        np.testing.assert_array_equal(
            point_data['weight'][()], [p['weight'] for p in params])
        np.testing.assert_array_equal(
            point_data['direction'][()], [p['uvw'] for p in params])
        np.testing.assert_array_equal(
            point_data['particle_type'][()],
            [p['particle_type'] for p in params])
        np.testing.assert_array_equal(
            point_data['batch'][()], [p['current_batch'] for p in params])
        np.testing.assert_array_equal(
            point_data['generation'][()],
            [p['current_generation'] for p in params])
        np.testing.assert_array_equal(
            point_data['id'][()], [p['particle_id'] for p in params])


def test_write_vtkhdf_single_file(particle_files, tmp_path):
    files, params = particle_files
    out = tmp_path / 'single.vtkhdf'
    openmc.Particle.write_vtkhdf(str(files[0]), out)

    with h5py.File(out, 'r') as f:
        root = f['VTKHDF']
        assert root['NumberOfPoints'][()] == [1]
        np.testing.assert_array_equal(root['Points'][()], [params[0]['xyz']])


def test_write_vtkhdf_no_files(tmp_path):
    with pytest.raises(ValueError):
        openmc.Particle.write_vtkhdf([], tmp_path / 'empty.vtkhdf')


def test_write_vtkhdf_read_with_vtk(particle_files, tmp_path):
    vtk = pytest.importorskip('vtk')
    files, params = particle_files
    out = tmp_path / 'lost_particles.vtkhdf'
    openmc.Particle.write_vtkhdf(files, out)

    reader = vtk.vtkHDFReader()
    reader.SetFileName(str(out))
    reader.Update()
    data = reader.GetOutput()

    assert data.GetNumberOfPoints() == len(files)
    assert data.GetNumberOfVerts() == len(files)
    point_data = data.GetPointData()
    names = {point_data.GetArrayName(i)
             for i in range(point_data.GetNumberOfArrays())}
    assert {'energy', 'weight', 'direction', 'particle_type', 'batch',
            'generation', 'id'} <= names
    energy = point_data.GetArray('energy')
    values = [energy.GetValue(i) for i in range(energy.GetNumberOfValues())]
    assert values == [p['energy'] for p in params]
