from pathlib import Path

import h5py
import numpy as np
import openmc
import pytest

from tests.testing_harness import config


@pytest.fixture
def sphere_model():
    openmc.reset_auto_ids()
    mat = openmc.Material()
    mat.add_nuclide('Zr90', 1.0)
    mat.set_density('g/cm3', 1.0)

    model = openmc.Model()
    sph = openmc.Sphere(r=25.0, boundary_type='vacuum')
    cell = openmc.Cell(fill=mat, region=-sph)
    model.geometry = openmc.Geometry([cell])

    model.settings.run_mode = 'fixed source'
    model.settings.batches = 2
    model.settings.particles = 50

    return model


def generate_track_file(model, **kwargs):
    # If running in MPI mode, setup proper keyword arguments for run()
    kwargs.setdefault('openmc_exec', config['exe'])
    if config['mpi']:
        kwargs['mpi_args'] = [config['mpiexec'], '-n', config['mpi_np']]
    model.run(**kwargs)

    if config['mpi'] and int(config['mpi_np']) > 1:
        # With MPI, we need to combine track files
        track_files = Path.cwd().glob('tracks_p*.h5')
        openmc.Tracks.combine(track_files, 'tracks.h5')
    else:
        track_file = Path('tracks.h5')
        assert track_file.is_file()


@pytest.mark.parametrize("particle", ["neutron", "photon"])
def test_tracks(sphere_model, particle, run_in_tmpdir):
    # Set track identifiers
    sphere_model.settings.track = [(1, 1, 1), (1, 1, 10), (2, 1, 15)]

    # Set source particle
    sphere_model.settings.source = openmc.IndependentSource(particle=particle)

    # Run OpenMC to generate tracks.h5 file
    generate_track_file(sphere_model)

    # Open track file and make sure we have correct number of tracks
    tracks = openmc.Tracks('tracks.h5')
    assert len(tracks) == len(sphere_model.settings.track)

    for track, identifier in zip(tracks, sphere_model.settings.track):
        # Check attributes on Track object
        assert isinstance(track, openmc.Track)
        assert track.identifier == identifier
        assert isinstance(track.particle_tracks, list)
        if particle == 'neutron':
            assert len(track.particle_tracks) == 1

        # Check attributes on ParticleTrack object
        particle_track = track.particle_tracks[0]
        assert isinstance(particle_track, openmc.ParticleTrack)

        assert particle_track.particle.name.lower() == particle
        assert isinstance(particle_track.states, np.ndarray)

        # Sanity checks on actual data
        for state in particle_track.states:
            assert np.linalg.norm([*state['r']]) <= 25.0001
            assert np.linalg.norm([*state['u']]) == pytest.approx(1.0)
            assert 0.0 <= state['E'] <= 20.0e6
            assert state['time'] >= 0.0
            assert 0.0 <= state['wgt'] <= 1.0
            assert state['cell_id'] == 1
            assert state['material_id'] == 1

        # Checks on 'sources' property
        sources = track.sources
        assert len(sources) == len(track.particle_tracks)
        x = sources[0]
        state = particle_track.states[0]
        assert x.r == (*state['r'],)
        assert x.u == (*state['u'],)
        assert x.E == state['E']
        assert x.time == state['time']
        assert x.wgt == state['wgt']
        assert x.particle == particle_track.particle


def test_max_tracks(sphere_model, run_in_tmpdir):
    # Set maximum number of tracks per process to write
    sphere_model.settings.max_tracks = expected_num_tracks = 10
    if config['mpi']:
        expected_num_tracks *= int(config['mpi_np'])

    # Run OpenMC to generate tracks.h5 file
    generate_track_file(sphere_model, tracks=True)

    # Open track file and make sure we have correct number of tracks
    tracks = openmc.Tracks('tracks.h5')
    assert len(tracks) == expected_num_tracks


def test_filter(sphere_model, run_in_tmpdir):
    # Set maximum number of tracks per process to write
    sphere_model.settings.max_tracks = 25
    sphere_model.settings.photon_transport = True

    # Run OpenMC to generate tracks.h5 file
    generate_track_file(sphere_model, tracks=True)

    tracks = openmc.Tracks('tracks.h5')
    for track in tracks:
        # Test filtering by particle
        matches = track.filter(particle='photon')
        for x in matches:
            assert x.particle == openmc.ParticleType.PHOTON

        # Test general state filter
        matches = track.filter(state_filter=lambda s: s['cell_id'] == 1)
        assert isinstance(matches, openmc.Track)
        assert matches.particle_tracks == track.particle_tracks
        matches = track.filter(state_filter=lambda s: s['cell_id'] == 2)
        assert matches.particle_tracks == []
        matches = track.filter(state_filter=lambda s: s['E'] < 0.0)
        assert matches.particle_tracks == []

    # Test filter method on Tracks
    matches = tracks.filter(particle='neutron')
    assert isinstance(matches, openmc.Tracks)
    assert matches == tracks
    matches = tracks.filter(state_filter=lambda s: s['E'] > 0.0)
    assert matches == tracks
    matches = tracks.filter(particle='bunnytron')
    assert matches == []


def test_write_to_vtk(sphere_model):
    vtk = pytest.importorskip('vtk')
    # Set maximum number of tracks per process to write
    sphere_model.settings.max_tracks = 25
    sphere_model.settings.photon_transport = True

    # Run OpenMC to generate tracks.h5 file
    generate_track_file(sphere_model, tracks=True)

    tracks = openmc.Tracks('tracks.h5')
    polydata = tracks.write_to_vtk('tracks.vtp')

    assert isinstance(polydata, vtk.vtkPolyData)
    assert Path('tracks.vtp').is_file()


def test_restart_track(run_in_tmpdir, sphere_model):
    # cut the sphere model in half with an improper boundary condition
    plane = openmc.XPlane(x0=-1.0)
    for cell in sphere_model.geometry.get_all_cells().values():
        cell.region &= +plane

    # generate lost particle files
    with pytest.raises(RuntimeError, match='Maximum number of lost particles has been reached.'):
        sphere_model.run(output=False, threads=1)

    lost_particle_files = list(Path.cwd().glob('particle_*.h5'))
    assert len(lost_particle_files) > 0
    particle_file = lost_particle_files[0]
     # restart the lost particle with tracks enabled
    sphere_model.run(tracks=True, restart_file=particle_file)
    tracks_file = Path('tracks.h5')
    assert tracks_file.is_file()

    # check that the last track of the file matches the lost particle file
    tracks = openmc.Tracks(tracks_file)
    initial_state = tracks[0].particle_tracks[0].states[0]
    restart_r = np.array(initial_state['r'])
    restart_u = np.array(initial_state['u'])

    with h5py.File(particle_file, 'r') as lost_particle_file:
        lost_r = np.array(lost_particle_file['xyz'][()])
        lost_u = np.array(lost_particle_file['uvw'][()])

    pytest.approx(restart_r, lost_r)
    pytest.approx(restart_u, lost_u)
