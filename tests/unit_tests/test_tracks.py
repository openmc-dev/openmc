from pathlib import Path

import numpy as np
import openmc
import pytest

from tests.testing_harness import config


@pytest.fixture
def sphere_model():
    mat = openmc.Material()
    mat.add_nuclide('Zr90', 1.0)
    mat.set_density('g/cm3', 1.0)

    model = openmc.Model()
    sph = openmc.Sphere(r=25.0, boundary_type='vacuum')
    cell = openmc.Cell(fill=mat, region=-sph)
    model.geometry = openmc.Geometry([cell])

    model.settings.run_mode = 'fixed source'
    model.settings.batches = 1
    model.settings.particles = 100
    model.settings.track = [(1, 1, 1), (1, 1, 10), (1, 1, 75)]

    return model


def test_tracks(sphere_model, run_in_tmpdir):
    # If running in MPI mode, setup proper keyword arguments for run()
    kwargs = {'openmc_exec': config['exe']}
    if config['mpi']:
        kwargs['mpi_args'] = [config['mpiexec'], '-n', config['mpi_np']]
    sphere_model.run(**kwargs)

    if config['mpi'] and int(config['mpi_np']) > 1:
        # With MPI, we need to combine track files
        track_files = Path.cwd().glob('tracks_p*.h5')
        openmc.TrackFile.combine(track_files, 'tracks.h5')
    else:
        track_file = Path('tracks.h5')
        assert track_file.is_file()

    # Open track file and make sure we have correct number of tracks
    tracks = openmc.TrackFile('tracks.h5')
    assert len(tracks) == len(sphere_model.settings.track)

    for track, identifier in zip(tracks, sphere_model.settings.track):
        # Check attributes on Track object
        assert isinstance(track, openmc.Track)
        assert track.identifier == identifier
        assert isinstance(track.particles, list)
        assert len(track.particles) == 1

        # Check attributes on ParticleTrack object
        particle_track = track.particles[0]
        assert isinstance(particle_track, openmc.ParticleTrack)
        assert particle_track.particle == openmc.ParticleType.NEUTRON
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
        assert len(sources) == 1
        x = sources[0]
        state = particle_track.states[0]
        assert x.r == (*state['r'],)
        assert x.u == (*state['u'],)
        assert x.E == state['E']
        assert x.time == state['time']
        assert x.wgt == state['wgt']
        assert x.particle == particle_track.particle
