from collections import namedtuple

import h5py

from .source import SourceParticle, ParticleType


TrackHistory = namedtuple('TrackHistory', ['particle', 'states'])


class TrackFile:
    """Particle track output

    Parameters
    ----------
    filepath : str or pathlib.Path
        Path of file to load

    Attributes
    ----------
    sources : list
        List of :class:`SourceParticle` representing each primary/secondary
        particle
    tracks : list
        List of tuples containing (particle type, array of track states)

    """

    def __init__(self, filepath):
        # Read data from track file
        with h5py.File(filepath, 'r') as fh:
            offsets = fh.attrs['offsets']
            tracks = fh['tracks'][()]
            particles = fh.attrs['particles']

        # Construct list of track histories
        tracks_list = []
        for particle, start, end in zip(particles, offsets[:-1], offsets[1:]):
            ptype = ParticleType(particle)
            tracks_list.append(TrackHistory(ptype, tracks[start:end]))
        self.tracks = tracks_list

    def __repr__(self):
        return f'<TrackFile: {len(self.tracks)} particles>'

    def plot(self):
        """Produce a 3D plot of particle tracks"""
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        for _, states in self.tracks:
            r = states['r']
            ax.plot3D(r['x'], r['y'], r['z'])
        ax.set_xlabel('x [cm]')
        ax.set_ylabel('y [cm]')
        ax.set_zlabel('z [cm]')
        plt.show()

    @property
    def sources(self):
        sources = []
        for track_history in self.tracks:
            particle_type = ParticleType(track_history.particle)
            state = track_history.states[0]
            sources.append(
                SourceParticle(
                    r=state['r'], u=state['u'], E=state['E'],
                    time=state['time'], wgt=state['wgt'],
                    particle=particle_type
                )
            )
        return sources
