from collections import namedtuple

import h5py

from .checkvalue import check_filetype_version
from .source import SourceParticle, ParticleType


ParticleTrack = namedtuple('ParticleTrack', ['particle', 'states'])

_VERSION_TRACK = 3


def _identifier(dset_name):
    """Return (batch, gen, particle) tuple given dataset name"""
    _, batch, gen, particle = dset_name.split('_')
    return (int(batch), int(gen), int(particle))


class Track:
    """Tracks resulting from a single source particle

    Parameters
    ----------
    dset : h5py.Dataset
        Dataset to read track data from

    Attributes
    ----------
    identifier : tuple
        Tuple of (batch, generation, particle number)
    particles : list
        List of tuples containing (particle type, array of track states)
    sources : list
        List of :class:`SourceParticle` representing each primary/secondary
        particle

    """

    def __init__(self, dset):
        tracks = dset[()]
        offsets = dset.attrs['offsets']
        particles = dset.attrs['particles']
        self.identifier = _identifier(dset.name)

        # Construct list of track histories
        tracks_list = []
        for particle, start, end in zip(particles, offsets[:-1], offsets[1:]):
            ptype = ParticleType(particle)
            tracks_list.append(ParticleTrack(ptype, tracks[start:end]))
        self.particles = tracks_list

    def __repr__(self):
        return f'<Track {self.identifier}: {len(self.particles)} particles>'

    def plot(self, axes=None):
        """Produce a 3D plot of particle tracks

        Parameters
        ----------
        axes : matplotlib.pyplot.Axes, optional
            Axes for plot
        """
        import matplotlib.pyplot as plt

        # Setup axes is one wasn't passed
        if axes is None:
            fig = plt.figure()
            ax = plt.axes(projection='3d')
            ax.set_xlabel('x [cm]')
            ax.set_ylabel('y [cm]')
            ax.set_zlabel('z [cm]')
        else:
            ax = axes

        # Plot each particle track
        for _, states in self.particles:
            r = states['r']
            ax.plot3D(r['x'], r['y'], r['z'])

        if axes is None:
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


class TrackFile(list):
    """Collection of particle tracks

    Parameters
    ----------
    filepath : str or pathlib.Path
        Path of file to load

    Attributes
    ----------
    tracks : list
        List of :class:`Track` objects

    """

    def __init__(self, filepath):
        # Read data from track file
        with h5py.File(filepath, 'r') as fh:
            # Check filetype and version
            check_filetype_version(fh, 'track', _VERSION_TRACK)

            for dset_name in sorted(fh, key=_identifier):
                dset = fh[dset_name]
                self.append(Track(dset))

    def plot(self):
        """Produce a 3D plot of particle tracks"""
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        for track in self:
            track.plot(ax)
        plt.show()
