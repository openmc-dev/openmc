from collections import namedtuple

import h5py

from .checkvalue import check_filetype_version
from .source import SourceParticle, ParticleType


ParticleTrack = namedtuple('ParticleTrack', ['particle', 'states'])
ParticleTrack.__doc__ = """\
Particle track information

Parameters
----------
particle : openmc.ParticleType
    Type of the particle
states : numpy.ndarray
    Structured array containing each state of the particle. The structured array
    contains the following fields: ``r`` (position; each direction in [cm]),
    ``u`` (direction), ``E`` (energy in [eV]), ``time`` (time in [s]), ``wgt``
    (weight), ``cell_id`` (cell ID) , and ``material_id`` (material ID).

"""

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
        axes : matplotlib.axes.Axes, optional
            Axes for plot

        Returns
        -------
        axes : matplotlib.axes.Axes
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

        return ax

    @property
    def sources(self):
        sources = []
        for particle_track in self.particles:
            particle_type = ParticleType(particle_track.particle)
            state = particle_track.states[0]
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

    This class behaves like a list and can be indexed using the normal subscript
    notation. Each element in the list is a :class:`openmc.Track` object.

    Parameters
    ----------
    filepath : str or pathlib.Path
        Path of file to load

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
        """Produce a 3D plot of particle tracks

        Returns
        -------
        matplotlib.axes.Axes
            Axes for plot

        """
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.set_xlabel('x [cm]')
        ax.set_ylabel('y [cm]')
        ax.set_zlabel('z [cm]')
        for track in self:
            track.plot(ax)
        return ax

    @staticmethod
    def combine(track_files, path='tracks.h5'):
        """Combine multiple track files into a single track file

        Parameters
        ----------
        track_files : list of path-like
            Paths to track files to combine
        path : path-like
            Path of combined track file to create

        """
        with h5py.File(path, 'w') as h5_out:
            for i, fname in enumerate(track_files):
                with h5py.File(fname, 'r') as h5_in:
                    # Copy file attributes for first file
                    if i == 0:
                        h5_out.attrs['filetype'] = h5_in.attrs['filetype']
                        h5_out.attrs['version'] = h5_in.attrs['version']

                    # Copy each 'track_*' dataset from input file
                    for dset in h5_in:
                        h5_in.copy(dset, h5_out)
