from collections import namedtuple
from collections.abc import Sequence

import h5py

from .checkvalue import check_filetype_version
from .source import SourceParticle, ParticleType

from pathlib import Path

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
    (weight), ``cell_id`` (cell ID) , ``cell_instance`` (cell instance), and
    ``material_id`` (material ID).

"""
def _particle_track_repr(self):
    return f"<ParticleTrack: {self.particle}, {len(self.states)} states>"
ParticleTrack.__repr__ = _particle_track_repr


_VERSION_TRACK = 3


def _identifier(dset_name):
    """Return (batch, gen, particle) tuple given dataset name"""
    _, batch, gen, particle = dset_name.split('_')
    return (int(batch), int(gen), int(particle))


class Track(Sequence):
    """Tracks resulting from a single source particle

    This class stores information for all tracks resulting from a primary source
    particle and any secondary particles that it created. The track for each
    primary/secondary particle is stored in the :attr:`particle_tracks`
    attribute.

    .. versionadded:: 0.13.1

    Parameters
    ----------
    dset : h5py.Dataset
        Dataset to read track data from

    Attributes
    ----------
    identifier : tuple
        Tuple of (batch, generation, particle number)
    particle_tracks : list
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
        self.particle_tracks = tracks_list

    def __repr__(self):
        return f'<Track {self.identifier}: {len(self.particle_tracks)} particles>'

    def __getitem__(self, index):
        return self.particle_tracks[index]

    def __len__(self):
        return len(self.particle_tracks)

    def filter(self, particle=None, state_filter=None):
        """Filter particle tracks by given criteria

        Parameters
        ----------
        particle : {'neutron', 'photon', 'electron', 'positron'}
            Matching particle type
        state_filter : function
            Function that takes a state (structured datatype) and returns a bool
            depending on some criteria.

        Returns
        -------
        Track
            New instance with only matching :class:`openmc.ParticleTrack` objects

        Examples
        --------
        Get all particle tracks for photons:

        >>> track.filter(particle='photon')

        Get all particle tracks that entered cell with ID=15:

        >>> track.filter(state_filter=lambda s: s['cell_id'] == 15)

        Get all particle tracks in entered material with ID=2:

        >>> track.filter(state_filter=lambda s: s['material_id'] == 2)

        See Also
        --------
        openmc.ParticleTrack

        """
        matching = []
        for t in self:
            # Check for matching particle
            if particle is not None:
                if t.particle.name.lower() != particle:
                    continue

            # Apply arbitrary state filter
            match = True
            if state_filter is not None:
                for state in t.states:
                    if state_filter(state):
                        break
                else:
                    match = False

            if match:
                matching.append(t)

        # Return new Track instance with only matching particle tracks
        track = type(self).__new__(type(self))
        track.identifier = self.identifier
        track.particle_tracks = matching
        return track

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
        for _, states in self:
            r = states['r']
            ax.plot3D(r['x'], r['y'], r['z'])

        return ax

    @property
    def sources(self):
        sources = []
        for particle_track in self:
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


class Tracks(list):
    """Collection of particle tracks

    This class behaves like a list and can be indexed using the normal subscript
    notation. Each element in the list is a :class:`openmc.Track` object.

    .. versionadded:: 0.13.1

    Parameters
    ----------
    filepath : str or pathlib.Path
        Path of file to load

    """

    def __init__(self, filepath='tracks.h5'):
        # Read data from track file
        with h5py.File(filepath, 'r') as fh:
            # Check filetype and version
            check_filetype_version(fh, 'track', _VERSION_TRACK)

            for dset_name in sorted(fh, key=_identifier):
                dset = fh[dset_name]
                self.append(Track(dset))

    def filter(self, particle=None, state_filter=None):
        """Filter tracks by given criteria

        Parameters
        ----------
        particle : {'neutron', 'photon', 'electron', 'positron'}
            Matching particle type
        state_filter : function
            Function that takes a state (structured datatype) and returns a bool
            depending on some criteria.

        Returns
        -------
        Tracks
            List of :class:`openmc.Track` objects

        See Also
        --------
        openmc.Track.filter

        """
        # Create a new Tracks instance but avoid call to __init__
        matching = type(self).__new__(type(self))

        # Append matching Track objects
        for track in self:
            if track.filter(particle, state_filter):
                matching.append(track)
        return matching

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

    def write_to_vtk(self, filename=Path('tracks.vtp')):
        """Creates a VTP file of the tracks

        Parameters
        ----------
        filename : path-like
            Name of the VTP file to write.

        Returns
        -------
        vtk.vtkPolyData
            the VTK vtkPolyData object produced
        """

        import vtk

        # Initialize data arrays and offset.
        points = vtk.vtkPoints()
        cells = vtk.vtkCellArray()

        point_offset = 0
        for particle in self:
            for pt in particle.particle_tracks:
                for state in pt.states:
                    points.InsertNextPoint(state['r'])

                # Create VTK line and assign points to line.
                n = pt.states.size
                line = vtk.vtkPolyLine()
                line.GetPointIds().SetNumberOfIds(n)
                for i in range(n):
                    line.GetPointIds().SetId(i, point_offset + i)
                point_offset += n

                # Add line to cell array
                cells.InsertNextCell(line)

        data = vtk.vtkPolyData()
        data.SetPoints(points)
        data.SetLines(cells)

        writer = vtk.vtkXMLPPolyDataWriter()
        if vtk.vtkVersion.GetVTKMajorVersion() > 5:
            writer.SetInputData(data)
        else:
            writer.SetInput(data)
        writer.SetFileName(str(filename))  # SetFileName requires a string
        writer.Write()

        return data

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
