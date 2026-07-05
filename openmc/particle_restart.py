from collections.abc import Iterable
import os

import h5py
import numpy as np

import openmc.checkvalue as cv
from .checkvalue import PathLike
from .particle_type import ParticleType

_VERSION_PARTICLE_RESTART = 2


class Particle:
    """Information used to restart a specific particle that caused a simulation to
    fail.

    Parameters
    ----------
    filename : str
        Path to the particle restart file

    Attributes
    ----------
    current_batch : int
        The batch containing the particle
    generations_per_batch : int
        Number of generations per batch
    current_generation : int
        The generation containing the particle
    n_particles : int
        Number of particles per generation
    run_mode : int
        Type of simulation (criticality or fixed source)
    id : long
        Identifier of the particle
    type : openmc.ParticleType
        Particle type
    weight : float
        Weight of the particle
    energy : float
        Energy of the particle in eV
    xyz : list of float
        Position of the particle
    uvw : list of float
        Directional cosines of the particle

    """

    def __init__(self, filename):
        with h5py.File(filename, 'r') as f:

            # Ensure filetype and version are correct
            cv.check_filetype_version(f, 'particle restart', _VERSION_PARTICLE_RESTART)

            self.current_batch = f['current_batch'][()]
            self.current_generation = f['current_generation'][()]
            self.energy = f['energy'][()]
            self.generations_per_batch = f['generations_per_batch'][()]
            self.id = f['id'][()]
            self.type = ParticleType(f['type'][()])
            self.n_particles = f['n_particles'][()]
            self.run_mode = f['run_mode'][()].decode()
            self.uvw = f['uvw'][()]
            self.weight = f['weight'][()]
            self.xyz = f['xyz'][()]

    @staticmethod
    def write_vtkhdf(
        particle_files: PathLike | Iterable[PathLike],
        filename: PathLike = 'lost_particles.vtkhdf'
    ):
        """Write lost particle files to a single VTK HDF PolyData file.

        Each lost particle is written as a vertex with point data arrays for
        energy, weight, direction, particle type, batch, generation, and
        particle ID. The resulting file can be opened directly in ParaView to
        visualize where particles are being lost. Only h5py and numpy are
        required to write the file.

        .. versionadded:: 0.15.4

        Parameters
        ----------
        particle_files : path-like or iterable of path-like
            Path(s) to lost particle restart files
        filename : path-like
            Name of the VTK HDF file to write

        Examples
        --------
        >>> openmc.Particle.write_vtkhdf(
        ...     ['particle_10_542.h5', 'particle_12_178.h5'],
        ...     'lost_particles.vtkhdf'
        ... )

        """
        if isinstance(particle_files, (str, os.PathLike)):
            particle_files = [particle_files]
        particles = [Particle(f) for f in particle_files]
        if not particles:
            raise ValueError('At least one particle file must be provided.')

        n = len(particles)
        points = np.array([p.xyz for p in particles], dtype=np.float64)

        with h5py.File(filename, 'w') as f:
            root = f.create_group('VTKHDF')
            root.attrs['Version'] = (2, 1)
            ascii_type = 'PolyData'.encode('ascii')
            root.attrs.create(
                'Type', ascii_type,
                dtype=h5py.string_dtype('ascii', len(ascii_type)))

            root.create_dataset('NumberOfPoints', data=(n,), dtype='i8')
            root.create_dataset('Points', data=points)

            # Represent each lost particle as a vertex cell
            vertices = root.create_group('Vertices')
            vertices.create_dataset('NumberOfCells', data=(n,), dtype='i8')
            vertices.create_dataset(
                'NumberOfConnectivityIds', data=(n,), dtype='i8')
            vertices.create_dataset(
                'Connectivity', data=np.arange(n, dtype=np.int64))
            vertices.create_dataset(
                'Offsets', data=np.arange(n + 1, dtype=np.int64))

            # The remaining topologies required by VTK HDF PolyData are empty
            for group_name in ('Lines', 'Polygons', 'Strips'):
                group = root.create_group(group_name)
                group.create_dataset('NumberOfCells', data=(0,), dtype='i8')
                group.create_dataset(
                    'NumberOfConnectivityIds', data=(0,), dtype='i8')
                group.create_dataset(
                    'Connectivity', data=np.empty(0, dtype=np.int64))
                group.create_dataset(
                    'Offsets', data=np.zeros(1, dtype=np.int64))

            point_data = root.create_group('PointData')
            point_data.create_dataset('energy', data=np.array(
                [p.energy for p in particles], dtype=np.float64))
            point_data.create_dataset('weight', data=np.array(
                [p.weight for p in particles], dtype=np.float64))
            point_data.create_dataset('direction', data=np.array(
                [p.uvw for p in particles], dtype=np.float64))
            point_data.create_dataset('particle_type', data=np.array(
                [int(p.type) for p in particles], dtype=np.int64))
            point_data.create_dataset('batch', data=np.array(
                [p.current_batch for p in particles], dtype=np.int64))
            point_data.create_dataset('generation', data=np.array(
                [p.current_generation for p in particles], dtype=np.int64))
            point_data.create_dataset('id', data=np.array(
                [p.id for p in particles], dtype=np.int64))
