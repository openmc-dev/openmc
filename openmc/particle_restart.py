import h5py

import openmc.checkvalue as cv

_VERSION_PARTICLE_RESTART = 2

class Particle(object):
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
        self._f = h5py.File(filename, 'r')

        # Ensure filetype and version are correct
        cv.check_filetype_version(self._f, 'particle restart',
                                  _VERSION_PARTICLE_RESTART)

    @property
    def current_batch(self):
        return self._f['current_batch'].value

    @property
    def current_generation(self):
        return self._f['current_generation'].value

    @property
    def energy(self):
        return self._f['energy'].value

    @property
    def generations_per_batch(self):
        return self._f['generations_per_batch'].value

    @property
    def id(self):
        return self._f['id'].value

    @property
    def n_particles(self):
        return self._f['n_particles'].value

    @property
    def run_mode(self):
        return self._f['run_mode'].value.decode()

    @property
    def uvw(self):
        return self._f['uvw'].value

    @property
    def weight(self):
        return self._f['weight'].value

    @property
    def xyz(self):
        return self._f['xyz'].value
