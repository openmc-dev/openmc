import struct


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
    gen_per_batch : int
        Number of generations per batch
    current_gen : int
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
        Energy of the particle in MeV
    xyz : list of float
        Position of the particle
    uvw : list of float
        Directional cosines of the particle

    """

    def __init__(self, filename):
        import h5py
        self._f = h5py.File(filename, 'r')

        # Ensure filetype and revision are correct
        if 'filetype' not in self._f or self._f[
                'filetype'].value.decode() != 'particle restart':
            raise IOError('{} is not a particle restart file.'.format(filename))
        if self._f['revision'].value != 1:
            raise IOError('Particle restart file has a file revision of {} '
                          'which is not consistent with the revision this '
                          'version of OpenMC expects ({}).'.format(
                              self._f['revision'].value, 1))

    @property
    def current_batch(self):
        return self._f['current_batch'].value

    @property
    def current_gen(self):
        return self._f['current_gen'].value

    @property
    def energy(self):
        return self._f['energy'].value

    @property
    def gen_per_batch(self):
        return self._f['gen_per_batch'].value

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
