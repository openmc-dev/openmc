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
    filetype : int
        Integer indicating the file type
    revision : int
        Revision of the particle restart format
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

        # Read all metadata
        self._read_data()

    def _read_data(self):
        # Read filetype
        self.filetype = self._get_int(path='filetype')[0]

        # Read statepoint revision
        self.revision = self._get_int(path='revision')[0]

        # Read current batch
        self.current_batch = self._get_int(path='current_batch')[0]

        # Read run information
        self.gen_per_batch = self._get_int(path='gen_per_batch')[0]
        self.current_gen = self._get_int(path='current_gen')[0]
        self.n_particles = self._get_long(path='n_particles')[0]
        self.run_mode = self._get_int(path='run_mode')[0]

        # Read particle properties
        self.id = self._get_long(path='id')[0]
        self.weight = self._get_double(path='weight')[0]
        self.energy = self._get_double(path='energy')[0]
        self.xyz = self._get_double(3, path='xyz')
        self.uvw = self._get_double(3, path='uvw')

    def _get_int(self, n=1, path=None):
        return [int(v) for v in self._f[path].value]

    def _get_long(self, n=1, path=None):
        return [int(v) for v in self._f[path].value]

    def _get_float(self, n=1, path=None):
        return [float(v) for v in self._f[path].value]

    def _get_double(self, n=1, path=None):
        return [float(v) for v in self._f[path].value]

    def _get_string(self, n=1, path=None):
        return str(self._f[path].value)
