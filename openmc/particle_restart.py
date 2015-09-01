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
        self.filetype = self._f['filetype'].value

        # Read statepoint revision
        self.revision = self._f['revision'].value

        # Read current batch
        self.current_batch = self._f['current_batch'].value

        # Read run information
        self.gen_per_batch = self._f['gen_per_batch'].value
        self.current_gen = self._f['current_gen'].value
        self.n_particles = self._f['n_particles'].value
        self.run_mode = self._f['run_mode'].value

        # Read particle properties
        self.id = self._f['id'].value
        self.weight = self._f['weight'].value
        self.energy = self._f['energy'].value
        self.xyz = self._f['xyz'].value
        self.uvw = self._f['uvw'].value
