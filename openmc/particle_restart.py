import h5py

import openmc.checkvalue as cv

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
    type : int
        Particle type (1 = neutron, 2 = photon, 3 = electron, 4 = positron)
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
        with h5py.File(filename, "r") as f:

            # Ensure filetype and version are correct
            cv.check_filetype_version(f, "particle restart", _VERSION_PARTICLE_RESTART)

            self.current_batch = f["current_batch"][()]
            self.current_generation = f["current_generation"][()]
            self.energy = f["energy"][()]
            self.generations_per_batch = f["generations_per_batch"][()]
            self.id = f["id"][()]
            self.type = f["type"][()]
            self.n_particles = f["n_particles"][()]
            self.run_mode = f["run_mode"][()].decode()
            self.uvw = f["uvw"][()]
            self.weight = f["weight"][()]
            self.xyz = f["xyz"][()]
