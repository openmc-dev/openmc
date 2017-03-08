.. _io_particle_restart:

============================
Particle Restart File Format
============================

The current version of the particle restart file format is 2.0.

**/**

:Attributes: - **filetype** (*char[]*) -- String indicating the type of file.
             - **version** (*int[2]*) -- Major and minor version of the particle
               restart file format.
             - **openmc_version** (*int[3]*) -- Major, minor, and release
               version number for OpenMC.
             - **git_sha1** (*char[40]*) -- Git commit SHA-1 hash.

:Datasets: - **current_batch** (*int*) -- The number of batches already
             simulated.
           - **generations_per_batch** (*int*) -- Number of generations per
             batch.
           - **current_generation** (*int*) -- The number of generations already
             simulated.
           - **n_particles** (*int8_t*) -- Number of particles used per
             generation.
           - **run_mode** (*char[]*) -- Run mode used, either 'fixed source',
             'eigenvalue', or 'particle restart'.
           - **id** (*int8_t*) -- Unique identifier of the particle.
           - **weight** (*double*) -- Weight of the particle.
           - **energy** (*double*) -- Energy of the particle in eV for
             continuous-energy mode, or the energy group of the particle for
             multi-group mode.
           - **xyz** (*double[3]*) -- Position of the particle.
           - **uvw** (*double[3]*) -- Direction of the particle.
