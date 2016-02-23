.. _usersguide_particle_restart:

============================
Particle Restart File Format
============================

The current revision of the particle restart file format is 1.

**/filetype** (*char[]*)

    String indicating the type of file.

**/revision** (*int*)

    Revision of the particle restart file format. Any time a change is made in
    the format, this integer is incremented.

**/current_batch** (*int*)

    The number of batches already simulated.

**/gen_per_batch** (*int*)

    Number of generations per batch.

**/current_gen** (*int*)

    The number of generations already simulated.

**/n_particles** (*int8_t*)

    Number of particles used per generation.

**/run_mode** (*int*)

    Run mode used. A value of 1 indicates a fixed-source run and a value of 2
    indicates an eigenvalue run.

**/id** (*int8_t*)

    Unique identifier of the particle.

**/weight** (*double*)

    Weight of the particle.

**/energy** (*double*)

    Energy of the particle in MeV for continuous-energy mode, or the energy
    group of the particle for multi-group mode.

**/xyz** (*double[3]*)

    Position of the particle.

**/uvw** (*double[3]*)

    Direction of the particle.
