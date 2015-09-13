.. _usersguide_track:

=================
Track File Format
=================

The current revision of the particle track file format is 1.

**/filetype** (*int*)

    Flags what type of file this is. A value of -1 indicates a statepoint file,
    a value of -2 indicates a particle restart file, a value of -3 indicates a
    source file, and a value of -4 indicates a track file.

**/revision** (*int*)

    Revision of the track file format. Any time a change is made in the format,
    this integer is incremented.

**/n_particles** (*int*)

    Number of particles for which tracks are recorded.

**/n_coords** (*int[]*)

    Number of coordinates for each particle.

*do i = 1, n_particles*

    **/coordinates_i** (*double[][3]*)

        (x,y,z) coordinates for the *i*-th particle.
