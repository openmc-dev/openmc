.. _usersguide_track:

=================
Track File Format
=================

The current revision of the particle track file format is 1.

**/filetype** (*char[]*)

    String indicating the type of file.

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
