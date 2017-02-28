.. _io_track:

=================
Track File Format
=================

The current revision of the particle track file format is 2.0.

**/**

:Attributes: - **filetype** (*char[]*) -- String indicating the type of file.
             - **version** (*int[2]*) -- Major and minor version of the track
               file format.
             - **n_particles** (*int*) -- Number of particles for which tracks
               are recorded.
             - **n_coords** (*int[]*) -- Number of coordinates for each
               particle.

:Datasets:
           - **coordinates_<i>** (*double[][3]*) -- (x,y,z) coordinates for the
             *i*-th particle.
