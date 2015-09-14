.. _usersguide_source:

==================
Source File Format
==================

Normally, source data is stored in a state point file. However, it is possible
to request that the source be written separately, in which case the format used
is that documented here.

**/filetype** (*int*)

    Flags what type of file this is. A value of -1 indicates a statepoint file,
    a value of -2 indicates a particle restart file, a value of -3 indicates a
    source file, and a value of -4 indicates a track file.

**/source_bank** (Compound type)

    Source bank information for each particle. The compound type has fields
    ``wgt``, ``xyz``, ``uvw``, and ``E`` which represent the weight, position,
    direction, and energy of the source particle, respectively.
