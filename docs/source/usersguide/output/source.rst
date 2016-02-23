.. _usersguide_source:

==================
Source File Format
==================

Normally, source data is stored in a state point file. However, it is possible
to request that the source be written separately, in which case the format used
is that documented here.

**/filetype** (*char[]*)

    String indicating the type of file.

**/source_bank** (Compound type)

    Source bank information for each particle. The compound type has fields
    ``wgt``, ``xyz``, ``uvw``, ``E``, and ``delayed_group``, which
    represent the weight, position, direction, energy, energy group, and
    delayed_group of the source particle, respectively.
