.. _io_source:

==================
Source File Format
==================

Normally, source data is stored in a state point file. However, it is possible
to request that the source be written separately, in which case the format used
is that documented here.

**/**

:Attributes: - **filetype** (*char[]*) -- String indicating the type of file.

:Datasets:

           - **source_bank** (Compound type) -- Source bank information for each
             particle. The compound type has fields ``wgt``, ``xyz``, ``uvw``,
             ``E``, ``delayed_group``, and ``particle``, which represent the
             weight, position, direction, energy, energy group, delayed group,
             and type of the source particle, respectively.
