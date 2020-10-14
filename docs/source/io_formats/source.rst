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
             particle. The compound type has fields ``r``, ``u``, ``E``,
             ``wgt``, ``delayed_group``, and ``particle``, which represent the
             position, direction, energy, weight, delayed group, and particle
             type (0=neutron, 1=photon, 2=electron, 3=positron), respectively.
