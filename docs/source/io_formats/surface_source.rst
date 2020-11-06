.. _io_surface_source:

==========================
Surface Source File Format
==========================

When surface source writing is triggered, a separate source file is written with
only the sources on specified surfaces. The format is identical to source file.

**/**

:Attributes: - **filetype** (*char[]*) -- String indicating the type of file.

:Datasets:

           - **source_bank** (Compound type) -- Source bank information for each
             particle. The compound type has fields ``r``, ``u``, ``E``,
             ``wgt``, ``delayed_group``, ``surf_id`` and ``particle``,
             which represent the position, direction, energy, weight,
             delayed group, surface ID, and particle type (0=neutron, 1=photon,
             2=electron, 3=positron), respectively.
