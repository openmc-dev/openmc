.. _io_source:

==================
Source File Format
==================

The current version of the source file format is 18.3. Source files share their
version number with the :ref:`state point file format <io_statepoint>`.

Normally, source data is stored in a state point file. However, it is possible
to request that the source be written separately, in which case the format used
is that documented here.

When surface source writing is triggered, a source file named
``surface_source.h5`` is written with only the sources on specified surfaces,
following the same format.

**/**

:Attributes: - **filetype** (*char[]*) -- String indicating the type of file.
             - **version** (*int[2]*) -- Major and minor version of the source
               file format.
             - **birth_mesh_id** (*int*) -- ID of the mesh used to tag fission
               sites with the ``<birth_mesh>`` settings element. Only present
               together with the ``birth_mesh_bin`` dataset.

:Datasets:

           - **source_bank** (Compound type) -- Source bank information for each
             particle. The compound type has fields ``r``, ``u``, ``E``,
             ``time``, ``wgt``, ``delayed_group``, ``surf_id`` and ``particle``,
             which represent the position, direction, energy, time, weight,
             delayed group, surface ID, and particle type (PDG number),
             respectively.
           - **birth_mesh_bin** (*int[]*) -- Fission bank birth tags, aligned
             by index with ``source_bank``, as described for the
             :ref:`state point file <io_statepoint>`. Only present in source
             point files of an eigenvalue calculation for which the
             ``<birth_mesh>`` settings element was specified; surface source
             files and the initial source file never contain this dataset.
             MCPL source files do not carry the tags. When a file containing
             this dataset is read, the dataset must have the same length as
             ``source_bank``.
