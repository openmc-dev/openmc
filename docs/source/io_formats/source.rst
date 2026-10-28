.. _io_source:

==================
Source File Format
==================

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
             - **energy_mode** (*char[]*) -- Whether ``E`` holds energy group
               indices ('multi-group') or energies in [eV]
               ('continuous-energy'). A file without this attribute is read as
               'continuous-energy'. When the file is read by
               :class:`openmc.FileSource`, any other value is an error, and so
               is a 'multi-group' file in a continuous-energy simulation.
               OpenMC writes this attribute, with the value 'multi-group', only
               to source point files written by a multigroup simulation.
               Surface source files and files written by a continuous-energy
               simulation, by :func:`openmc.write_source_file` or by
               :meth:`openmc.ParticleList.export_to_hdf5` do not have this
               attribute, so their ``E`` is read as an energy in [eV].

:Datasets:

           - **source_bank** (Compound type) -- Source bank information for each
             particle. The compound type has fields ``r``, ``u``, ``E``,
             ``time``, ``wgt``, ``delayed_group``, ``surf_id`` and ``particle``,
             which represent the position, direction, energy, time, weight,
             delayed group, surface ID, and particle type (PDG number),
             respectively. When the **energy_mode** attribute is
             'multi-group', ``E`` is the energy group index of the particle
             (0 for the highest-energy group) rather than its energy in [eV].
