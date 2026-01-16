.. _io_collision_track:

===========================
Collision Track File Format
===========================

When collision tracking is enabled with ``mcpl=false`` (the default), OpenMC
writes binary data to an HDF5 file named ``collision_track.h5``. The same data
may also be written after each batch when multiple files are requested
(``collision_track.N.h5``) or when the run is performed in parallel. The file
contains the information needed to reconstruct each recorded collision.

The current revision of the collision track file format is 1.0.

**/**

:Attributes:
            - **filetype** (*char[]*) -- String indicating the type of file.
              For collision-track files the value is ``"collision_track"``.

:Datasets:

          - **collision_track_bank** (Compound type) -- Collision information
            for each stored event. Each entry in the dataset corresponds to one
            collision and contains the following fields:

              - ``r`` (*double[3]*) -- Position of the collision in [cm].
              - ``u`` (*double[3]*) -- Direction unit vector immediately after the collision.
              - ``E`` (*double*) -- Incident particle energy before the collision in [eV].
              - ``dE`` (*double*) -- Energy loss over the collision (:math:`E_\text{before} - E_\text{after}`) in [eV].
              - ``time`` (*double*) -- Time of the collision in [s].
              - ``wgt`` (*double*) -- Particle weight at the collision.
              - ``event_mt`` (*int*) -- ENDF MT number identifying the reaction.
              - ``delayed_group`` (*int*) -- Delayed neutron group index (non-zero for delayed events).
              - ``cell_id`` (*int*) -- ID of the cell in which the collision occurred.
              - ``nuclide_id`` (*int*) -- ZA identifier of the nuclide (ZZZAAAM format).
              - ``material_id`` (*int*) -- ID of the material containing the collision site.
              - ``universe_id`` (*int*) -- ID of the universe containing the collision site.
              - ``n_collision`` (*int*) -- Collision counter for the particle history.
              - ``particle`` (*int*) -- Particle type (0=neutron, 1=photon, 2=electron, 3=positron).
              - ``parent_id`` (*int64*) -- Unique ID of the parent particle.
              - ``progeny_id`` (*int64*) -- Progeny ID of the particle.

In an MPI run, OpenMC writes the combined dataset by gathering collision-track
entries from all ranks before flushing them to disk, so the final file appears
as though it were produced serially.
