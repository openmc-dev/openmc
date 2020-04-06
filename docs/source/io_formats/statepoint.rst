.. _io_statepoint:

=======================
State Point File Format
=======================

The current version of the statepoint file format is 17.0.

**/**

:Attributes: - **filetype** (*char[]*) -- String indicating the type of file.
             - **version** (*int[2]*) -- Major and minor version of the
               statepoint file format.
             - **openmc_version** (*int[3]*) -- Major, minor, and release
               version number for OpenMC.
             - **git_sha1** (*char[40]*) -- Git commit SHA-1 hash.
             - **date_and_time** (*char[]*) -- Date and time the summary was
               written.
             - **path** (*char[]*) -- Path to directory containing input files.
             - **tallies_present** (*int*) -- Flag indicating whether tallies
               are present (1) or not (0).
             - **source_present** (*int*) -- Flag indicating whether the source
               bank is present (1) or not (0).

:Datasets: - **seed** (*int8_t*) -- Pseudo-random number generator seed.
           - **energy_mode** (*char[]*) -- Energy mode of the run, either
             'continuous-energy' or 'multi-group'.
           - **run_mode** (*char[]*) -- Run mode used, either 'eigenvalue' or
             'fixed source'.
           - **n_particles** (*int8_t*) -- Number of particles used per generation.
           - **n_batches** (*int*) -- Number of batches to simulate.
           - **current_batch** (*int*) -- The number of batches already simulated.
           - **n_inactive** (*int*) -- Number of inactive batches. Only present
             when `run_mode` is 'eigenvalue'.
           - **generations_per_batch** (*int*) -- Number of generations per
             batch. Only present when `run_mode` is 'eigenvalue'.
           - **k_generation** (*double[]*) -- k-effective for each generation
             simulated.
           - **entropy** (*double[]*) -- Shannon entropy for each generation
             simulated.
           - **k_col_abs** (*double*) -- Sum of product of collision/absorption
             estimates of k-effective.
           - **k_col_tra** (*double*) -- Sum of product of
             collision/track-length estimates of k-effective.
           - **k_abs_tra** (*double*) -- Sum of product of
             absorption/track-length estimates of k-effective.
           - **k_combined** (*double[2]*) -- Mean and standard deviation of a
             combined estimate of k-effective.
           - **n_realizations** (*int*) -- Number of realizations for global
             tallies.
           - **global_tallies** (*double[][2]*) -- Accumulated sum and
             sum-of-squares for each global tally.
           - **source_bank** (Compound type) -- Source bank information for each
             particle. The compound type has fields ``wgt``, ``xyz``, ``uvw``,
             ``E``, ``g``, and ``delayed_group``, which represent the weight,
             position, direction, energy, energy group, and delayed_group of the
             source particle, respectively. Only present when `run_mode` is
             'eigenvalue'.

**/tallies/**

:Attributes: - **n_tallies** (*int*) -- Number of user-defined tallies.
             - **ids** (*int[]*) -- User-defined unique ID of each tally.

**/tallies/meshes/**

:Attributes: - **n_meshes** (*int*) -- Number of meshes in the problem.
             - **ids** (*int[]*) -- User-defined unique ID of each mesh.

**/tallies/meshes/mesh <uid>/**

:Datasets: - **type** (*char[]*) -- Type of mesh.
           - **dimension** (*int*) -- Number of mesh cells in each dimension.
           - **lower_left** (*double[]*) -- Coordinates of lower-left corner of
             mesh.
           - **upper_right** (*double[]*) -- Coordinates of upper-right corner
             of mesh.
           - **width** (*double[]*) -- Width of each mesh cell in each
             dimension.
           - **Unstructured Mesh Only:**
              - **volumes** (*double[]*) -- Volume of each mesh cell.
              - **centroids** (*double[]*) -- Location of the mesh cell
                centroids.

**/tallies/filters/**

:Attributes: - **n_filters** (*int*) -- Number of filters in the problem.
             - **ids** (*int[]*) -- User-defined unique ID of each filter.

**/tallies/filters/filter <uid>/**

:Datasets: - **type** (*char[]*) -- Type of the j-th filter. Can be 'universe',
             'material', 'cell', 'cellborn', 'surface', 'mesh', 'energy',
             'energyout', 'distribcell', 'mu', 'polar', 'azimuthal',
             'delayedgroup', or 'energyfunction'.
           - **n_bins** (*int*) -- Number of bins for the j-th filter. Not
             present for 'energyfunction' filters.
           - **bins** (*int[]* or *double[]*) -- Value for each filter bin of
             this type. Not present for 'energyfunction' filters.
           - **energy** (*double[]*) -- Energy grid points for energyfunction
             interpolation. Only used for 'energyfunction' filters.
           - **y** (*double[]*) -- Interpolant values for energyfunction
             interpolation. Only used for 'energyfunction' filters.

**/tallies/derivatives/derivative <id>/**

:Datasets: - **independent variable** (*char[]*) -- Independent variable of
             tally derivative.
           - **material** (*int*) -- ID of the perturbed material.
           - **nuclide** (*char[]*) -- Alias of the perturbed nuclide.
           - **estimator** (*char[]*) -- Type of tally estimator, either
             'analog', 'tracklength', or 'collision'.

**/tallies/tally <uid>/**

:Attributes:
             - **internal** (*int*) -- Flag indicating the presence of tally
               data (0) or absence of tally data (1). All user defined
               tallies will have a value of 0 unless otherwise instructed.

:Datasets: - **n_realizations** (*int*) -- Number of realizations.
           - **n_filters** (*int*) -- Number of filters used.
           - **filters** (*int[]*) -- User-defined unique IDs of the filters on
             the tally
           - **nuclides** (*char[][]*) -- Array of nuclides to tally. Note that
             if no nuclide is specified in the user input, a single 'total'
             nuclide appears here.
           - **derivative** (*int*) -- ID of the derivative applied to the
             tally.
           - **n_score_bins** (*int*) -- Number of scoring bins for a single
             nuclide.
           - **score_bins** (*char[][]*) -- Values of specified scores.
           - **results** (*double[][][2]*) -- Accumulated sum and sum-of-squares
             for each bin of the i-th tally. The first dimension represents
             combinations of filter bins, the second dimensions represents
             scoring bins, and the third dimension has two entries for the sum
             and the sum-of-squares.

**/runtime/**

All values are given in seconds and are measured on the master process.

:Datasets: - **total initialization** (*double*) -- Time spent reading inputs,
             allocating arrays, etc.
           - **reading cross sections** (*double*) -- Time spent loading cross
             section libraries (this is a subset of initialization).
           - **simulation** (*double*) -- Time spent between initialization and
             finalization.
           - **transport** (*double*) -- Time spent transporting particles.
           - **inactive batches** (*double*) -- Time spent in the inactive
             batches (including non-transport activities like communcating
             sites).
           - **active batches** (*double*) -- Time spent in the active batches
             (including non-transport activities like communicating sites).
           - **synchronizing fission bank** (*double*) -- Time spent sampling
             source particles from fission sites and communicating them to other
             processes for load balancing.
           - **sampling source sites** (*double*) -- Time spent sampling source
             particles from fission sites.
           - **SEND-RECV source sites** (*double*) -- Time spent communicating
             source sites between processes for load balancing.
           - **accumulating tallies** (*double*) -- Time spent communicating
             tally results and evaluating their statistics.
