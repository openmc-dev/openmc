.. _usersguide_statepoint:

=======================
State Point File Format
=======================

The current revision of the statepoint file format is 15.

**/filetype** (*char[]*)

    String indicating the type of file.

**/revision** (*int*)

    Revision of the state point file format. Any time a change is made in the
    format, this integer is incremented.

**/version_major** (*int*)

    Major version number for OpenMC

**/version_minor** (*int*)

    Minor version number for OpenMC

**/version_release** (*int*)

    Release version number for OpenMC

**/date_and_time** (*char[]*)

    Date and time the state point was written.

**/path** (*char[]*)

    Absolute path to directory containing input files.

**/seed** (*int8_t*)

    Pseudo-random number generator seed.

**/run_CE** (*int*)

    Flag to denote continuous-energy or multi-group mode. A value of 1
    indicates a continuous-energy run while a value of 0 indicates a
    multi-group run.

**/run_mode** (*char[]*)

    Run mode used. A value of 1 indicates a fixed-source run and a value of 2
    indicates an eigenvalue run.

**/n_particles** (*int8_t*)

    Number of particles used per generation.

**/n_batches** (*int*)

    Number of batches to simulate.

**/current_batch** (*int*)

    The number of batches already simulated.

if run_mode == 'k-eigenvalue':

    **/n_inactive** (*int*)

        Number of inactive batches.

    **/gen_per_batch** (*int*)

        Number of generations per batch.

    **/k_generation** (*double[]*)

        k-effective for each generation simulated.

    **/entropy** (*double[]*)

        Shannon entropy for each generation simulated

    **/k_col_abs** (*double*)

        Sum of product of collision/absorption estimates of k-effective

    **/k_col_tra** (*double*)

        Sum of product of collision/track-length estimates of k-effective

    **/k_abs_tra** (*double*)

        Sum of product of absorption/track-length estimates of k-effective

    **/k_combined** (*double[2]*)

        Mean and standard deviation of a combined estimate of k-effective

    **/cmfd_on** (*int*)

        Flag indicating whether CMFD is on (1) or off (0).

    if (cmfd_on)

        **/cmfd/indices** (*int[4]*)

            Indices for cmfd mesh (i,j,k,g)

        **/cmfd/k_cmfd** (*double[]*)

            CMFD eigenvalues

        **/cmfd/cmfd_src** (*double[][][][]*)

            CMFD fission source

        **/cmfd/cmfd_entropy** (*double[]*)

            CMFD estimate of Shannon entropy

        **/cmfd/cmfd_balance** (*double[]*)

            RMS of the residual neutron balance equation on CMFD mesh

        **/cmfd/cmfd_dominance** (*double[]*)

            CMFD estimate of dominance ratio

        **/cmfd/cmfd_srccmp** (*double[]*)

            RMS comparison of difference between OpenMC and CMFD fission source

**/tallies/n_meshes** (*int*)

    Number of meshes in tallies.xml file

**/tally/meshes/ids** (*int[]*)

    Internal unique ID of each mesh.

**/tally/meshes/keys** (*int[]*)

    User-identified unique ID of each mesh.

**/tallies/meshes/mesh <uid>/type** (*char[]*)

    Type of mesh.

**/tallies/meshes/mesh <uid>/dimension** (*int*)

    Number of mesh cells in each dimension.

**/tallies/meshes/mesh <uid>/lower_left** (*double[]*)

    Coordinates of lower-left corner of mesh.

**/tallies/meshes/mesh <uid>/upper_right** (*double[]*)

    Coordinates of upper-right corner of mesh.

**/tallies/meshes/mesh <uid>/width** (*double[]*)

    Width of each mesh cell in each dimension.

**/tallies/n_tallies** (*int*)

    Number of user-defined tallies.

**/tallies/ids** (*int[]*)

    Internal unique ID of each tally.

**/tallies/keys** (*int[]*)

    User-identified unique ID of each tally.

**/tallies/tally <uid>/estimator** (*char[]*)

    Type of tally estimator, either 'analog', 'tracklength', or 'collision'.

**/tallies/tally <uid>/n_realizations** (*int*)

    Number of realizations.

**/tallies/tally <uid>/n_filters** (*int*)

    Number of filters used.

**/tallies/tally <uid>/filter <j>/type** (*char[]*)

    Type of the j-th filter. Can be 'universe', 'material', 'cell', 'cellborn',
    'surface', 'mesh', 'energy', 'energyout', or 'distribcell'.

**/tallies/tally <uid>/filter <j>/n_bins** (*int*)

    Number of bins for the j-th filter.

**/tallies/tally <uid>/filter <j>/bins** (*int[]* or *double[]*)

    Value for each filter bin of this type.

**/tallies/tally <uid>/nuclides** (*char[][]*)

    Array of nuclides to tally. Note that if no nuclide is specified in the user
    input, a single 'total' nuclide appears here.

**/tallies/tally <uid>/n_score_bins** (*int*)

    Number of scoring bins for a single nuclide. In general, this can be greater
    than the number of user-specified scores since each score might have
    multiple scoring bins, e.g., scatter-PN.

**/tallies/tally <uid>/score_bins** (*char[][]*)

    Values of specified scores.

**/tallies/tally <uid>/n_user_scores** (*int*)

    Number of scores without accounting for those added by expansions,
    e.g. scatter-PN.

**/tallies/tally <uid>/moment_orders** (*char[][]*)

    Tallying moment orders for Legendre and spherical harmonic tally expansions
    (*e.g.*, 'P2', 'Y1,2', etc.).

**/tallies/tally <uid>/results** (Compound type)

    Accumulated sum and sum-of-squares for each bin of the i-th tally. This is a
    two-dimensional array, the first dimension of which represents combinations
    of filter bins and the second dimensions of which represents scoring
    bins. Each element of the array has fields 'sum' and 'sum_sq'.

**/source_present** (*int*)

    Flag indicated if source bank is present in the file

**/n_realizations** (*int*)

    Number of realizations for global tallies.

**/n_global_tallies** (*int*)

    Number of global tally scores.

**/global_tallies** (Compound type)

    Accumulated sum and sum-of-squares for each global tally. The compound type
    has fields named ``sum`` and ``sum_sq``.

**/tallies_present** (*int*)

    Flag indicated if tallies are present in the file.

if (run_mode == 'k-eigenvalue' and source_present > 0)

    **/source_bank** (Compound type)

        Source bank information for each particle. The compound type has fields
        ``wgt``, ``xyz``, ``uvw``, ``E``, ``g``, and ``delayed_group``, which
        represent the weight, position, direction, energy, energy group, and
        delayed_group of the source particle, respectively.

**/runtime/total initialization** (*double*)

    Time (in seconds on the master process) spent reading inputs, allocating
    arrays, etc.

**/runtime/reading cross sections** (*double*)

    Time (in seconds on the master process) spent loading cross section
    libraries (this is a subset of initialization).

**/runtime/simulation** (*double*)

    Time (in seconds on the master process) spent between initialization and
    finalization.

**/runtime/transport** (*double*)

    Time (in seconds on the master process) spent transporting particles.

**/runtime/inactive batches** (*double*)

    Time (in seconds on the master process) spent in the inactive batches
    (including non-transport activities like communcating sites).

**/runtime/active batches** (*double*)

    Time (in seconds on the master process) spent in the active batches
    (including non-transport activities like communicating sites).

**/runtime/synchronizing fission bank** (*double*)

    Time (in seconds on the master process) spent sampling source particles
    from fission sites and communicating them to other processes for load
    balancing.

**/runtime/sampling source sites** (*double*)

    Time (in seconds on the master process) spent sampling source particles
    from fission sites.

**/runtime/SEND-RECV source sites** (*double*)

    Time (in seconds on the master process) spent communicating source sites
    between processes for load balancing.

**/runtime/accumulating tallies** (*double*)

    Time (in seconds on the master process) spent communicating tally results
    and evaluating their statistics.

**/runtime/CMFD** (*double*)

    Time (in seconds on the master process) spent evaluating CMFD.

**/runtime/CMFD building matrices** (*double*)

    Time (in seconds on the master process) spent buliding CMFD matrices.

**/runtime/CMFD solving matrices** (*double*)

    Time (in seconds on the master process) spent solving CMFD matrices.

**/runtime/total** (*double*)

    Total time spent (in seconds on the master process) in the program.
