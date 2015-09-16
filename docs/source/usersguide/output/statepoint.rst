.. _usersguide_statepoint:

=======================
State Point File Format
=======================

The current revision of the statepoint file format is 13.

**/filetype** (*int*)

    Flags what type of file this is. A value of -1 indicates a statepoint file,
    a value of -2 indicates a particle restart file, a value of -3 indicates a
    source file, and a value of -4 indicates a track file.

**/revision** (*int*)

    Revision of the state point file format. Any time a change is made in the
    format, this integer is incremented.

**/version_major** (*int*)

    Major version number for OpenMC

**/version_minor** (*int*)

    Minor version number for OpenMC

**/version_release** (*int*)

    Release version number for OpenMC

**/time_stamp** (*char[]*)

    Date and time the state point was written.

**/path** (*char[]*)

    Absolute path to directory containing input files.

**/seed** (*int8_t*)

    Pseudo-random number generator seed.

**/run_mode** (*int*)

    Run mode used. A value of 1 indicates a fixed-source run and a value of 2
    indicates an eigenvalue run.

**/n_particles** (*int8_t*)

    Number of particles used per generation.

**/n_batches** (*int*)

    Number of batches to simulate.

**/current_batch** (*int*)

    The number of batches already simulated.

if (run_mode == MODE_EIGENVALUE)

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

    User-identified unique ID of each mesh

*do i = 1, n_meshes*

    **/tallies/meshes/mesh i/id** (*int*)

        Unique identifier of the mesh.

    **/tallies/meshes/mesh i/type** (*int*)

        Type of mesh.

    **/tallies/meshes/mesh i/n_dimension** (*int*)

        Number of dimensions for mesh (2 or 3).

    **/tallies/meshes/mesh i/dimension** (*int*)

        Number of mesh cells in each dimension.

    **/tallies/meshes/mesh i/lower_left** (*double[]*)

        Coordinates of lower-left corner of mesh.

    **/tallies/meshes/mesh i/upper_right** (*double[]*)

        Coordinates of upper-right corner of mesh.

    **/tallies/meshes/mesh i/width** (*double[]*)

        Width of each mesh cell in each dimension.

**/tallies/n_tallies** (*int*)

    Number of user-defined tallies.

**/tallies/ids** (*int[]*)

    Internal unique ID of each tally.

**/tallies/keys** (*int[]*)

    User-identified unique ID of each tally.

*do i = 1, n_tallies*

    **/tallies/tally i/estimator** (*int*)

        Type of tally estimator: analog (1) or tracklength (2).

    **/tallies/tally i/n_realizations** (*int*)

        Number of realizations.

    **/tallies/tally i/n_filters** (*int*)

        Number of filters used.

    *do j = 1, tallies(i) % n_filters*

        **/tallies/tally i/filter j/type** (*int*)

            Type of tally filter.

        **/tallies/tally i/filter j/offset** (*int*)

            Filter offset (used for distribcell).

        **/tallies/tally i/filter j/n_bins** (*int*)

            Number of bins for filter.

        **/tallies/tally i/filter j/bins** (*int[]* or *double[]*)

            Value for each filter bin of this type.

    **/tallies/tally i/n_nuclides** (*int*)

        Number of nuclide bins. If none are specified, this is just one.

    **/tallies/tally i/nuclides** (*int[]*)

        Values of specified nuclide bins (ZAID identifiers)

    **/tallies/tally i/n_score_bins** (*int*)

        Number of scoring bins.

    **/tallies/tally i/score_bins** (*int*)

        Values of specified scoring bins (e.g. SCORE_FLUX).

    **/tallies/tally i/n_user_score_bins** (*int*)

        Number of scoring bins without accounting for those added by
        expansions, e.g. scatter-PN.

    **/tallies/tally i/moment_orders** (*char[][]*)

        Tallying moment orders for Legendre and spherical harmonic tally
        expansions (*e.g.*, 'P2', 'Y1,2', etc.).

**/source_present** (*int*)

    Flag indicated if source bank is present in the file

**/n_realizations** (*int*)

    Number of realizations for global tallies.

**/n_global_tallies** (*int*)

    Number of global tally scores.

**/global_tallies** (Compound type)

    Accumulated sum and sum-of-squares for each global tally. The compound type
    has fields named ``sum`` and ``sum_sq``.

**tallies_present** (*int*)

    Flag indicated if tallies are present in the file.

*do i = 1, n_tallies*

**/tallies/tally i/results** (Compound type)

    Accumulated sum and sum-of-squares for each bin of the tally i-th tally

if (run_mode == MODE_EIGENVALUE and source_present)

    **/source_bank** (Compound type)

        Source bank information for each particle. The compound type has fields
        ``wgt``, ``xyz``, ``uvw``, and ``E`` which represent the weight,
        position, direction, and energy of the source particle, respectively.
