.. _devguide_statepoint:

======================================
State Point Binary File Specifications
======================================

-----------
Revision 11
-----------

**integer(4) FILETYPE_STATEPOINT**

    Flags whether this file is a statepoint file or a particle restart file.

**integer(4) REVISION_STATEPOINT**

    Revision of the binary state point file. Any time a change is made in the
    format of the state-point file, this integer is incremented.

**integer(4) VERSION_MAJOR**

    Major version number for OpenMC

**integer(4) VERSION_MINOR**

    Minor version number for OpenMC

**integer(4) VERSION_RELEASE**

    Release version number for OpenMC

**character(19) time_stamp**

    Date and time the state point was written.

**character(255) path**

    Absolute path to directory containing input files.

**integer(8) seed**

    Pseudo-random number generator seed.

**integer(4) run_mode**

    run mode used. The modes are described in constants.F90.

**integer(8) n_particles**

    Number of particles used per generation.

**integer(4) n_batches**

    Total number of batches (active + inactive).

**integer(4) current_batch**

    The number of batches already simulated.

if (run_mode == MODE_EIGENVALUE)

    **integer(4) n_inactive**

        Number of inactive batches

    **integer(4) gen_per_batch**

        Number of generations per batch for criticality calculations

    *do i = 1, current_batch \* gen_per_batch*

        **real(8) k_generation(i)**

             k-effective for the i-th total generation

    *do i = 1, current_batch \* gen_per_batch*

        **real(8) entropy(i)**

            Shannon entropy for the i-th total generation

    **real(8) k_col_abs**

        Sum of product of collision/absorption estimates of k-effective

    **real(8) k_col_tra**

        Sum of product of collision/track-length estimates of k-effective

    **real(8) k_abs_tra**

        Sum of product of absorption/track-length estimates of k-effective

    **real(8) k_combined(2)**

        Mean and standard deviation of a combined estimate of k-effective

    **integer(4) cmfd_on**

        Flag that cmfd is on

    if (cmfd_on)

        **integer(4) cmfd % indices**

            Indices for cmfd mesh (i,j,k,g)

        **real(8) cmfd % k_cmfd(1:current_batch)**

            CMFD eigenvalues

        **real(8) cmfd % src(1:G,1:I,1:J,1:K)**

            CMFD fission source

        **real(8) cmfd % entropy(1:current_batch)**

            CMFD estimate of Shannon entropy

        **real(8) cmfd % balance(1:current_batch)**

            RMS of the residual neutron balance equation on CMFD mesh

        **real(8) cmfd % dom(1:current_batch)**

            CMFD estimate of dominance ratio

        **real(8) cmfd % scr_cmp(1:current_batch)**

            RMS comparison of difference between OpenMC and CMFD fission source

**integer(4) n_meshes**

    Number of meshes in tallies.xml file

*do i = 1, n_meshes*

    **integer(4) meshes(i) % id**

        Unique ID of mesh.

    **integer(4) meshes(i) % type**

        Type of mesh.

    **integer(4) meshes(i) % n_dimension**

        Number of dimensions for mesh (2 or 3).

    **integer(4) meshes(i) % dimension(:)**

        Number of mesh cells in each dimension.

    **real(8) meshes(i) % lower_left(:)**

        Coordinates of lower-left corner of mesh.

    **real(8) meshes(i) % upper_right(:)**

        Coordinates of upper-right corner of mesh.

    **real(8) meshes(i) % width(:)**

        Width of each mesh cell in each dimension.

**integer(4) n_tallies**

*do i = 1, n_tallies*

    **integer(4) tallies(i) % id**

        Unique ID of tally.

    **integer(4) tallies(i) % n_realizations**

        Number of realizations for the i-th tally.

    **integer(4) size(tallies(i) % scores, 1)**

        Total number of score bins for the i-th tally

    **integer(4) size(tallies(i) % scores, 2)**

        Total number of filter bins for the i-th tally

    **integer(4) tallies(i) % n_filters**

    *do j = 1, tallies(i) % n_filters*

        **integer(4) tallies(i) % filter(j) % type**

            Type of tally filter.

        **integer(4) tallies(i) % filter(j) % n_bins**

            Number of bins for filter.

        **integer(4)/real(8) tallies(i) % filter(j) % bins(:)**

            Value for each filter bin of this type.

    **integer(4) tallies(i) % n_nuclide_bins**

        Number of nuclide bins. If none are specified, this is just one.

    *do j = 1, tallies(i) % n_nuclide_bins*

        **integer(4) tallies(i) % nuclide_bins(j)**

            Values of specified nuclide bins

    **integer(4) tallies(i) % n_score_bins**

        Number of scoring bins.

    *do j = 1, tallies(i) % n_score_bins*

        **integer(4) tallies(i) % score_bins(j)**

            Values of specified scoring bins (e.g. SCORE_FLUX).
            
    *do j = 1, tallies(i) % n_score_bins*

        **integer(4) tallies(i) % scatt_order(j)**

            Scattering Order specified scoring bins.
            
    **integer(4) tallies(i) % n_score_bins**

        Number of scoring bins without accounting for those added by
        the scatter-pn command.

**integer(4) source_present**

    Flag indicated if source bank is present in the file

**integer(4) n_realizations**

    Number of realizations for global tallies.

**integer(4) N_GLOBAL_TALLIES**

    Number of global tally scores

*do i = 1, N_GLOBAL_TALLIES*

    **real(8) global_tallies(i) % sum**

        Accumulated sum for the i-th global tally

    **real(8) global_tallies(i) % sum_sq**

        Accumulated sum of squares for the i-th global tally

**integer(4) tallies_on**

    Flag indicated if tallies are present in the file.

if (tallies_on > 0)

    *do i = 1, n_tallies*

        *do k = 1, size(tallies(i) % scores, 2)*

            *do j = 1, size(tallies(i) % scores, 1)*

                **real(8) tallies(i) % scores(j,k) % sum**
            
                    Accumulated sum for the j-th score and k-th filter of the
                    i-th tally

                **real(8) tallies(i) % scores(j,k) % sum_sq**

                    Accumulated sum of squares for the j-th score and k-th
                    filter of the i-th tally

if (run_mode == MODE_EIGENVALUE and source_present)

    *do i = 1, n_particles*

        **real(8) source_bank(i) % wgt**

            Weight of the i-th source particle

        **real(8) source_bank(i) % xyz(1:3)**

            Coordinates of the i-th source particle.

        **real(8) source_bank(i) % uvw(1:3)**

            Direction of the i-th source particle

        **real(8) source_bank(i) % E**

            Energy of the i-th source particle.

-----------
Revision 10 
-----------

**integer(4) FILETYPE_STATEPOINT**

    Flags whether this file is a statepoint file or a particle restart file.

**integer(4) REVISION_STATEPOINT**

    Revision of the binary state point file. Any time a change is made in the
    format of the state-point file, this integer is incremented.

**integer(4) VERSION_MAJOR**

    Major version number for OpenMC

**integer(4) VERSION_MINOR**

    Minor version number for OpenMC

**integer(4) VERSION_RELEASE**

    Release version number for OpenMC

**character(19) time_stamp**

    Date and time the state point was written.

**character(255) path**

    Absolute path to directory containing input files.

**integer(8) seed**

    Pseudo-random number generator seed.

**integer(4) run_mode**

    run mode used. The modes are described in constants.F90.

**integer(8) n_particles**

    Number of particles used per generation.

**integer(4) n_batches**

    Total number of batches (active + inactive).

**integer(4) current_batch**

    The number of batches already simulated.

if (run_mode == MODE_EIGENVALUE)

    **integer(4) n_inactive**

        Number of inactive batches

    **integer(4) gen_per_batch**

        Number of generations per batch for criticality calculations

    *do i = 1, current_batch \* gen_per_batch*

        **real(8) k_generation(i)**

             k-effective for the i-th total generation

    *do i = 1, current_batch \* gen_per_batch*

        **real(8) entropy(i)**

            Shannon entropy for the i-th total generation

    **real(8) k_col_abs**

        Sum of product of collision/absorption estimates of k-effective

    **real(8) k_col_tra**

        Sum of product of collision/track-length estimates of k-effective

    **real(8) k_abs_tra**

        Sum of product of absorption/track-length estimates of k-effective

    **real(8) k_combined(2)**

        Mean and standard deviation of a combined estimate of k-effective

    **integer(4) cmfd_on**

        Flag that cmfd is on

    if (cmfd_on)

        **integer(4) cmfd % indices**

            Indices for cmfd mesh (i,j,k,g)

        **real(8) cmfd % k_cmfd(1:current_batch)**

            CMFD eigenvalues

        **real(8) cmfd % src(1:I,1:J,1:K,1:G)**

            CMFD fission source

        **real(8) cmfd % entropy(1:current_batch)**

            CMFD estimate of Shannon entropy

        **real(8) cmfd % balance(1:current_batch)**

            RMS of the residual neutron balance equation on CMFD mesh

        **real(8) cmfd % dom(1:current_batch)**

            CMFD estimate of dominance ratio

        **real(8) cmfd % scr_cmp(1:current_batch)**

            RMS comparison of difference between OpenMC and CMFD fission source

**integer(4) n_meshes**

    Number of meshes in tallies.xml file

*do i = 1, n_meshes*

    **integer(4) meshes(i) % id**

        Unique ID of mesh.

    **integer(4) meshes(i) % type**

        Type of mesh.

    **integer(4) meshes(i) % n_dimension**

        Number of dimensions for mesh (2 or 3).

    **integer(4) meshes(i) % dimension(:)**

        Number of mesh cells in each dimension.

    **real(8) meshes(i) % lower_left(:)**

        Coordinates of lower-left corner of mesh.

    **real(8) meshes(i) % upper_right(:)**

        Coordinates of upper-right corner of mesh.

    **real(8) meshes(i) % width(:)**

        Width of each mesh cell in each dimension.

**integer(4) n_tallies**

*do i = 1, n_tallies*

    **integer(4) tallies(i) % id**

        Unique ID of tally.

    **integer(4) tallies(i) % n_realizations**

        Number of realizations for the i-th tally.

    **integer(4) size(tallies(i) % scores, 1)**

        Total number of score bins for the i-th tally

    **integer(4) size(tallies(i) % scores, 2)**

        Total number of filter bins for the i-th tally

    **integer(4) tallies(i) % n_filters**

    *do j = 1, tallies(i) % n_filters*

        **integer(4) tallies(i) % filter(j) % type**

            Type of tally filter.

        **integer(4) tallies(i) % filter(j) % n_bins**

            Number of bins for filter.

        **integer(4)/real(8) tallies(i) % filter(j) % bins(:)**

            Value for each filter bin of this type.

    **integer(4) tallies(i) % n_nuclide_bins**

        Number of nuclide bins. If none are specified, this is just one.

    *do j = 1, tallies(i) % n_nuclide_bins*

        **integer(4) tallies(i) % nuclide_bins(j)**

            Values of specified nuclide bins

    **integer(4) tallies(i) % n_score_bins**

        Number of scoring bins.

    *do j = 1, tallies(i) % n_score_bins*

        **integer(4) tallies(i) % score_bins(j)**

            Values of specified scoring bins (e.g. SCORE_FLUX).
            
    *do j = 1, tallies(i) % n_score_bins*

        **integer(4) tallies(i) % scatt_order(j)**

            Scattering Order specified scoring bins.
            
    **integer(4) tallies(i) % n_score_bins**

        Number of scoring bins without accounting for those added by
        the scatter-pn command.

**integer(4) n_realizations**

    Number of realizations for global tallies.

**integer(4) N_GLOBAL_TALLIES**

    Number of global tally scores

*do i = 1, N_GLOBAL_TALLIES*

    **real(8) global_tallies(i) % sum**

        Accumulated sum for the i-th global tally

    **real(8) global_tallies(i) % sum_sq**

        Accumulated sum of squares for the i-th global tally

**integer(4) tallies_on**

    Flag indicated if tallies are present in the file.

if (tallies_on > 0)

    *do i = 1, n_tallies*

        *do k = 1, size(tallies(i) % scores, 2)*

            *do j = 1, size(tallies(i) % scores, 1)*

                **real(8) tallies(i) % scores(j,k) % sum**
            
                    Accumulated sum for the j-th score and k-th filter of the
                    i-th tally

                **real(8) tallies(i) % scores(j,k) % sum_sq**

                    Accumulated sum of squares for the j-th score and k-th
                    filter of the i-th tally

if (run_mode == MODE_EIGENVALUE)

    *do i = 1, n_particles*

        **real(8) source_bank(i) % wgt**

            Weight of the i-th source particle

        **real(8) source_bank(i) % xyz(1:3)**

            Coordinates of the i-th source particle.

        **real(8) source_bank(i) % uvw(1:3)**

            Direction of the i-th source particle

        **real(8) source_bank(i) % E**

            Energy of the i-th source particle.

----------
Revision 9
----------

**integer(4) FILETYPE_STATEPOINT**

    Flags whether this file is a statepoint file or a particle restart file.

**integer(4) REVISION_STATEPOINT**

    Revision of the binary state point file. Any time a change is made in the
    format of the state-point file, this integer is incremented.

**integer(4) VERSION_MAJOR**

    Major version number for OpenMC

**integer(4) VERSION_MINOR**

    Minor version number for OpenMC

**integer(4) VERSION_RELEASE**

    Release version number for OpenMC

**character(19) time_stamp**

    Date and time the state point was written.

**character(255) path**

    Absolute path to directory containing input files.

**integer(8) seed**

    Pseudo-random number generator seed.

**integer(4) run_mode**

    run mode used. The modes are described in constants.F90.

**integer(8) n_particles**

    Number of particles used per generation.

**integer(4) n_batches**

    Total number of batches (active + inactive).

**integer(4) current_batch**

    The number of batches already simulated.

if (run_mode == MODE_EIGENVALUE)

    **integer(4) n_inactive**

        Number of inactive batches

    **integer(4) gen_per_batch**

        Number of generations per batch for criticality calculations

    *do i = 1, current_batch \* gen_per_batch*

        **real(8) k_generation(i)**

             k-effective for the i-th total generation

    *do i = 1, current_batch \* gen_per_batch*

        **real(8) entropy(i)**

            Shannon entropy for the i-th total generation

    **real(8) k_col_abs**

        Sum of product of collision/absorption estimates of k-effective

    **real(8) k_col_tra**

        Sum of product of collision/track-length estimates of k-effective

    **real(8) k_abs_tra**

        Sum of product of absorption/track-length estimates of k-effective

    **real(8) k_combined(2)**

        Mean and standard deviation of a combined estimate of k-effective

**integer(4) n_meshes**

    Number of meshes in tallies.xml file

*do i = 1, n_meshes*

    **integer(4) meshes(i) % id**

        Unique ID of mesh.

    **integer(4) meshes(i) % type**

        Type of mesh.

    **integer(4) meshes(i) % n_dimension**

        Number of dimensions for mesh (2 or 3).

    **integer(4) meshes(i) % dimension(:)**

        Number of mesh cells in each dimension.

    **real(8) meshes(i) % lower_left(:)**

        Coordinates of lower-left corner of mesh.

    **real(8) meshes(i) % upper_right(:)**

        Coordinates of upper-right corner of mesh.

    **real(8) meshes(i) % width(:)**

        Width of each mesh cell in each dimension.

**integer(4) n_tallies**

*do i = 1, n_tallies*

    **integer(4) tallies(i) % id**

        Unique ID of tally.

    **integer(4) tallies(i) % n_realizations**

        Number of realizations for the i-th tally.

    **integer(4) size(tallies(i) % scores, 1)**

        Total number of score bins for the i-th tally

    **integer(4) size(tallies(i) % scores, 2)**

        Total number of filter bins for the i-th tally

    **integer(4) tallies(i) % n_filters**

    *do j = 1, tallies(i) % n_filters*

        **integer(4) tallies(i) % filter(j) % type**

            Type of tally filter.

        **integer(4) tallies(i) % filter(j) % n_bins**

            Number of bins for filter.

        **integer(4)/real(8) tallies(i) % filter(j) % bins(:)**

            Value for each filter bin of this type.

    **integer(4) tallies(i) % n_nuclide_bins**

        Number of nuclide bins. If none are specified, this is just one.

    *do j = 1, tallies(i) % n_nuclide_bins*

        **integer(4) tallies(i) % nuclide_bins(j)**

            Values of specified nuclide bins

    **integer(4) tallies(i) % n_score_bins**

        Number of scoring bins.

    *do j = 1, tallies(i) % n_score_bins*

        **integer(4) tallies(i) % score_bins(j)**

            Values of specified scoring bins (e.g. SCORE_FLUX).
            
    *do j = 1, tallies(i) % n_score_bins*

        **integer(4) tallies(i) % scatt_order(j)**

            Scattering Order specified scoring bins.
            
    **integer(4) tallies(i) % n_score_bins**

        Number of scoring bins without accounting for those added by
        the scatter-pn command.

**integer(4) n_realizations**

    Number of realizations for global tallies.

**integer(4) N_GLOBAL_TALLIES**

    Number of global tally scores

*do i = 1, N_GLOBAL_TALLIES*

    **real(8) global_tallies(i) % sum**

        Accumulated sum for the i-th global tally

    **real(8) global_tallies(i) % sum_sq**

        Accumulated sum of squares for the i-th global tally

**integer(4) tallies_on**

    Flag indicated if tallies are present in the file.

if (tallies_on > 0)

    *do i = 1, n_tallies*

        *do k = 1, size(tallies(i) % scores, 2)*

            *do j = 1, size(tallies(i) % scores, 1)*

                **real(8) tallies(i) % scores(j,k) % sum**
            
                    Accumulated sum for the j-th score and k-th filter of the
                    i-th tally

                **real(8) tallies(i) % scores(j,k) % sum_sq**

                    Accumulated sum of squares for the j-th score and k-th
                    filter of the i-th tally

if (run_mode == MODE_EIGENVALUE)

    *do i = 1, n_particles*

        **real(8) source_bank(i) % wgt**

            Weight of the i-th source particle

        **real(8) source_bank(i) % xyz(1:3)**

            Coordinates of the i-th source particle.

        **real(8) source_bank(i) % uvw(1:3)**

            Direction of the i-th source particle

        **real(8) source_bank(i) % E**

            Energy of the i-th source particle.

----------
Revision 8
----------

**integer(4) REVISION_STATEPOINT**

    Revision of the binary state point file. Any time a change is made in the
    format of the state-point file, this integer is incremented.

**integer(4) VERSION_MAJOR**

    Major version number for OpenMC

**integer(4) VERSION_MINOR**

    Minor version number for OpenMC

**integer(4) VERSION_RELEASE**

    Release version number for OpenMC

**character(19) time_stamp**

    Date and time the state point was written.

**character(255) path**

    Absolute path to directory containing input files.

**integer(8) seed**

    Pseudo-random number generator seed.

**integer(4) run_mode**

    run mode used. The modes are described in constants.F90.

**integer(8) n_particles**

    Number of particles used per generation.

**integer(4) n_batches**

    Total number of batches (active + inactive).

**integer(4) current_batch**

    The number of batches already simulated.

if (run_mode == MODE_EIGENVALUE)

    **integer(4) n_inactive**

        Number of inactive batches

    **integer(4) gen_per_batch**

        Number of generations per batch for criticality calculations

    *do i = 1, current_batch*

        **real(8) k_batch(i)**

             k-effective for the i-th batch

    *do i = 1, current_batch \* gen_per_batch*

        **real(8) entropy(i)**

            Shannon entropy for the i-th batch

    **real(8) k_col_abs**

        Sum of product of collision/absorption estimates of k-effective

    **real(8) k_col_tra**

        Sum of product of collision/track-length estimates of k-effective

    **real(8) k_abs_tra**

        Sum of product of absorption/track-length estimates of k-effective

    **real(8) k_combined(2)**

        Mean and standard deviation of a combined estimate of k-effective

**integer(4) n_meshes**

    Number of meshes in tallies.xml file

*do i = 1, n_meshes*

    **integer(4) meshes(i) % id**

        Unique ID of mesh.

    **integer(4) meshes(i) % type**

        Type of mesh.

    **integer(4) meshes(i) % n_dimension**

        Number of dimensions for mesh (2 or 3).

    **integer(4) meshes(i) % dimension(:)**

        Number of mesh cells in each dimension.

    **real(8) meshes(i) % lower_left(:)**

        Coordinates of lower-left corner of mesh.

    **real(8) meshes(i) % upper_right(:)**

        Coordinates of upper-right corner of mesh.

    **real(8) meshes(i) % width(:)**

        Width of each mesh cell in each dimension.

**integer(4) n_tallies**

*do i = 1, n_tallies*

    **integer(4) tallies(i) % id**

        Unique ID of tally.

    **integer(4) tallies(i) % n_realizations**

        Number of realizations for the i-th tally.

    **integer(4) size(tallies(i) % scores, 1)**

        Total number of score bins for the i-th tally

    **integer(4) size(tallies(i) % scores, 2)**

        Total number of filter bins for the i-th tally

    **integer(4) tallies(i) % n_filters**

    *do j = 1, tallies(i) % n_filters*

        **integer(4) tallies(i) % filter(j) % type**

            Type of tally filter.

        **integer(4) tallies(i) % filter(j) % n_bins**

            Number of bins for filter.

        **integer(4)/real(8) tallies(i) % filter(j) % bins(:)**

            Value for each filter bin of this type.

    **integer(4) tallies(i) % n_nuclide_bins**

        Number of nuclide bins. If none are specified, this is just one.

    *do j = 1, tallies(i) % n_nuclide_bins*

        **integer(4) tallies(i) % nuclide_bins(j)**

            Values of specified nuclide bins

    **integer(4) tallies(i) % n_score_bins**

        Number of scoring bins.

    *do j = 1, tallies(i) % n_score_bins*

        **integer(4) tallies(i) % score_bins(j)**

            Values of specified scoring bins (e.g. SCORE_FLUX).
            
    *do j = 1, tallies(i) % n_score_bins*

        **integer(4) tallies(i) % scatt_order(j)**

            Scattering Order specified scoring bins.
            
    **integer(4) tallies(i) % n_score_bins**

        Number of scoring bins without accounting for those added by
        the scatter-pn command.

**integer(4) n_realizations**

    Number of realizations for global tallies.

**integer(4) N_GLOBAL_TALLIES**

    Number of global tally scores

*do i = 1, N_GLOBAL_TALLIES*

    **real(8) global_tallies(i) % sum**

        Accumulated sum for the i-th global tally

    **real(8) global_tallies(i) % sum_sq**

        Accumulated sum of squares for the i-th global tally

**integer(4) tallies_on**

    Flag indicated if tallies are present in the file.

if (tallies_on > 0)

    *do i = 1, n_tallies*

        *do k = 1, size(tallies(i) % scores, 2)*

            *do j = 1, size(tallies(i) % scores, 1)*

                **real(8) tallies(i) % scores(j,k) % sum**
            
                    Accumulated sum for the j-th score and k-th filter of the
                    i-th tally

                **real(8) tallies(i) % scores(j,k) % sum_sq**

                    Accumulated sum of squares for the j-th score and k-th
                    filter of the i-th tally

if (run_mode == MODE_EIGENVALUE)

    *do i = 1, n_particles*

        **real(8) source_bank(i) % wgt**

            Weight of the i-th source particle

        **real(8) source_bank(i) % xyz(1:3)**

            Coordinates of the i-th source particle.

        **real(8) source_bank(i) % uvw(1:3)**

            Direction of the i-th source particle

        **real(8) source_bank(i) % E**

            Energy of the i-th source particle.

----------
Revision 7
----------

**integer(4) REVISION_STATEPOINT**

    Revision of the binary state point file. Any time a change is made in the
    format of the state-point file, this integer is incremented.

**integer(4) VERSION_MAJOR**

    Major version number for OpenMC

**integer(4) VERSION_MINOR**

    Minor version number for OpenMC

**integer(4) VERSION_RELEASE**

    Release version number for OpenMC

**character(19) time_stamp**

    Date and time the state point was written.

**character(255) path**

    Absolute path to directory containing input files.

**integer(8) seed**

    Pseudo-random number generator seed.

**integer(4) run_mode**

    run mode used. The modes are described in constants.F90.

**integer(8) n_particles**

    Number of particles used per generation.

**integer(4) n_batches**

    Total number of batches (active + inactive).

**integer(4) current_batch**

    The number of batches already simulated.

if (run_mode == MODE_EIGENVALUE)

    **integer(4) n_inactive**

        Number of inactive batches

    **integer(4) gen_per_batch**

        Number of generations per batch for criticality calculations

    *do i = 1, current_batch*

        **real(8) k_batch(i)**

             k-effective for the i-th batch

    *do i = 1, current_batch \* gen_per_batch*

        **real(8) entropy(i)**

            Shannon entropy for the i-th batch

**integer(4) n_meshes**

    Number of meshes in tallies.xml file

*do i = 1, n_meshes*

    **integer(4) meshes(i) % id**

        Unique ID of mesh.

    **integer(4) meshes(i) % type**

        Type of mesh.

    **integer(4) meshes(i) % n_dimension**

        Number of dimensions for mesh (2 or 3).

    **integer(4) meshes(i) % dimension(:)**

        Number of mesh cells in each dimension.

    **real(8) meshes(i) % lower_left(:)**

        Coordinates of lower-left corner of mesh.

    **real(8) meshes(i) % upper_right(:)**

        Coordinates of upper-right corner of mesh.

    **real(8) meshes(i) % width(:)**

        Width of each mesh cell in each dimension.

**integer(4) n_tallies**

*do i = 1, n_tallies*

    **integer(4) tallies(i) % id**

        Unique ID of tally.

    **integer(4) tallies(i) % n_realizations**

        Number of realizations for the i-th tally.

    **integer(4) size(tallies(i) % scores, 1)**

        Total number of score bins for the i-th tally

    **integer(4) size(tallies(i) % scores, 2)**

        Total number of filter bins for the i-th tally

    **integer(4) tallies(i) % n_filters**

    *do j = 1, tallies(i) % n_filters*

        **integer(4) tallies(i) % filter(j) % type**

            Type of tally filter.

        **integer(4) tallies(i) % filter(j) % n_bins**

            Number of bins for filter.

        **integer(4)/real(8) tallies(i) % filter(j) % bins(:)**

            Value for each filter bin of this type.

    **integer(4) tallies(i) % n_nuclide_bins**

        Number of nuclide bins. If none are specified, this is just one.

    *do j = 1, tallies(i) % n_nuclide_bins*

        **integer(4) tallies(i) % nuclide_bins(j)**

            Values of specified nuclide bins

    **integer(4) tallies(i) % n_score_bins**

        Number of scoring bins.

    *do j = 1, tallies(i) % n_score_bins*

        **integer(4) tallies(i) % score_bins(j)**

            Values of specified scoring bins (e.g. SCORE_FLUX).
            
    *do j = 1, tallies(i) % n_score_bins*

        **integer(4) tallies(i) % scatt_order(j)**

            Scattering Order specified scoring bins.
            
    **integer(4) tallies(i) % n_score_bins**

        Number of scoring bins without accounting for those added by
        the scatter-pn command.

**integer(4) n_realizations**

    Number of realizations for global tallies.

**integer(4) N_GLOBAL_TALLIES**

    Number of global tally scores

*do i = 1, N_GLOBAL_TALLIES*

    **real(8) global_tallies(i) % sum**

        Accumulated sum for the i-th global tally

    **real(8) global_tallies(i) % sum_sq**

        Accumulated sum of squares for the i-th global tally

**integer(4) tallies_on**

    Flag indicated if tallies are present in the file.

if (tallies_on > 0)

    *do i = 1, n_tallies*

        *do k = 1, size(tallies(i) % scores, 2)*

            *do j = 1, size(tallies(i) % scores, 1)*

                **real(8) tallies(i) % scores(j,k) % sum**
            
                    Accumulated sum for the j-th score and k-th filter of the
                    i-th tally

                **real(8) tallies(i) % scores(j,k) % sum_sq**

                    Accumulated sum of squares for the j-th score and k-th
                    filter of the i-th tally

if (run_mode == MODE_EIGENVALUE)

    *do i = 1, n_particles*

        **real(8) source_bank(i) % wgt**

            Weight of the i-th source particle

        **real(8) source_bank(i) % xyz(1:3)**

            Coordinates of the i-th source particle.

        **real(8) source_bank(i) % uvw(1:3)**

            Direction of the i-th source particle

        **real(8) source_bank(i) % E**

            Energy of the i-th source particle.

----------
Revision 6
----------

**integer(4) REVISION_STATEPOINT**

    Revision of the binary state point file. Any time a change is made in the
    format of the state-point file, this integer is incremented.

**integer(4) VERSION_MAJOR**

    Major version number for OpenMC

**integer(4) VERSION_MINOR**

    Minor version number for OpenMC

**integer(4) VERSION_RELEASE**

    Release version number for OpenMC

**character(19) time_stamp**

    Date and time the state point was written.

**character(255) path**

    Absolute path to directory containing input files.

**integer(8) seed**

    Pseudo-random number generator seed.

**integer(4) run_mode**

    run mode used. The modes are described in constants.F90.

**integer(8) n_particles**

    Number of particles used per generation.

**integer(4) n_batches**

    Total number of batches (active + inactive).

**integer(4) current_batch**

    The number of batches already simulated.

if (run_mode == MODE_EIGENVALUE)

    **integer(4) n_inactive**

        Number of inactive batches

    **integer(4) gen_per_batch**

        Number of generations per batch for criticality calculations

    *do i = 1, current_batch*

        **real(8) k_batch(i)**

             k-effective for the i-th batch

    *do i = 1, current_batch*

        **real(8) entropy(i)**

            Shannon entropy for the i-th batch

**integer(4) n_meshes**

    Number of meshes in tallies.xml file

*do i = 1, n_meshes*

    **integer(4) meshes(i) % id**

        Unique ID of mesh.

    **integer(4) meshes(i) % type**

        Type of mesh.

    **integer(4) meshes(i) % n_dimension**

        Number of dimensions for mesh (2 or 3).

    **integer(4) meshes(i) % dimension(:)**

        Number of mesh cells in each dimension.

    **real(8) meshes(i) % lower_left(:)**

        Coordinates of lower-left corner of mesh.

    **real(8) meshes(i) % upper_right(:)**

        Coordinates of upper-right corner of mesh.

    **real(8) meshes(i) % width(:)**

        Width of each mesh cell in each dimension.

**integer(4) n_tallies**

*do i = 1, n_tallies*

    **integer(4) tallies(i) % id**

        Unique ID of tally.

    **integer(4) tallies(i) % n_realizations**

        Number of realizations for the i-th tally.

    **integer(4) size(tallies(i) % scores, 1)**

        Total number of score bins for the i-th tally

    **integer(4) size(tallies(i) % scores, 2)**

        Total number of filter bins for the i-th tally

    **integer(4) tallies(i) % n_filters**

    *do j = 1, tallies(i) % n_filters*

        **integer(4) tallies(i) % filter(j) % type**

            Type of tally filter.

        **integer(4) tallies(i) % filter(j) % n_bins**

            Number of bins for filter.

        **integer(4)/real(8) tallies(i) % filter(j) % bins(:)**

            Value for each filter bin of this type.

    **integer(4) tallies(i) % n_nuclide_bins**

        Number of nuclide bins. If none are specified, this is just one.

    *do j = 1, tallies(i) % n_nuclide_bins*

        **integer(4) tallies(i) % nuclide_bins(j)**

            Values of specified nuclide bins

    **integer(4) tallies(i) % n_score_bins**

        Number of scoring bins.

    *do j = 1, tallies(i) % n_score_bins*

        **integer(4) tallies(i) % score_bins(j)**

            Values of specified scoring bins (e.g. SCORE_FLUX).

**integer(4) n_realizations**

    Number of realizations for global tallies.

**integer(4) N_GLOBAL_TALLIES**

    Number of global tally scores

*do i = 1, N_GLOBAL_TALLIES*

    **real(8) global_tallies(i) % sum**

        Accumulated sum for the i-th global tally

    **real(8) global_tallies(i) % sum_sq**

        Accumulated sum of squares for the i-th global tally

**integer(4) tallies_on**

    Flag indicated if tallies are present in the file.

if (tallies_on > 0)

    *do i = 1, n_tallies*

        *do k = 1, size(tallies(i) % scores, 2)*

            *do j = 1, size(tallies(i) % scores, 1)*

                **real(8) tallies(i) % scores(j,k) % sum**
            
                    Accumulated sum for the j-th score and k-th filter of the
                    i-th tally

                **real(8) tallies(i) % scores(j,k) % sum_sq**

                    Accumulated sum of squares for the j-th score and k-th
                    filter of the i-th tally

if (run_mode == MODE_EIGENVALUE)

    *do i = 1, n_particles*

        **real(8) source_bank(i) % wgt**

            Weight of the i-th source particle

        **real(8) source_bank(i) % xyz(1:3)**

            Coordinates of the i-th source particle.

        **real(8) source_bank(i) % uvw(1:3)**

            Direction of the i-th source particle

        **real(8) source_bank(i) % E**

            Energy of the i-th source particle.

----------
Revision 5
----------

**integer(4) REVISION_STATEPOINT**

    Revision of the binary state point file. Any time a change is made in the
    format of the state-point file, this integer is incremented.

**integer(4) VERSION_MAJOR**

    Major version number for OpenMC

**integer(4) VERSION_MINOR**

    Minor version number for OpenMC

**integer(4) VERSION_RELEASE**

    Release version number for OpenMC

**character(19) time_stamp**

    Date and time the state point was written.

**integer(8) seed**

    Pseudo-random number generator seed.

**integer(4) run_mode**

    run mode used. The modes are described in constants.F90.

**integer(8) n_particles**

    Number of particles used per generation.

**integer(4) n_batches**

    Total number of batches (active + inactive).

**integer(4) current_batch**

    The number of batches already simulated.

if (run_mode == MODE_EIGENVALUE)

    **integer(4) n_inactive**

        Number of inactive batches

    **integer(4) gen_per_batch**

        Number of generations per batch for criticality calculations

    *do i = 1, current_batch*

        **real(8) k_batch(i)**

             k-effective for the i-th batch

    *do i = 1, current_batch*

        **real(8) entropy(i)**

            Shannon entropy for the i-th batch

**integer(4) n_meshes**

    Number of meshes in tallies.xml file

*do i = 1, n_meshes*

    **integer(4) meshes(i) % type**

        Type of mesh.

    **integer(4) meshes(i) % n_dimension**

        Number of dimensions for mesh (2 or 3).

    **integer(4) meshes(i) % dimension(:)**

        Number of mesh cells in each dimension.

    **real(8) meshes(i) % lower_left(:)**

        Coordinates of lower-left corner of mesh.

    **real(8) meshes(i) % upper_right(:)**

        Coordinates of upper-right corner of mesh.

    **real(8) meshes(i) % width(:)**

        Width of each mesh cell in each dimension.

**integer(4) n_tallies**

*do i = 1, n_tallies*

    **integer(4) tallies(i) % n_realizations**

        Number of realizations for the i-th tally.

    **integer(4) size(tallies(i) % scores, 1)**

        Total number of score bins for the i-th tally

    **integer(4) size(tallies(i) % scores, 2)**

        Total number of filter bins for the i-th tally

    **integer(4) tallies(i) % n_filters**

    *do j = 1, tallies(i) % n_filters*

        **integer(4) tallies(i) % filter(j) % type**

            Type of tally filter.

        **integer(4) tallies(i) % filter(j) % n_bins**

            Number of bins for filter.

        **integer(4)/real(8) tallies(i) % filter(j) % bins(:)**

            Value for each filter bin of this type.

    **integer(4) tallies(i) % n_nuclide_bins**

        Number of nuclide bins. If none are specified, this is just one.

    *do j = 1, tallies(i) % n_nuclide_bins*

        **integer(4) tallies(i) % nuclide_bins(j)**

            Values of specified nuclide bins

    **integer(4) tallies(i) % n_score_bins**

        Number of scoring bins.

    *do j = 1, tallies(i) % n_score_bins*

        **integer(4) tallies(i) % score_bins(j)**

            Values of specified scoring bins (e.g. SCORE_FLUX).

**integer(4) n_realizations**

    Number of realizations for global tallies.

**integer(4) N_GLOBAL_TALLIES**

    Number of global tally scores

*do i = 1, N_GLOBAL_TALLIES*

    **real(8) global_tallies(i) % sum**

        Accumulated sum for the i-th global tally

    **real(8) global_tallies(i) % sum_sq**

        Accumulated sum of squares for the i-th global tally

**integer(4) tallies_on**

    Flag indicated if tallies are present in the file.

if (tallies_on > 0)

    *do i = 1, n_tallies*

        *do k = 1, size(tallies(i) % scores, 2)*

            *do j = 1, size(tallies(i) % scores, 1)*

                **real(8) tallies(i) % scores(j,k) % sum**
            
                    Accumulated sum for the j-th score and k-th filter of the
                    i-th tally

                **real(8) tallies(i) % scores(j,k) % sum_sq**

                    Accumulated sum of squares for the j-th score and k-th
                    filter of the i-th tally

if (run_mode == MODE_EIGENVALUE)

    *do i = 1, n_particles*

        **real(8) source_bank(i) % wgt**

            Weight of the i-th source particle

        **real(8) source_bank(i) % xyz(1:3)**

            Coordinates of the i-th source particle.

        **real(8) source_bank(i) % uvw(1:3)**

            Direction of the i-th source particle

        **real(8) source_bank(i) % E**

            Energy of the i-th source particle.

----------
Revision 4
----------

**integer(4) REVISION_STATEPOINT**

    Revision of the binary state point file. Any time a change is made in the
    format of the state-point file, this integer is incremented.

**integer(4) VERSION_MAJOR**

    Major version number for OpenMC

**integer(4) VERSION_MINOR**

    Minor version number for OpenMC

**integer(4) VERSION_RELEASE**

    Release version number for OpenMC

**character(19) time_stamp**

    Date and time the state point was written.

**integer(8) seed**

    Pseudo-random number generator seed.

**integer(4) run_mode**

    run mode used. The modes are described in constants.F90.

**integer(8) n_particles**

    Number of particles used per generation.

**integer(4) n_batches**

    Total number of batches (active + inactive).

**integer(4) current_batch**

    The number of batches already simulated.

if (run_mode == MODE_EIGENVALUE)

    **integer(4) n_inactive**

        Number of inactive batches

    **integer(4) gen_per_batch**

        Number of generations per batch for criticality calculations

    *do i = 1, current_batch*

        **real(8) k_batch(i)**

             k-effective for the i-th batch

    *do i = 1, current_batch*

        **real(8) entropy(i)**

            Shannon entropy for the i-th batch

**integer(4) n_meshes**

    Number of meshes in tallies.xml file

*do i = 1, n_meshes*

    **integer(4) meshes(i) % type**

        Type of mesh.

    **integer(4) meshes(i) % n_dimension**

        Number of dimensions for mesh (2 or 3).

    **integer(4) meshes(i) % dimension(:)**

        Number of mesh cells in each dimension.

    **real(8) meshes(i) % lower_left(:)**

        Coordinates of lower-left corner of mesh.

    **real(8) meshes(i) % upper_right(:)**

        Coordinates of upper-right corner of mesh.

    **real(8) meshes(i) % width(:)**

        Width of each mesh cell in each dimension.

**integer(4) n_tallies**

*do i = 1, n_tallies*

    **integer(4) size(tallies(i) % scores, 1)**

        Total number of score bins for the i-th tally

    **integer(4) size(tallies(i) % scores, 2)**

        Total number of filter bins for the i-th tally

    **integer(4) tallies(i) % n_filters**

    *do j = 1, tallies(i) % n_filters*

        **integer(4) tallies(i) % filter(j) % type**

            Type of tally filter.

        **integer(4) tallies(i) % filter(j) % n_bins**

            Number of bins for filter.

        **integer(4)/real(8) tallies(i) % filter(j) % bins(:)**

            Value for each filter bin of this type.

    **integer(4) tallies(i) % n_nuclide_bins**

        Number of nuclide bins. If none are specified, this is just one.

    *do j = 1, tallies(i) % n_nuclide_bins*

        **integer(4) tallies(i) % nuclide_bins(j)**

            Values of specified nuclide bins

    **integer(4) tallies(i) % n_score_bins**

        Number of scoring bins.

    *do j = 1, tallies(i) % n_score_bins*

        **integer(4) tallies(i) % score_bins(j)**

            Values of specified scoring bins (e.g. SCORE_FLUX).

**integer(4) N_GLOBAL_TALLIES**

    Number of global tally scores

*do i = 1, N_GLOBAL_TALLIES*

    **real(8) global_tallies(i) % sum**

        Accumulated sum for the i-th global tally

    **real(8) global_tallies(i) % sum_sq**

        Accumulated sum of squares for the i-th global tally

**integer(4) tallies_on**

    Flag indicated if tallies are present in the file.

if (tallies_on > 0)

    **integer(4) n_realizations**

        Number of realizations for tally random variables.

    *do i = 1, n_tallies*

        *do k = 1, size(tallies(i) % scores, 2)*

            *do j = 1, size(tallies(i) % scores, 1)*

                **real(8) tallies(i) % scores(j,k) % sum**
            
                    Accumulated sum for the j-th score and k-th filter of the
                    i-th tally

                **real(8) tallies(i) % scores(j,k) % sum_sq**

                    Accumulated sum of squares for the j-th score and k-th
                    filter of the i-th tally

if (run_mode == MODE_EIGENVALUE)

    *do i = 1, n_particles*

        **real(8) source_bank(i) % wgt**

            Weight of the i-th source particle

        **real(8) source_bank(i) % xyz(1:3)**

            Coordinates of the i-th source particle.

        **real(8) source_bank(i) % uvw(1:3)**

            Direction of the i-th source particle

        **real(8) source_bank(i) % E**

            Energy of the i-th source particle.

----------
Revision 3
----------

**integer(4) REVISION_STATEPOINT**

    Revision of the binary state point file. Any time a change is made in the
    format of the state-point file, this integer is incremented.

**integer(4) VERSION_MAJOR**

    Major version number for OpenMC

**integer(4) VERSION_MINOR**

    Minor version number for OpenMC

**integer(4) VERSION_RELEASE**

    Release version number for OpenMC

**character(19) time_stamp**

    Date and time the state point was written.

**integer(8) seed**

    Pseudo-random number generator seed.

**integer(4) run_mode**

    run mode used. The modes are described in constants.F90.

**integer(8) n_particles**

    Number of particles used per generation.

**integer(4) n_batches**

    Total number of batches (active + inactive).

**integer(4) current_batch**

    The number of batches already simulated.

if (run_mode == MODE_EIGENVALUE)

    **integer(4) n_inactive**

        Number of inactive batches

    **integer(4) gen_per_batch**

        Number of generations per batch for criticality calculations

    *do i = 1, current_batch*

        **real(8) k_batch(i)**

             k-effective for the i-th batch

    *do i = 1, current_batch*

        **real(8) entropy(i)**

            Shannon entropy for the i-th batch

**integer(4) N_GLOBAL_TALLIES**

    Number of global tally scores

*do i = 1, N_GLOBAL_TALLIES*

    **real(8) global_tallies(i) % sum**

        Accumulated sum for the i-th global tally

    **real(8) global_tallies(i) % sum_sq**

        Accumulated sum of squares for the i-th global tally

**integer(4) n_meshes**

    Number of meshes in tallies.xml file

*do i = 1, n_meshes*

    **integer(4) meshes(i) % type**

        Type of mesh.

    **integer(4) meshes(i) % n_dimension**

        Number of dimensions for mesh (2 or 3).

    **integer(4) meshes(i) % dimension(:)**

        Number of mesh cells in each dimension.

    **real(8) meshes(i) % lower_left(:)**

        Coordinates of lower-left corner of mesh.

    **real(8) meshes(i) % upper_right(:)**

        Coordinates of upper-right corner of mesh.

    **real(8) meshes(i) % width(:)**

        Width of each mesh cell in each dimension.

**integer(4) n_tallies**

*do i = 1, n_tallies*

    **integer(4) size(tallies(i) % scores, 1)**

        Total number of score bins for the i-th tally

    **integer(4) size(tallies(i) % scores, 2)**

        Total number of filter bins for the i-th tally

    **integer(4) tallies(i) % n_filters**

    *do j = 1, tallies(i) % n_filters*

        **integer(4) tallies(i) % filter(j) % type**

            Type of tally filter.

        **integer(4) tallies(i) % filter(j) % n_bins**

            Number of bins for filter.

        **integer(4)/real(8) tallies(i) % filter(j) % bins(:)**

            Value for each filter bin of this type.

    **integer(4) tallies(i) % n_nuclide_bins**

        Number of nuclide bins. If none are specified, this is just one.

    *do j = 1, tallies(i) % n_nuclide_bins*

        **integer(4) tallies(i) % nuclide_bins(j)**

            Values of specified nuclide bins

    **integer(4) tallies(i) % n_score_bins**

        Number of scoring bins.

    *do j = 1, tallies(i) % n_score_bins*

        **integer(4) tallies(i) % score_bins(j)**

            Values of specified scoring bins (e.g. SCORE_FLUX).

**integer(4) tallies_on**

    Flag indicated if tallies are present in the file.

if (tallies_on > 0)

    *do i = 1, n_tallies*

        *do k = 1, size(tallies(i) % scores, 2)*

            *do j = 1, size(tallies(i) % scores, 1)*

                **real(8) tallies(i) % scores(j,k) % sum**
            
                    Accumulated sum for the j-th score and k-th filter of the
                    i-th tally

                **real(8) tallies(i) % scores(j,k) % sum_sq**

                    Accumulated sum of squares for the j-th score and k-th
                    filter of the i-th tally

if (run_mode == MODE_EIGENVALUE)

    *do i = 1, n_particles*

        **real(8) source_bank(i) % wgt**

            Weight of the i-th source particle

        **real(8) source_bank(i) % xyz(1:3)**

            Coordinates of the i-th source particle.

        **real(8) source_bank(i) % uvw(1:3)**

            Direction of the i-th source particle

        **real(8) source_bank(i) % E**

            Energy of the i-th source particle.

----------
Revision 2
----------

**integer(4) REVISION_STATEPOINT**

    Revision of the binary state point file. Any time a change is made in the
    format of the state-point file, this integer is incremented.

**integer(4) VERSION_MAJOR**

    Major version number for OpenMC

**integer(4) VERSION_MINOR**

    Minor version number for OpenMC

**integer(4) VERSION_RELEASE**

    Release version number for OpenMC

**integer(4) run_mode**

    run mode used. The modes are described in constants.F90.

**integer(8) n_particles**

    Number of particles used per generation.

**integer(4) n_batches**

    Total number of batches (active + inactive).

**integer(4) n_inactive**

    Number of inactive batches

**integer(4) gen_per_batch**

    Number of generations per batch for criticality calculations

**integer(4) current_batch**

    The number of batches already simulated.

*do i = 1, current_batch*

    **real(8) k_batch(i)**

        k-effective for the i-th batch

    if (entropy_on)

        **real(8) entropy(i)**

            Shannon entropy for the i-th batch

**integer(4) N_GLOBAL_TALLIES**

    Number of global tally scores

*do i = 1, N_GLOBAL_TALLIES*

    **real(8) global_tallies(i) % sum**

        Accumulated sum for the i-th global tally

*do i = 1, N_GLOBAL_TALLIES*

    **real(8) global_tallies(i) % sum_sq**

        Accumulated sum of squares for the i-th global tally

**integer(4) n_tallies**

*do i = 1, n_tallies*

    **integer(4) size(tallies(i) % scores, 1)**

        Total number of score bins for the i-th tally

    **integer(4) size(tallies(i) % scores, 2)**

        Total number of filter bins for the i-th tally

*do i = 1, n_tallies*

    *do k = 1, size(tallies(i) % scores, 2)*

        *do j = 1, size(tallies(i) % scores, 1)*

            **real(8) tallies(i) % scores(j,k) % sum**
            
                Accumulated sum for the j-th score and k-th filter of the i-th
                tally

            **real(8) tallies(i) % scores(j,k) % sum_sq**

                Accumulated sum of squares for the j-th score and k-th filter of
                the i-th tally

----------
Revision 1
----------

**integer(4) REVISION_STATEPOINT**

    Revision of the binary state point file. Any time a change is made in the
    format of the state-point file, this integer is incremented.

**integer(4) VERSION_MAJOR**

    Major version number for OpenMC

**integer(4) VERSION_MINOR**

    Minor version number for OpenMC

**integer(4) VERSION_RELEASE**

    Release version number for OpenMC

**integer(4) run_mode**

    run mode used. The modes are described in constants.F90.

**integer(8) n_particles**

    Number of particles used per generation.

**integer(4) n_batches**

    Total number of batches (active + inactive).

**integer(4) n_inactive**

    Number of inactive batches

**integer(4) gen_per_batch**

    Number of generations per batch for criticality calculations

**integer(4) current_batch**

    The number of batches already simulated.

*do i = 1, current_batch*

    **real(8) k_batch(i)**

        k-effective for the i-th batch

    if (entropy_on)

        **real(8) entropy(i)**

            Shannon entropy for the i-th batch

**integer(4) N_GLOBAL_TALLIES**

    Number of global tally scores

*do i = 1, N_GLOBAL_TALLIES*

    **real(8) global_tallies(i) % sum**

        Accumulated sum for the i-th global tally

*do i = 1, N_GLOBAL_TALLIES*

    **real(8) global_tallies(i) % sum_sq**

        Accumulated sum of squares for the i-th global tally

**integer(4) n_tallies**

*do i = 1, n_tallies*

    **integer(4) size(tallies(i) % scores, 1)**

        Total number of score bins for the i-th tally

    **integer(4) size(tallies(i) % scores, 2)**

        Total number of filter bins for the i-th tally

    *do k = 1, size(tallies(i) % scores, 2)*

        *do j = 1, size(tallies(i) % scores, 1)*

            **real(8) tallies(i) % scores(j,k) % sum**
            
                Accumulated sum for the j-th score and k-th filter of the i-th
                tally

    *do k = 1, size(tallies(i) % scores, 2)*

        *do j = 1, size(tallies(i) % scores, 1)*

            **real(8) tallies(i) % scores(j,k) % sum_sq**

                Accumulated sum of squares for the j-th score and k-th filter of
                the i-th tally
