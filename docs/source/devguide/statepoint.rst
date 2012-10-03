.. _devguide_statepoint:

======================================
State Point Binary File Specifications
======================================

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

        **real(8) entropy for the i-th batch**

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

        **real(8) entropy for the i-th batch**

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
