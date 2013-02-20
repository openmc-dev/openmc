module finalize

# ifdef PETSC
  use cmfd_output,    only: finalize_cmfd
# endif
  use global
  use output,         only: print_runtime, print_results, write_tallies
  use tally,          only: tally_statistics

#ifdef MPI
  use mpi
#endif

#ifdef HDF5
  use hdf5_interface, only: hdf5_finalize
#endif

  implicit none

contains

!===============================================================================
! FINALIZE_RUN does all post-simulation tasks such as calculating tally
! statistics and writing out tallies
!===============================================================================

  subroutine finalize_run()

    ! Start finalization timer
    call time_finalize % start()

    if (run_mode /= MODE_PLOTTING) then
      ! Calculate statistics for tallies and write to tallies.out
      if (master) call tally_statistics()
      if (master .and. run_mode == MODE_EIGENVALUE) &
           call calculate_combined_keff()
      if (output_tallies) then
        if (master) call write_tallies()
      end if
    end if

#ifdef PETSC
    ! finalize cmfd
    if (cmfd_run) call finalize_cmfd()
#endif

    ! stop timers and show timing statistics
    call time_finalize % stop()
    call time_total % stop()
    if (master .and. (run_mode /= MODE_PLOTTING .and. &
         run_mode /= MODE_TALLIES)) then
      call print_runtime()
      call print_results()
    end if

    ! deallocate arrays
    call free_memory()

#ifdef HDF5
    ! Close HDF5 interface and release memory
    call hdf5_finalize()
#endif

#ifdef MPI
    ! If MPI is in use and enabled, terminate it
    call MPI_FINALIZE(mpi_err)
#endif

  end subroutine finalize_run

!===============================================================================
! CALCULATE_COMBINED_KEFF calculates a minimum variance estimate of k-effective
! based on a linear combination of the collision, absorption, and tracklength
! estimates. The theory behind this can be found in M. Halperin, "Almost
! linearly-optimum combination of unbiased estimates," J. Am. Stat. Assoc., 56,
! 36-43 (1961), doi:10.1080/01621459.1961.10482088. The implementation here
! follows that described in T. Urbatsch et al., "Estimation and interpretation
! of keff confidence intervals in MCNP," Nucl. Technol., 111, 169-182 (1995).
!===============================================================================

  subroutine calculate_combined_keff()

    integer :: i 
    integer :: n
    real(8) :: k_col, k_abs, k_tra
    real(8) :: s_col, s_abs, s_tra
    real(8) :: s_col_abs
    real(8) :: s_col_tra
    real(8) :: s_abs_tra
    real(8) :: s_ij, s_ik, s_jk
    real(8) :: s_kk, s_jj
    real(8) :: s_1i
    real(8) :: k_i, k_j
    real(8) :: f, g
    real(8) :: S(3)

    ! Initialize variables
    n = n_realizations
    g = ZERO
    S = ZERO

    ! Copy estimates of k-effective and its variance (not variance of the mean)
    k_col = global_tallies(K_COLLISION) % sum
    k_abs = global_tallies(K_ABSORPTION) % sum
    k_tra = global_tallies(K_TRACKLENGTH) % sum
    s_col = global_tallies(K_COLLISION) % sum_sq**2 * n
    s_abs = global_tallies(K_ABSORPTION) % sum_sq**2 * n
    s_tra = global_tallies(K_TRACKLENGTH) % sum_sq**2 * n

    ! Calculate covariances based on sums with Bessel's correction
    s_col_abs = (k_col_abs - n * k_col * k_abs)/(n - 1)
    s_col_tra = (k_col_tra - n * k_col * k_tra)/(n - 1)
    s_abs_tra = (k_abs_tra - n * k_abs * k_tra)/(n - 1)

    do i = 1, 3
      if (i == 1) then
        ! i = collision, j = absorption, k = tracklength
        k_i  = k_col
        k_j  = k_abs
        s_1i = s_col
        s_jj = s_abs
        s_kk = s_tra
        s_ij = s_col_abs
        s_ik = s_col_tra
        s_jk = s_abs_tra
      elseif (i == 2) then
        ! i = absortion, j = tracklength, k = collision
        k_i  = k_abs
        k_j  = k_tra
        s_1i = s_col_abs
        s_jj = s_tra
        s_kk = s_col
        s_ij = s_abs_tra
        s_ik = s_col_abs
        s_jk = s_col_tra
      elseif (i == 3) then
        ! i = tracklength, j = collision, k = absorption
        k_i  = k_tra
        k_j  = k_col
        s_1i = s_col_tra
        s_jj = s_col
        s_kk = s_abs
        s_ij = s_col_tra
        s_ik = s_abs_tra
        s_jk = s_col_abs
      end if

      ! Calculate weighting
      f = s_jj*(s_kk - s_ik) - s_kk*s_ij + s_jk*(s_ij + s_ik - s_jk)

      ! Add to S sums for variance of combined estimate
      S(1) = S(1) + f*s_1i
      S(2) = S(2) + (s_jj + s_kk - TWO*s_jk)*k_i*k_i
      S(3) = S(3) + (s_kk + s_ij - s_jk - s_ik)*k_i*k_j

      ! Add to sum for combined k-effective
      k_combined(1) = k_combined(1) + f*k_i
      g = g + f
    end do

    ! Complete calculations of S sums
    S = (n - 1)*S
    S(1) = (n - 1)**2 * S(1)

    ! Calculate combined estimate of k-effective
    k_combined(1) = k_combined(1) / g

    ! Calculate standard deviation of combined estimate
    g = (n - 1)**2 * g
    k_combined(2) = sqrt(S(1)/(g*n*(n-3)) * (ONE + n*((S(2) - TWO*S(3))/g)))

  end subroutine calculate_combined_keff

end module finalize
