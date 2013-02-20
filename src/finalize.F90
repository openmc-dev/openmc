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

    integer :: l        ! loop index
    integer :: i, j, k  ! indices referring to collision, absorption, or track
    integer :: n        ! number of realizations
    real(8) :: kv(3)    ! vector of k-effective estimates
    real(8) :: cov(3,3) ! sample covariance matrix
    real(8) :: f        ! weighting factor
    real(8) :: g        ! sum of weighting factors
    real(8) :: S(3)     ! sums used for variance calculation

    ! Initialize variables
    n = n_realizations
    g = ZERO
    S = ZERO

    ! Copy estimates of k-effective and its variance (not variance of the mean)
    kv(1) = global_tallies(K_COLLISION) % sum
    kv(2) = global_tallies(K_ABSORPTION) % sum
    kv(3) = global_tallies(K_TRACKLENGTH) % sum
    cov(1,1) = global_tallies(K_COLLISION) % sum_sq**2 * n
    cov(2,2) = global_tallies(K_ABSORPTION) % sum_sq**2 * n
    cov(3,3) = global_tallies(K_TRACKLENGTH) % sum_sq**2 * n

    ! Calculate covariances based on sums with Bessel's correction
    cov(1,2) = (k_col_abs - n * kv(1) * kv(2))/(n - 1)
    cov(1,3) = (k_col_tra - n * kv(1) * kv(3))/(n - 1)
    cov(2,3) = (k_abs_tra - n * kv(2) * kv(3))/(n - 1)
    cov(2,1) = cov(1,2)
    cov(3,1) = cov(1,3)
    cov(3,2) = cov(2,3)

    do l = 1, 3
      ! Permutations of estimates
      if (l == 1) then
        ! i = collision, j = absorption, k = tracklength
        i = 1
        j = 2
        k = 3
      elseif (l == 2) then
        ! i = absortion, j = tracklength, k = collision
        i = 2
        j = 3
        k = 1
      elseif (l == 3) then
        ! i = tracklength, j = collision, k = absorption
        i = 3
        j = 1
        k = 2
      end if

      ! Calculate weighting
      f = cov(j,j)*(cov(k,k) - cov(i,k)) - cov(k,k)*cov(i,j) + &
           cov(j,k)*(cov(i,j) + cov(i,k) - cov(j,k))

      ! Add to S sums for variance of combined estimate
      S(1) = S(1) + f * cov(1,l)
      S(2) = S(2) + (cov(j,j) + cov(k,k) - TWO*cov(j,k))*kv(l)*kv(l)
      S(3) = S(3) + (cov(k,k) + cov(i,j) - cov(j,k) - cov(i,k))*kv(l)*kv(j)

      ! Add to sum for combined k-effective
      k_combined(1) = k_combined(1) + f * kv(l)
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
