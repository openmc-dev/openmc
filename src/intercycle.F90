module intercycle

  use ISO_FORTRAN_ENV

  use global
  use error,  only: warning
  use output, only: write_message

#ifdef MPI
  use mpi
#endif

contains

!===============================================================================
! SHANNON_ENTROPY calculates the Shannon entropy of the fission source
! distribution to assess source convergence
!===============================================================================

  subroutine shannon_entropy()

    integer :: i              ! x-index for entropy mesh
    integer :: j              ! y-index for entropy mesh
    integer :: k              ! z-index for entropy mesh
    integer :: m              ! index for bank sites
    integer(8) :: total_bank  ! total # of fission bank sites
    integer, save :: n_box    ! total # of boxes on mesh
    integer, save :: n        ! # of boxes in each dimension
    real(8), save :: width(3) ! width of box in each dimension
    logical :: outside_box    ! were there sites outside entropy box?

    ! On the first pass through this subroutine, we need to determine how big
    ! the entropy mesh should be in each direction and then allocate a
    ! three-dimensional array to store the fraction of source sites in each mesh
    ! box

    if (.not. allocated(entropy_p)) then
       ! determine number of boxes in each direction
       n = ceiling((n_particles/20)**(1.0/3.0))
       n_box = n*n*n

       ! determine width
       width = (entropy_upper_right - entropy_lower_left)/n

       ! allocate p
       allocate(entropy_p(n,n,n))
    end if

    ! initialize p
    entropy_p = ZERO
    outside_box = .false.

    ! loop over fission sites and count how many are in each mesh box
    FISSION_SITES: do m = 1, int(n_bank,4)
       ! determine indices for entropy mesh box
       i = int((fission_bank(m) % xyz(1) - entropy_lower_left(1))/width(1)) + 1
       j = int((fission_bank(m) % xyz(2) - entropy_lower_left(2))/width(2)) + 1 
       k = int((fission_bank(m) % xyz(3) - entropy_lower_left(3))/width(3)) + 1

       ! if outside mesh, skip particle
       if (i < 1 .or. i > n .or. j < 1 .or. &
            j > n .or. k < 1 .or. k > n) then
          outside_box = .true.
          cycle
       end if

       ! add to appropriate mesh box
       entropy_p(i,j,k) = entropy_p(i,j,k) + 1
    end do FISSION_SITES

    ! display warning message if there were sites outside entropy box
    if (outside_box) then
       message = "Fission source site(s) outside of entropy box."
       call warning()
    end if

#ifdef MPI
    ! collect values from all processors
    if (master) then
       call MPI_REDUCE(MPI_IN_PLACE, entropy_p, n_box, MPI_REAL8, MPI_SUM, &
            0, MPI_COMM_WORLD, mpi_err)
    else
       call MPI_REDUCE(entropy_p, entropy_p, n_box, MPI_REAL8, MPI_SUM, &
            0, MPI_COMM_WORLD, mpi_err)
    end if

    ! determine total number of bank sites
    call MPI_REDUCE(n_bank, total_bank, 1, MPI_INTEGER8, MPI_SUM, 0, &
         MPI_COMM_WORLD, mpi_err)
#else
    total_bank = n_bank
#endif

    ! sum values to obtain shannon entropy
    if (master) then
       entropy_p = entropy_p / total_bank
       entropy = -sum(entropy_p * log(entropy_p)/log(2.0), entropy_p > ZERO)
    end if

  end subroutine shannon_entropy

!===============================================================================
! CALCULATE_KEFF calculates the single cycle estimate of keff as well as the
! mean and standard deviation of the mean for active cycles and displays them
!===============================================================================

  subroutine calculate_keff(i_cycle)

    integer, intent(in) :: i_cycle ! index of current cycle

    integer(8)    :: total_bank ! total number of source sites
    integer       :: n          ! active cycle number
    real(8)       :: k_cycle    ! single cycle estimate of keff
    real(8), save :: k_sum      ! accumulated keff
    real(8), save :: k_sum_sq   ! accumulated keff**2

    message = "Calculate cycle keff..."
    call write_message(8)

    ! initialize sum and square of sum at beginning of run
    if (i_cycle == 1) then
       k_sum = ZERO
       k_sum_sq = ZERO
    end if

#ifdef MPI
    ! Collect number bank sites onto master process
    call MPI_REDUCE(n_bank, total_bank, 1, MPI_INTEGER8, MPI_SUM, 0, &
         & MPI_COMM_WORLD, mpi_err)
#else
    total_bank = n_bank
#endif

    ! Collect statistics and print output
    if (master) then
       ! Since the creation of bank sites was originally weighted by the last
       ! cycle keff, we need to multiply by that keff to get the current cycle's
       ! value

       k_cycle = real(total_bank)/real(n_particles)*keff

       if (i_cycle > n_inactive) then
          ! Active cycle number
          n = i_cycle - n_inactive

          ! Accumulate cycle estimate of k
          k_sum =    k_sum    + k_cycle
          k_sum_sq = k_sum_sq + k_cycle*k_cycle

          ! Determine mean and standard deviation of mean
          keff = k_sum/n
          keff_std = sqrt((k_sum_sq/n - keff*keff)/n)

          ! Display output for this cycle
          if (i_cycle > n_inactive+1) then
             if (entropy_on) then
                write(UNIT=OUTPUT_UNIT, FMT=103) i_cycle, k_cycle, entropy, &
                     keff, keff_std
             else
                write(UNIT=OUTPUT_UNIT, FMT=101) i_cycle, k_cycle, keff, keff_std
             end if
          else
             if (entropy_on) then
                write(UNIT=OUTPUT_UNIT, FMT=102) i_cycle, k_cycle, entropy
             else
                write(UNIT=OUTPUT_UNIT, FMT=100) i_cycle, k_cycle
             end if
          end if
       else
          ! Display output for inactive cycle
          if (entropy_on) then
             write(UNIT=OUTPUT_UNIT, FMT=102) i_cycle, k_cycle, entropy
          else
             write(UNIT=OUTPUT_UNIT, FMT=100) i_cycle, k_cycle
          end if
          keff = k_cycle
       end if
    end if

#ifdef MPI
    ! Broadcast new keff value to all processors
    call MPI_BCAST(keff, 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
#endif

100 format (2X,I5,2X,F8.5)
101 format (2X,I5,2X,F8.5,5X,F8.5," +/-",F8.5)
102 format (2X,I5,2X,F8.5,3X,F8.5)
103 format (2X,I5,2X,F8.5,3X,F8.5,3X,F8.5," +/-",F8.5)

  end subroutine calculate_keff

end module intercycle
