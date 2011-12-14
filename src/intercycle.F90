module intercycle

  use global
  use error, only: warning

#ifdef MPI
  use mpi
#endif

contains

!===============================================================================
! SHANNON_ENTROPY calculates the Shannon entropy of the fission source
! distribution to assess source convergence
!===============================================================================

  subroutine shannon_entropy()

    integer :: i_bank         ! index for bank sites
    integer :: i              ! x-index for entropy mesh
    integer :: j              ! y-index for entropy mesh
    integer :: k              ! z-index for entropy mesh
    integer, save :: n_box    ! total # of boxes on mesh
    integer, save :: n        ! # of boxes in each dimension
    real(8), save :: width(3) ! width of box in each dimension

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

    ! loop over fission sites and count how many are in each mesh box
    FISSION_SITES: do i_bank = 1, n_bank
       ! determine indices for entropy mesh box
       i = (fission_bank(i_bank) % xyz(1) - entropy_lower_left(1))/n
       j = (fission_bank(i_bank) % xyz(2) - entropy_lower_left(2))/n
       k = (fission_bank(i_bank) % xyz(3) - entropy_lower_left(3))/n

       ! if outside mesh, skip particle
       if (i < 1 .or. i > n .or. j < 1 .or. &
            j > n .or. k < 1 .or. k > n) then
          message = "Fission source site outside of entropy box."
          call warning()
          cycle
       end if

       ! add to appropriate mesh box
       entropy_p(i,j,k) = entropy_p(i,j,k) + 1
    end do FISSION_SITES

    ! normalize to number of fission sites
    entropy_p = entropy_p/n

    ! collect values from all processors
#ifdef MPI
    call MPI_REDUCE(MPI_IN_PLACE, entropy_p, n_box, MPI_REAL8, MPI_SUM, &
         0, MPI_COMM_WORLD, mpi_err)
#endif

    ! sum values to obtain shannon entropy
    if (master) entropy = sum(entropy_p * log(entropy_p)/log(2.0))

  end subroutine shannon_entropy

end module intercycle
