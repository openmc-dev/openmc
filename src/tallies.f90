module tallies

  use global
  use mpi
  use output, only: message, error

  implicit none

contains

!=====================================================================
! CALCULATE_KEFF
!=====================================================================

  subroutine calculate_keff(i_cycle)

    integer, intent(in) :: i_cycle ! index of current cycle

    integer :: total_bank
    integer :: n
    integer :: ierr
    real(8) :: kcoll   ! keff collision estimator         
    real(8) :: ktemp   ! MPI-reduced keff and stdev
    real(8), save :: k1 = 0. ! accumulated keff
    real(8), save :: k2 = 0. ! accumulated keff**2
    real(8) :: std     ! stdev of keff over active cycles
    character(250) :: msg

    msg = "Calculate cycle keff..."
    call message(msg, 8)

    ! set k1 and k2 at beginning of run
    if (i_cycle == 1) then
       k1 = ZERO
       k2 = ZERO
    end if

    ! Collect number bank sites onto master process
    call MPI_REDUCE(n_bank, total_bank, 1, MPI_INTEGER, MPI_SUM, 0, &
         & MPI_COMM_WORLD, ierr)

    ! Collect statistics and print output
    if (master) then
       kcoll = real(total_bank)/real(n_particles)*keff
       if (i_cycle > n_inactive) then
          n = i_cycle - n_inactive
          k1 = k1 + kcoll
          k2 = k2 + kcoll**2
          keff = k1/n
          std  = sqrt((k2/n-keff**2)/n)
          if (i_cycle > n_inactive+1) then
             write(6,101) i_cycle, kcoll, keff, std
          else
             write(6,100) i_cycle, kcoll
          end if
       else
          write(6,100) i_cycle, kcoll
          keff = kcoll
       end if
    end if
    call MPI_BCAST(keff, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

100 format (2X,I4,3X,F7.5)
101 format (2X,I4,3X,F7.5,10X,F7.5,2X,F7.5)

  end subroutine calculate_keff

end module tallies
