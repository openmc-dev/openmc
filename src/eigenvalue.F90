module eigenvalue

  use, intrinsic :: ISO_C_BINDING

  use message_passing
  use settings
  use simulation_header
  use tally_header

  implicit none

  interface
    subroutine calculate_generation_keff() bind(C)
    end subroutine

    subroutine calculate_average_keff() bind(C)
    end subroutine

    function openmc_get_keff(k_combined) result(err) bind(C)
      import C_INT, C_DOUBLE
      real(C_DOUBLE), intent(out) :: k_combined(2)
      integer(C_INT) :: err
    end function
  end interface

contains

#ifdef _OPENMP
!===============================================================================
! JOIN_BANK_FROM_THREADS joins threadprivate fission banks into a single fission
! bank that can be sampled. Note that this operation is necessarily sequential
! to preserve the order of the bank when using varying numbers of threads.
!===============================================================================

  subroutine join_bank_from_threads()

    integer(8) :: total ! total number of fission bank sites
    integer    :: i     ! loop index for threads

    ! Initialize the total number of fission bank sites
    total = 0

!$omp parallel

    ! Copy thread fission bank sites to one shared copy
!$omp do ordered schedule(static)
    do i = 1, n_threads
!$omp ordered
      master_fission_bank(total+1:total+n_bank) = fission_bank(1:n_bank)
      total = total + n_bank
!$omp end ordered
    end do
!$omp end do

    ! Make sure all threads have made it to this point
!$omp barrier

    ! Now copy the shared fission bank sites back to the master thread's copy.
    if (thread_id == 0) then
      n_bank = total
      fission_bank(1:n_bank) = master_fission_bank(1:n_bank)
    else
      n_bank = 0
    end if

!$omp end parallel

  end subroutine join_bank_from_threads
#endif

end module eigenvalue
