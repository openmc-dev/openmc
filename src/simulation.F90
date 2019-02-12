module simulation

  use, intrinsic :: ISO_C_BINDING

  use nuclide_header,  only: micro_xs, n_nuclides
  use photon_header,   only: micro_photon_xs, n_elements

  implicit none
  private

contains

!===============================================================================
! INITIALIZE_SIMULATION
!===============================================================================

  subroutine simulation_init_f() bind(C)

!$omp parallel
    ! Allocate array for microscopic cross section cache
    allocate(micro_xs(n_nuclides))
    allocate(micro_photon_xs(n_elements))
!$omp end parallel

  end subroutine

!===============================================================================
! FINALIZE_SIMULATION calculates tally statistics, writes tallies, and displays
! execution time and results
!===============================================================================

  subroutine simulation_finalize_f() bind(C)
    ! Free up simulation-specific memory
!$omp parallel
    deallocate(micro_xs, micro_photon_xs)
!$omp end parallel
  end subroutine

end module simulation
