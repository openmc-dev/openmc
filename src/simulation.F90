module simulation

  use, intrinsic :: ISO_C_BINDING

  use nuclide_header,  only: micro_xs, n_nuclides
  use photon_header,   only: micro_photon_xs, n_elements
  use tally_filter_header, only: filter_matches, n_filters, filter_match_pointer

  implicit none
  private

contains

!===============================================================================
! INITIALIZE_SIMULATION
!===============================================================================

  subroutine simulation_init_f() bind(C)

    integer :: i

!$omp parallel
    ! Allocate array for microscopic cross section cache
    allocate(micro_xs(n_nuclides))
    allocate(micro_photon_xs(n_elements))

    ! Allocate array for matching filter bins
    allocate(filter_matches(n_filters))
    do i = 1, n_filters
      filter_matches(i) % ptr = filter_match_pointer(i - 1)
    end do
!$omp end parallel

  end subroutine

!===============================================================================
! FINALIZE_SIMULATION calculates tally statistics, writes tallies, and displays
! execution time and results
!===============================================================================

  subroutine simulation_finalize_f() bind(C)

    ! Free up simulation-specific memory
!$omp parallel
    deallocate(micro_xs, micro_photon_xs, filter_matches)
!$omp end parallel

  end subroutine

end module simulation
