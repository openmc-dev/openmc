module simulation

  use, intrinsic :: ISO_C_BINDING

#ifdef _OPENMP
  use omp_lib
#endif

  use bank_header,     only: source_bank
  use constants,       only: ZERO
  use error,           only: fatal_error
  use material_header, only: n_materials, materials
  use message_passing
  use nuclide_header,  only: micro_xs, n_nuclides
  use photon_header,   only: micro_photon_xs, n_elements
  use settings
  use simulation_header
  use tally_header
  use tally_filter_header, only: filter_matches, n_filters, filter_match_pointer

  implicit none
  private

contains

!===============================================================================
! INITIALIZE_SIMULATION
!===============================================================================

  subroutine simulation_init_f() bind(C)

    integer :: i

    ! Set up material nuclide index mapping
    do i = 1, n_materials
      call materials(i) % init_nuclide_index()
    end do

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

    integer :: i       ! loop index

    ! Free up simulation-specific memory
    do i = 1, n_materials
      deallocate(materials(i) % mat_nuclide_index)
    end do
!$omp parallel
    deallocate(micro_xs, micro_photon_xs, filter_matches)
!$omp end parallel

  end subroutine

!===============================================================================
! ALLOCATE_BANKS allocates memory for the fission and source banks
!===============================================================================

  subroutine allocate_banks() bind(C)

    integer :: alloc_err  ! allocation error code

    ! Allocate source bank
    if (allocated(source_bank)) deallocate(source_bank)
    allocate(source_bank(work), STAT=alloc_err)

    ! Check for allocation errors
    if (alloc_err /= 0) then
      call fatal_error("Failed to allocate source bank.")
    end if

    if (run_mode == MODE_EIGENVALUE) then

#ifdef _OPENMP
      ! If OpenMP is being used, each thread needs its own private fission
      ! bank. Since the private fission banks need to be combined at the end of
      ! a generation, there is also a 'master_fission_bank' that is used to
      ! collect the sites from each thread.

      n_threads = omp_get_max_threads()

!$omp parallel
      thread_id = omp_get_thread_num()

      if (allocated(fission_bank)) deallocate(fission_bank)
      if (thread_id == 0) then
        allocate(fission_bank(3*work))
      else
        allocate(fission_bank(3*work/n_threads))
      end if
!$omp end parallel
      if (allocated(master_fission_bank)) deallocate(master_fission_bank)
      allocate(master_fission_bank(3*work), STAT=alloc_err)
#else
      if (allocated(fission_bank)) deallocate(fission_bank)
      allocate(fission_bank(3*work), STAT=alloc_err)
#endif

      ! Check for allocation errors
      if (alloc_err /= 0) then
        call fatal_error("Failed to allocate fission bank.")
      end if
    end if

  end subroutine allocate_banks

end module simulation
