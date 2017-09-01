module global

  use, intrinsic :: ISO_C_BINDING

  ! Inherit module variables from other modules
  use cmfd_header
  use geometry_header
  use material_header
  use mesh_header
  use mgxs_header
  use nuclide_header
  use plot_header
  use sab_header
  use settings
  use simulation_header
  use surface_header
  use tally_filter_header
  use tally_header
  use tally_derivative_header
  use timer_header
  use trigger_header
  use volume_header

  implicit none

contains

!===============================================================================
! FREE_MEMORY deallocates and clears  all global allocatable arrays in the
! program
!===============================================================================

  subroutine free_memory()

    integer :: i ! Loop Index

    ! Deallocate cells, surfaces, materials
    if (allocated(cells)) deallocate(cells)
    if (allocated(universes)) deallocate(universes)
    if (allocated(lattices)) deallocate(lattices)
    if (allocated(surfaces)) deallocate(surfaces)
    if (allocated(materials)) deallocate(materials)
    if (allocated(plots)) deallocate(plots)
    if (allocated(volume_calcs)) deallocate(volume_calcs)

    ! Deallocate geometry debugging information
    if (allocated(overlap_check_cnt)) deallocate(overlap_check_cnt)

    ! Deallocate cross section data, listings, and cache
    if (allocated(nuclides)) then
    ! First call the clear routines
      do i = 1, size(nuclides)
        call nuclides(i) % clear()
      end do
      deallocate(nuclides)
    end if
    if (allocated(libraries)) deallocate(libraries)

    if (allocated(res_scat_nuclides)) deallocate(res_scat_nuclides)

    if (allocated(nuclides_MG)) deallocate(nuclides_MG)

    if (allocated(macro_xs)) deallocate(macro_xs)

    if (allocated(sab_tables)) deallocate(sab_tables)

    ! Deallocate external source
    if (allocated(external_source)) deallocate(external_source)

    ! Deallocate k and entropy
    if (allocated(k_generation)) deallocate(k_generation)
    if (allocated(entropy)) deallocate(entropy)
    if (allocated(entropy_p)) deallocate(entropy_p)

    ! Deallocate tally-related arrays
    if (allocated(global_tallies)) deallocate(global_tallies)
    if (allocated(meshes)) deallocate(meshes)
    if (allocated(filters)) deallocate(filters)
    if (allocated(tallies)) deallocate(tallies)

    ! Deallocate fission and source bank and entropy
!$omp parallel
    if (allocated(fission_bank)) deallocate(fission_bank)
    if (allocated(tally_derivs)) deallocate(tally_derivs)
!$omp end parallel
#ifdef _OPENMP
    if (allocated(master_fission_bank)) deallocate(master_fission_bank)
#endif
    if (allocated(source_bank)) deallocate(source_bank)

    ! Deallocate array of work indices
    if (allocated(work_index)) deallocate(work_index)

    if (allocated(energy_bins)) deallocate(energy_bins)
    if (allocated(energy_bin_avg)) deallocate(energy_bin_avg)

    ! Deallocate CMFD
    call deallocate_cmfd(cmfd)

    ! Deallocate tally node lists
    call active_analog_tallies % clear()
    call active_tracklength_tallies % clear()
    call active_current_tallies % clear()
    call active_collision_tallies % clear()
    call active_surface_tallies % clear()
    call active_tallies % clear()

    ! Deallocate track_identifiers
    if (allocated(track_identifiers)) deallocate(track_identifiers)

    ! Deallocate dictionaries
    call cell_dict % clear()
    call universe_dict % clear()
    call lattice_dict % clear()
    call surface_dict % clear()
    call material_dict % clear()
    call mesh_dict % clear()
    call filter_dict % clear()
    call tally_dict % clear()
    call plot_dict % clear()
    call nuclide_dict % clear()
    call sab_dict % clear()
    call library_dict % clear()

    ! Clear statepoint and sourcepoint batch set
    call statepoint_batch % clear()
    call sourcepoint_batch % clear()

    ! Deallocate ufs
    if (allocated(source_frac)) deallocate(source_frac)

  end subroutine free_memory

!===============================================================================
! OVERALL_GENERATION determines the overall generation number
!===============================================================================

  pure function overall_generation() result(gen)
    integer :: gen
    gen = gen_per_batch*(current_batch - 1) + current_gen
  end function overall_generation

end module global
