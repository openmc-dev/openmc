module openmc_api

  use, intrinsic :: ISO_C_BINDING

  use hdf5, only: HID_T, h5tclose_f, h5close_f

  use constants,       only: K_BOLTZMANN
  use eigenvalue,      only: k_sum, openmc_get_keff
  use error
  use geometry,        only: find_cell
  use geometry_header
  use hdf5_interface
  use material_header
  use mesh_header
  use message_passing
  use nuclide_header
  use initialize,      only: openmc_init
  use particle_header, only: Particle
  use plot,            only: openmc_plot_geometry
  use random_lcg,      only: seed, openmc_set_seed
  use settings
  use simulation_header
  use tally_header
  use tally_filter_header
  use tally_filter
  use tally,           only: openmc_tally_set_type
  use simulation,      only: openmc_run
  use string,          only: to_f_string
  use timer_header
  use volume_calc,     only: openmc_calculate_volumes

  implicit none

  private
  public :: openmc_calculate_volumes
  public :: openmc_cell_get_id
  public :: openmc_cell_get_fill
  public :: openmc_cell_set_fill
  public :: openmc_cell_set_id
  public :: openmc_cell_set_temperature
  public :: openmc_energy_filter_get_bins
  public :: openmc_energy_filter_set_bins
  public :: openmc_extend_filters
  public :: openmc_extend_cells
  public :: openmc_extend_materials
  public :: openmc_extend_tallies
  public :: openmc_filter_get_id
  public :: openmc_filter_get_type
  public :: openmc_filter_set_id
  public :: openmc_filter_set_type
  public :: openmc_finalize
  public :: openmc_find
  public :: openmc_get_cell_index
  public :: openmc_get_keff
  public :: openmc_get_filter_index
  public :: openmc_get_filter_next_id
  public :: openmc_get_material_index
  public :: openmc_get_nuclide_index
  public :: openmc_get_tally_index
  public :: openmc_hard_reset
  public :: openmc_init
  public :: openmc_load_nuclide
  public :: openmc_material_add_nuclide
  public :: openmc_material_get_id
  public :: openmc_material_get_densities
  public :: openmc_material_set_density
  public :: openmc_material_set_densities
  public :: openmc_material_set_id
  public :: openmc_material_filter_get_bins
  public :: openmc_material_filter_set_bins
  public :: openmc_mesh_filter_set_mesh
  public :: openmc_nuclide_name
  public :: openmc_plot_geometry
  public :: openmc_reset
  public :: openmc_run
  public :: openmc_tally_get_id
  public :: openmc_tally_get_filters
  public :: openmc_tally_get_nuclides
  public :: openmc_tally_get_scores
  public :: openmc_tally_results
  public :: openmc_tally_set_filters
  public :: openmc_tally_set_id
  public :: openmc_tally_set_nuclides
  public :: openmc_tally_set_scores
  public :: openmc_tally_set_type

contains

!===============================================================================
! OPENMC_FINALIZE frees up memory by deallocating arrays and resetting global
! variables
!===============================================================================

  subroutine openmc_finalize() bind(C)

    integer :: err

    ! Clear results
    call openmc_reset()

    ! Reset global variables
    assume_separate = .false.
    check_overlaps = .false.
    confidence_intervals = .false.
    create_fission_neutrons = .true.
    energy_cutoff = ZERO
    energy_max_neutron = INFINITY
    energy_min_neutron = ZERO
    entropy_on = .false.
    gen_per_batch = 1
    keff = ONE
    legendre_to_tabular = .true.
    legendre_to_tabular_points = 33
    n_batch_interval = 1
    n_lost_particles = 0
    n_particles = 0
    n_source_points = 0
    n_state_points = 0
    n_tallies = 0
    output_summary = .true.
    output_tallies = .true.
    particle_restart_run = .false.
    pred_batches = .false.
    reduce_tallies = .true.
    res_scat_on = .false.
    res_scat_method = RES_SCAT_ARES
    res_scat_energy_min = 0.01_8
    res_scat_energy_max = 1000.0_8
    restart_run = .false.
    root_universe = -1
    run_CE = .true.
    run_mode = NONE
    satisfy_triggers = .false.
    seed = 1_8
    source_latest = .false.
    source_separate = .false.
    source_write = .true.
    survival_biasing = .false.
    temperature_default = 293.6_8
    temperature_method = TEMPERATURE_NEAREST
    temperature_multipole = .false.
    temperature_range = [ZERO, ZERO]
    temperature_tolerance = 10.0_8
    total_gen = 0
    trigger_on = .false.
    ufs = .false.
    urr_ptables_on = .true.
    verbosity = 7
    weight_cutoff = 0.25_8
    weight_survive = ONE
    write_all_tracks = .false.
    write_initial_source = .false.

    ! Deallocate arrays
    call free_memory()

    ! Release compound datatypes
    call h5tclose_f(hdf5_bank_t, err)

    ! Close FORTRAN interface.
    call h5close_f(err)

#ifdef MPI
    ! Free all MPI types
    call MPI_TYPE_FREE(MPI_BANK, err)
#endif

  end subroutine openmc_finalize

!===============================================================================
! OPENMC_FIND determines the ID or a cell or material at a given point in space
!===============================================================================

  function openmc_find(xyz, rtype, id, instance) result(err) bind(C)
    real(C_DOUBLE), intent(in)        :: xyz(3) ! Cartesian point
    integer(C_INT), intent(in), value :: rtype  ! 1 for cell, 2 for material
    integer(C_INT32_T), intent(out)   :: id
    integer(C_INT32_T), intent(out)   :: instance
    integer(C_INT) :: err

    logical :: found
    type(Particle) :: p

    call p % initialize()
    p % coord(1) % xyz(:) = xyz
    p % coord(1) % uvw(:) = [ZERO, ZERO, ONE]
    call find_cell(p, found)

    id = -1
    instance = -1
    err = E_UNASSIGNED

    if (found) then
      if (rtype == 1) then
        id = cells(p % coord(p % n_coord) % cell) % id
      elseif (rtype == 2) then
        if (p % material == MATERIAL_VOID) then
          id = 0
        else
          id = materials(p % material) % id
        end if
      end if
      instance = p % cell_instance - 1
      err = 0
    else
      err = E_GEOMETRY
      call set_errmsg("Could not find cell/material at position (" // &
           trim(to_str(xyz(1))) // "," // trim(to_str(xyz(2))) // "," // &
           trim(to_str(xyz(3))) // ").")
    end if

  end function openmc_find

!===============================================================================
! OPENMC_HARD_RESET reset tallies and timers as well as the pseudorandom
! generator state
!===============================================================================

  subroutine openmc_hard_reset() bind(C)
    integer :: err

    ! Reset all tallies and timers
    call openmc_reset()

    ! Reset total generations and keff guess
    keff = ONE
    total_gen = 0

    ! Reset the random number generator state
    err = openmc_set_seed(1_8)
  end subroutine openmc_hard_reset

!===============================================================================
! OPENMC_RESET resets tallies and timers
!===============================================================================

  subroutine openmc_reset() bind(C)
    integer :: i

    if (allocated(tallies)) then
      do i = 1, size(tallies)
        associate (t => tallies(i) % obj)
          t % active = .false.
          t % n_realizations = 0
          if (allocated(t % results)) then
            t % results(:, :, :) = ZERO
          end if
        end associate
      end do
    end if

    ! Reset global tallies
    n_realizations = 0
    if (allocated(global_tallies)) then
      global_tallies(:, :) = ZERO
    end if
    k_col_abs = ZERO
    k_col_tra = ZERO
    k_abs_tra = ZERO
    k_sum(:) = ZERO

    ! Clear active tally lists
    call active_analog_tallies % clear()
    call active_tracklength_tallies % clear()
    call active_current_tallies % clear()
    call active_collision_tallies % clear()
    call active_tallies % clear()

    ! Reset timers
    call time_total % reset()
    call time_total % reset()
    call time_initialize % reset()
    call time_read_xs % reset()
    call time_unionize % reset()
    call time_bank % reset()
    call time_bank_sample % reset()
    call time_bank_sendrecv % reset()
    call time_tallies % reset()
    call time_inactive % reset()
    call time_active % reset()
    call time_transport % reset()
    call time_finalize % reset()

  end subroutine openmc_reset

!===============================================================================
! FREE_MEMORY deallocates and clears  all global allocatable arrays in the
! program
!===============================================================================

  subroutine free_memory()

    use cmfd_header
    use mgxs_header
    use plot_header
    use sab_header
    use settings
    use source_header
    use surface_header
    use tally_derivative_header
    use trigger_header
    use volume_header

    call free_memory_geometry()
    call free_memory_surfaces()
    call free_memory_material()
    call free_memory_plot()
    call free_memory_volume()
    call free_memory_simulation()
    call free_memory_nuclide()
    call free_memory_settings()
    call free_memory_mgxs()
    call free_memory_sab()
    call free_memory_source()
    call free_memory_mesh()
    call free_memory_tally()
    call free_memory_tally_filter()
    call free_memory_tally_derivative()
    call free_memory_bank()

    ! Deallocate CMFD
    call deallocate_cmfd(cmfd)

  end subroutine free_memory

end module openmc_api
