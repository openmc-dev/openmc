module openmc_api

  use, intrinsic :: ISO_C_BINDING

  use hdf5, only: HID_T, h5tclose_f, h5close_f

  use constants,       only: K_BOLTZMANN
  use eigenvalue,      only: k_sum, openmc_get_keff
  use error
  use geometry,        only: find_cell
  use geometry_header, only: root_universe
  use global
  use hdf5_interface
  use message_passing
  use initialize,      only: openmc_init
  use input_xml,       only: assign_0K_elastic_scattering, check_data_version
  use particle_header, only: Particle
  use plot,            only: openmc_plot_geometry
  use random_lcg,      only: seed, initialize_prng
  use tally_header
  use tally_filter_header
  use tally_filter
  use tally,           only: openmc_tally_set_type
  use simulation,      only: openmc_run
  use string,          only: to_f_string
  use volume_calc,     only: openmc_calculate_volumes

  implicit none

  private
  public :: openmc_calculate_volumes
  public :: openmc_cell_get_id
  public :: openmc_cell_set_temperature
  public :: openmc_energy_filter_get_bins
  public :: openmc_energy_filter_set_bins
  public :: openmc_extend_filters
  public :: openmc_extend_tallies
  public :: openmc_filter_get_id
  public :: openmc_filter_get_type
  public :: openmc_filter_set_id
  public :: openmc_filter_set_type
  public :: openmc_finalize
  public :: openmc_find
  public :: openmc_get_cell
  public :: openmc_get_keff
  public :: openmc_get_material
  public :: openmc_get_nuclide
  public :: openmc_get_tally
  public :: openmc_hard_reset
  public :: openmc_init
  public :: openmc_load_nuclide
  public :: openmc_material_add_nuclide
  public :: openmc_material_get_id
  public :: openmc_material_get_densities
  public :: openmc_material_set_density
  public :: openmc_material_set_densities
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
  public :: openmc_tally_results
  public :: openmc_tally_set_filters
  public :: openmc_tally_set_id
  public :: openmc_tally_set_nuclides
  public :: openmc_tally_set_scores
  public :: openmc_tally_set_type

contains

!===============================================================================
! OPENMC_CELL_GET_ID returns the ID of a cell
!===============================================================================

  function openmc_cell_get_id(index, id) result(err) bind(C)
    integer(C_INT32_T), value       :: index
    integer(C_INT32_T), intent(out) :: id
    integer(C_INT) :: err

    if (index >= 1 .and. index <= size(cells)) then
      id = cells(index) % id
      err = 0
    else
      err = E_OUT_OF_BOUNDS
    end if
  end function openmc_cell_get_id

!===============================================================================
! OPENMC_CELL_SET_TEMPERATURE sets the temperature of a cell
!===============================================================================

  function openmc_cell_set_temperature(index, T, instance) result(err) bind(C)
    integer(C_INT32_T), value, intent(in)    :: index    ! cell index in cells
    real(C_DOUBLE), value, intent(in)        :: T        ! temperature
    integer(C_INT32_T), optional, intent(in) :: instance ! cell instance

    integer(C_INT) :: err     ! error code
    integer :: j              ! looping variable
    integer :: n              ! number of cell instances
    integer :: material_index ! material index in materials array
    integer :: num_nuclides   ! num nuclides in material
    integer :: nuclide_index  ! index of nuclide in nuclides array
    real(8) :: min_temp       ! min common-denominator avail temp
    real(8) :: max_temp       ! max common-denominator avail temp
    real(8) :: temp           ! actual temp we'll assign
    logical :: outside_low    ! lower than available data
    logical :: outside_high   ! higher than available data

    outside_low = .false.
    outside_high = .false.

    err = E_UNASSIGNED

    if (index >= 1 .and. index <= size(cells)) then

      ! error if the cell is filled with another universe
      if (cells(index) % fill /= NONE) then
        err = E_CELL_NO_MATERIAL
      else
        ! find which material is associated with this cell (material_index
        ! is the index into the materials array)
        if (present(instance)) then
          material_index = cells(index) % material(instance)
        else
          material_index = cells(index) % material(1)
        end if

        ! number of nuclides associated with this material
        num_nuclides = size(materials(material_index) % nuclide)

        min_temp = ZERO
        max_temp = INFINITY

        do j = 1, num_nuclides
          nuclide_index = materials(material_index) % nuclide(j)
          min_temp = max(min_temp, minval(nuclides(nuclide_index) % kTs))
          max_temp = min(max_temp, maxval(nuclides(nuclide_index) % kTs))
        end do

        ! adjust the temperature to be within bounds if necessary
        if (K_BOLTZMANN * T < min_temp) then
          outside_low = .true.
          temp = min_temp / K_BOLTZMANN
        else if (K_BOLTZMANN * T > max_temp) then
          outside_high = .true.
          temp = max_temp / K_BOLTZMANN
        else
          temp = T
        end if

        associate (c => cells(index))
          if (allocated(c % sqrtkT)) then
            n = size(c % sqrtkT)
            if (present(instance) .and. n > 1) then
              if (instance >= 0 .and. instance < n) then
                c % sqrtkT(instance + 1) = sqrt(K_BOLTZMANN * temp)
                err = 0
              end if
            else
              c % sqrtkT(:) = sqrt(K_BOLTZMANN * temp)
              err = 0
            end if
          end if
        end associate

        ! Assign error codes for outside of temperature bounds provided the
        ! temperature was changed correctly. This needs to be done after
        ! changing the temperature based on the logical structure above.
        if (err == 0) then
          if (outside_low) err = W_BELOW_MIN_BOUND
          if (outside_high) err = W_ABOVE_MAX_BOUND
        end if

      end if

    else
      err = E_OUT_OF_BOUNDS
    end if
  end function openmc_cell_set_temperature

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
    n_filters = 0
    n_meshes = 0
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
      err = E_CELL_NOT_FOUND
    end if

  end function openmc_find

!===============================================================================
! OPENMC_GET_CELL returns the index in the cells array of a cell with a given ID
!===============================================================================

  function openmc_get_cell(id, index) result(err) bind(C)
    integer(C_INT32_T), value :: id
    integer(C_INT32_T), intent(out) :: index
    integer(C_INT) :: err

    if (allocated(cells)) then
      if (cell_dict % has_key(id)) then
        index = cell_dict % get_key(id)
        err = 0
      else
        err = E_CELL_INVALID_ID
      end if
    else
      err = E_CELL_NOT_ALLOCATED
    end if
  end function openmc_get_cell

!===============================================================================
! OPENMC_GET_MATERIAL returns the index in the materials array of a material
! with a given ID
!===============================================================================

  function openmc_get_material(id, index) result(err) bind(C)
    integer(C_INT32_T), value :: id
    integer(C_INT32_T), intent(out) :: index
    integer(C_INT) :: err

    if (allocated(materials)) then
      if (material_dict % has_key(id)) then
        index = material_dict % get_key(id)
        err = 0
      else
        err = E_MATERIAL_INVALID_ID
      end if
    else
      err = E_MATERIAL_NOT_ALLOCATED
    end if
  end function openmc_get_material

!===============================================================================
! OPENMC_GET_NUCLIDE returns the index in the nuclides array of a nuclide
! with a given name
!===============================================================================

  function openmc_get_nuclide(name, index) result(err) bind(C)
    character(kind=C_CHAR), intent(in) :: name(*)
    integer(C_INT), intent(out) :: index
    integer(C_INT) :: err

    character(:), allocatable :: name_

    ! Copy array of C_CHARs to normal Fortran string
    name_ = to_f_string(name)

    if (allocated(nuclides)) then
      if (nuclide_dict % has_key(to_lower(name_))) then
        index = nuclide_dict % get_key(to_lower(name_))
        err = 0
      else
        err = E_NUCLIDE_NOT_LOADED
      end if
    else
      err = E_NUCLIDE_NOT_ALLOCATED
    end if
  end function openmc_get_nuclide

!===============================================================================
! OPENMC_HARD_RESET reset tallies and timers as well as the pseudorandom
! generator state
!===============================================================================

  subroutine openmc_hard_reset() bind(C)
    ! Reset all tallies and timers
    call openmc_reset()

    ! Reset total generations and keff guess
    keff = ONE
    total_gen = 0

    ! Reset the random number generator state
    seed = 1_8
    call initialize_prng()
  end subroutine openmc_hard_reset

!===============================================================================
! OPENMC_LOAD_NUCLIDE loads a nuclide from the cross section library
!===============================================================================

  function openmc_load_nuclide(name) result(err) bind(C)
    character(kind=C_CHAR), intent(in) :: name(*)
    integer(C_INT) :: err

    integer :: i_library
    integer :: n
    integer(HID_T) :: file_id
    integer(HID_T) :: group_id
    character(:), allocatable :: name_
    real(8) :: minmax(2) = [ZERO, INFINITY]
    type(VectorReal) :: temperature
    type(Nuclide), allocatable :: new_nuclides(:)

    ! Copy array of C_CHARs to normal Fortran string
    name_ = to_f_string(name)

    err = 0
    if (.not. nuclide_dict % has_key(to_lower(name_))) then
      if (library_dict % has_key(to_lower(name_))) then
        ! allocate extra space in nuclides array
        n = n_nuclides
        allocate(new_nuclides(n + 1))
        new_nuclides(1:n) = nuclides(:)
        call move_alloc(FROM=new_nuclides, TO=nuclides)
        n = n + 1

        i_library = library_dict % get_key(to_lower(name_))

        ! Open file and make sure version is sufficient
        file_id = file_open(libraries(i_library) % path, 'r')
        call check_data_version(file_id)

        ! Read nuclide data from HDF5
        group_id = open_group(file_id, name_)
        call nuclides(n) % from_hdf5(group_id, temperature, &
             temperature_method, temperature_tolerance, minmax, &
             master)
        call close_group(group_id)
        call file_close(file_id)

        ! Add entry to nuclide dictionary
        call nuclide_dict % add_key(to_lower(name_), n)
        n_nuclides = n

        ! Assign resonant scattering data
        if (res_scat_on) call assign_0K_elastic_scattering(nuclides(n))

        ! Initialize nuclide grid
        call nuclides(n) % init_grid(energy_min_neutron, &
             energy_max_neutron, n_log_bins)
      else
        err = E_NUCLIDE_NOT_IN_LIBRARY
      end if
    end if

  end function openmc_load_nuclide

!===============================================================================
! OPENMC_MATERIAL_ADD_NUCLIDE
!===============================================================================

  function openmc_material_add_nuclide(index, name, density) result(err) bind(C)
    integer(C_INT32_T), value, intent(in) :: index
    character(kind=C_CHAR) :: name(*)
    real(C_DOUBLE), value, intent(in) :: density
    integer(C_INT) :: err

    integer :: j, k, n
    real(8) :: awr
    integer, allocatable :: new_nuclide(:)
    real(8), allocatable :: new_density(:)
    character(:), allocatable :: name_

    name_ = to_f_string(name)

    err = E_UNASSIGNED
    if (index >= 1 .and. index <= size(materials)) then
      associate (m => materials(index))
        ! Check if nuclide is already in material
        do j = 1, size(m % nuclide)
          k = m % nuclide(j)
          if (nuclides(k) % name == name_) then
            awr = nuclides(k) % awr
            m % density = m % density + density - m % atom_density(j)
            m % density_gpcc = m % density_gpcc + (density - &
                 m % atom_density(j)) * awr * MASS_NEUTRON / N_AVOGADRO
            m % atom_density(j) = density
            err = 0
          end if
        end do

        ! If nuclide wasn't found, extend nuclide/density arrays
        if (err /= 0) then
          ! If nuclide hasn't been loaded, load it now
          err = openmc_load_nuclide(name)

          if (err == 0) then
            ! Extend arrays
            n = size(m % nuclide)
            allocate(new_nuclide(n + 1))
            new_nuclide(1:n) = m % nuclide
            call move_alloc(FROM=new_nuclide, TO=m % nuclide)

            allocate(new_density(n + 1))
            new_density(1:n) = m % atom_density
            call move_alloc(FROM=new_density, TO=m % atom_density)

            ! Append new nuclide/density
            k = nuclide_dict % get_key(to_lower(name_))
            m % nuclide(n + 1) = k
            m % atom_density(n + 1) = density
            m % density = m % density + density
            m % density_gpcc = m % density_gpcc + &
                 density * nuclides(k) % awr * MASS_NEUTRON / N_AVOGADRO
            m % n_nuclides = n + 1
          end if
        end if
      end associate
    else
      err = E_OUT_OF_BOUNDS
    end if

  end function openmc_material_add_nuclide

!===============================================================================
! OPENMC_MATERIAL_GET_DENSITIES returns an array of nuclide densities in a
! material
!===============================================================================

  function openmc_material_get_densities(index, nuclides, densities, n) &
       result(err) bind(C)
    integer(C_INT32_T), value :: index
    type(C_PTR),        intent(out) :: nuclides
    type(C_PTR),        intent(out) :: densities
    integer(C_INT),     intent(out) :: n
    integer(C_INT) :: err

    err = E_UNASSIGNED
    if (index >= 1 .and. index <= size(materials)) then
      associate (m => materials(index))
        if (allocated(m % atom_density)) then
          nuclides = C_LOC(m % nuclide(1))
          densities = C_LOC(m % atom_density(1))
          n = size(m % atom_density)
          err = 0
        end if
      end associate
    else
      err = E_OUT_OF_BOUNDS
    end if
  end function openmc_material_get_densities

!===============================================================================
! OPENMC_MATERIAL_GET_ID returns the ID of a material
!===============================================================================

  function openmc_material_get_id(index, id) result(err) bind(C)
    integer(C_INT32_T), value       :: index
    integer(C_INT32_T), intent(out) :: id
    integer(C_INT) :: err

    if (index >= 1 .and. index <= size(materials)) then
      id = materials(index) % id
      err = 0
    else
      err = E_OUT_OF_BOUNDS
    end if
  end function openmc_material_get_id

!===============================================================================
! OPENMC_MATERIAL_SET_DENSITY sets the total density of a material in atom/b-cm
!===============================================================================

  function openmc_material_set_density(index, density) result(err) bind(C)
    integer(C_INT32_T), value, intent(in) :: index
    real(C_DOUBLE), value, intent(in) :: density
    integer(C_INT) :: err

    err = E_UNASSIGNED
    if (index >= 1 .and. index <= size(materials)) then
      associate (m => materials(index))
        err = m % set_density(density)
      end associate
    else
      err = E_OUT_OF_BOUNDS
    end if
  end function openmc_material_set_density

!===============================================================================
! OPENMC_MATERIAL_SET_DENSITIES sets the densities for a list of nuclides in a
! material. If the nuclides don't already exist in the material, they will be
! added
!===============================================================================

  function openmc_material_set_densities(index, n, name, density) result(err) bind(C)
    integer(C_INT32_T), value, intent(in) :: index
    integer(C_INT), value, intent(in) :: n
    type(C_PTR),    intent(in) :: name(n)
    real(C_DOUBLE), intent(in) :: density(n)
    integer(C_INT) :: err

    integer :: i
    character(C_CHAR), pointer :: string(:)
    character(len=:, kind=C_CHAR), allocatable :: name_

    if (index >= 1 .and. index <= size(materials)) then
      associate (m => materials(index))
        ! If nuclide/density arrays are not correct size, reallocate
        if (n /= size(m % nuclide)) then
          deallocate(m % nuclide, m % atom_density, m % p0)
          allocate(m % nuclide(n), m % atom_density(n), m % p0(n))
        end if

        do i = 1, n
          ! Convert C string to Fortran string
          call c_f_pointer(name(i), string, [10])
          name_ = to_lower(to_f_string(string))

          if (.not. nuclide_dict % has_key(name_)) then
            err = openmc_load_nuclide(string)
            if (err < 0) return
          end if

          m % nuclide(i) = nuclide_dict % get_key(name_)
          m % atom_density(i) = density(i)
        end do
        m % n_nuclides = n

        ! Set isotropic flags to flags
        m % p0(:) = .false.

        ! Set total density to the sum of the vector
        err = m % set_density(sum(density))

        ! Assign S(a,b) tables
        call m % assign_sab_tables()
      end associate
    else
      err = E_OUT_OF_BOUNDS
    end if

  end function openmc_material_set_densities

!===============================================================================
! OPENMC_NUCLIDE_NAME returns the name of a nuclide with a given index
!===============================================================================

  function openmc_nuclide_name(index, name) result(err) bind(C)
    integer(C_INT), value, intent(in) :: index
    type(c_ptr), intent(out) :: name
    integer(C_INT) :: err

    character(C_CHAR), pointer :: name_

    err = E_UNASSIGNED
    if (allocated(nuclides)) then
      if (index >= 1 .and. index <= size(nuclides)) then
        name_ => nuclides(index) % name(1:1)
        name = C_LOC(name_)
        err = 0
      else
        err = E_OUT_OF_BOUNDS
      end if
    else
      err = E_NUCLIDE_NOT_ALLOCATED
    end if
  end function openmc_nuclide_name

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

end module openmc_api
