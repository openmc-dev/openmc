module openmc_api

  use, intrinsic :: ISO_C_BINDING

  use hdf5, only: HID_T

  use constants,       only: K_BOLTZMANN
  use eigenvalue,      only: k_sum
  use finalize,        only: openmc_finalize
  use geometry,        only: find_cell
  use global
  use hdf5_interface
  use message_passing, only: master
  use initialize,      only: openmc_init
  use input_xml,       only: assign_0K_elastic_scattering, check_data_version
  use particle_header, only: Particle
  use plot,            only: openmc_plot_geometry
  use simulation,      only: openmc_run
  use volume_calc,     only: openmc_calculate_volumes

  private
  public :: openmc_calculate_volumes
  public :: openmc_cell_set_temperature
  public :: openmc_finalize
  public :: openmc_find
  public :: openmc_init
  public :: openmc_material_add_nuclide
  public :: openmc_material_get_densities
  public :: openmc_material_set_density
  public :: openmc_load_nuclide
  public :: openmc_plot_geometry
  public :: openmc_reset
  public :: openmc_run
  public :: openmc_set_density
  public :: openmc_set_temperature
  public :: openmc_tally_results

contains

!===============================================================================
! OPENMC_CELL_SET_TEMPERATURE sets the temperature of a cell
!===============================================================================

  function openmc_cell_set_temperature(id, T) result(err) bind(C)
    integer(C_INT), value, intent(in) :: id  ! id of cell
    real(C_DOUBLE), value, intent(in) :: T
    integer(C_INT) :: err

    integer :: i

    err = -1
    if (allocated(cells)) then
      if (cell_dict % has_key(id)) then
        i = cell_dict % get_key(id)
        associate (c => cells(i))
          if (allocated(c % sqrtkT)) then
            c % sqrtkT(:) = sqrt(K_BOLTZMANN * T)
            err = 0
          end if
        end associate
      end if
    end if
  end function openmc_cell_set_temperature

!===============================================================================
! OPENMC_FIND determines the ID or a cell or material at a given point in space
!===============================================================================

  function openmc_find(xyz, rtype) result(id) bind(C)
    real(C_DOUBLE), intent(in)        :: xyz(3) ! Cartesian point
    integer(C_INT), intent(in), value :: rtype  ! 1 for cell, 2 for material
    integer(C_INT) :: id

    logical :: found
    type(Particle) :: p

    call p % initialize()
    p % coord(1) % xyz(:) = xyz
    p % coord(1) % uvw(:) = [ZERO, ZERO, ONE]
    call find_cell(p, found)

    id = -1
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
    end if
  end function openmc_find

!===============================================================================
! OPENMC_LOAD_NUCLIDE loads a nuclide from the cross section library
!===============================================================================

  function openmc_load_nuclide(name) result(err) bind(C)
    character(kind=C_CHAR) :: name(*)
    integer(C_INT) :: err

    integer :: n
    integer(HID_T) :: file_id
    integer(HID_T) :: group_id
    character(:), allocatable :: name_
    real(8) :: minmax(2) = [ZERO, INFINITY]
    type(VectorReal) :: temperature
    type(Nuclide), allocatable :: new_nuclides(:)

    ! Copy array of C_CHARs to normal Fortran string
    name_ = to_f_string(name)

    err = -1
    if (.not. nuclide_dict % has_key(to_lower(name_))) then
      if (library_dict % has_key(to_lower(name_))) then
        ! allocate extra space in nuclides array
        n = n_nuclides_total
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
        n_nuclides_total = n

        ! Assign resonant scattering data
        if (res_scat_on) call assign_0K_elastic_scattering(nuclides(n))

        ! Initialize nuclide grid
        call nuclides(n) % init_grid(energy_min_neutron, &
             energy_max_neutron, n_log_bins)

        err = 0
      else
        err = -2
      end if
    end if

  end function openmc_load_nuclide

!===============================================================================
! OPENMC_MATERIAL_ADD_NUCLIDE
!===============================================================================

  function openmc_material_add_nuclide(id, name, density) result(err) bind(C)
    integer(C_INT), value, intent(in) :: id
    character(kind=C_CHAR) :: name(*)
    real(C_DOUBLE), value, intent(in) :: density
    integer(C_INT) :: err

    integer :: i, j, k, n
    integer :: err2
    real(8) :: awr
    integer, allocatable :: new_nuclide(:)
    real(8), allocatable :: new_density(:)
    character(:), allocatable :: name_

    name_ = to_f_string(name)

    err = -1
    if (allocated(materials)) then
      if (material_dict % has_key(id)) then
        i = material_dict % get_key(id)
        associate (m => materials(i))
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
            err2 = openmc_load_nuclide(name)

            if (err2 /= -2) then
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

              err = 0
            end if
          end if
        end associate
      end if
    end if

  end function openmc_material_add_nuclide

!===============================================================================
! OPENMC_MATERIAL_GET_DENSITIES returns an array of nuclide densities in a
! material
!===============================================================================

  function openmc_material_get_densities(id, ptr) result(n) bind(C)
    integer(C_INT), intent(in), value :: id
    type(C_PTR),    intent(out) :: ptr
    integer(C_INT) :: n

    ptr = C_NULL_PTR
    n = 0
    if (allocated(materials)) then
      if (material_dict % has_key(id)) then
        i = material_dict % get_key(id)
        associate (m => materials(i))
          if (allocated(m % atom_density)) then
            ptr = C_LOC(m % atom_density(1))
            n = size(m % atom_density)
          end if
        end associate
      end if
    end if
  end function openmc_material_get_densities

!===============================================================================
! OPENMC_MATERIAL_SET_DENSITY sets the total density of a material in atom/b-cm
!===============================================================================

  function openmc_material_set_density(id, density) result(err) bind(C)
    integer(C_INT), value, intent(in) :: id
    real(C_DOUBLE), value, intent(in) :: density
    integer(C_INT) :: err

    integer :: i

    err = -1
    if (allocated(materials)) then
      if (material_dict % has_key(id)) then
        i = material_dict % get_key(id)
        associate (m => materials(i))
          err = m % set_density(density, nuclides)
        end associate
      end if
    end if
  end function openmc_material_set_density

!===============================================================================
! OPENMC_RESET resets all tallies
!===============================================================================

  subroutine openmc_reset() bind(C)
    integer :: i

    do i = 1, size(tallies)
      tallies(i) % n_realizations = 0
      if (allocated(tallies(i) % results)) then
        tallies(i) % results(:, :, :) = ZERO
      end if
    end do

    ! Reset global tallies
    n_realizations = 0
    if (allocated(global_tallies)) then
      global_tallies(:, :) = ZERO
    end if
    k_col_abs = ZERO
    k_col_tra = ZERO
    k_abs_tra = ZERO
    k_sum(:) = ZERO

    ! Turn off tally flags
    tallies_on = .false.
    active_batches = .false.

    ! Clear active tally lists
    call active_analog_tallies % clear()
    call active_tracklength_tallies % clear()
    call active_current_tallies % clear()
    call active_collision_tallies % clear()
    call active_tallies % clear()

  end subroutine openmc_reset

!===============================================================================
! OPENMC_SET_DENSITY sets the density of a material at a given point
!===============================================================================

  function openmc_set_density(xyz, density) result(err) bind(C)
    real(C_DOUBLE), intent(in) :: xyz(3)
    real(C_DOUBLE), intent(in), value :: density
    integer(C_INT) :: err

    logical :: found
    type(Particle) :: p

    call p % initialize()
    p % coord(1) % xyz(:) = xyz
    p % coord(1) % uvw(:) = [ZERO, ZERO, ONE]
    call find_cell(p, found)

    err = -1
    if (found) then
      if (p % material /= MATERIAL_VOID) then
        associate (m => materials(p % material))
          err = m % set_density(density, nuclides)
        end associate
      end if
    end if
  end function openmc_set_density

!===============================================================================
! OPENMC_SET_TEMPERATURE sets the temperature of a cell at a given point
!===============================================================================

  function openmc_set_temperature(xyz, T) result(err) bind(C)
    real(C_DOUBLE), intent(in) :: xyz(3)
    real(C_DOUBLE), intent(in), value :: T
    integer(C_INT) :: err

    logical :: found
    type(Particle) :: p

    call p % initialize()
    p % coord(1) % xyz(:) = xyz
    p % coord(1) % uvw(:) = [ZERO, ZERO, ONE]
    call find_cell(p, found)

    err = -1
    if (found) then
      associate (c => cells(p % coord(p % n_coord) % cell))
        if (size(c % sqrtkT) > 1) then
          c % sqrtkT(p % cell_instance) = sqrt(K_BOLTZMANN * T)
        else
          c % sqrtkT(1) = sqrt(K_BOLTZMANN * T)
        end if
        err = 0
      end associate
    end if
  end function openmc_set_temperature

!===============================================================================
! OPENMC_TALLY_RESULTS returns a pointer to a tally results array along with its
! shape. This allows a user to obtain in-memory tally results from Python
! directly.
!===============================================================================

  subroutine openmc_tally_results(i, ptr, shape_) bind(C)
    integer(C_INT), intent(in), value :: i
    type(C_PTR),    intent(out) :: ptr
    integer(C_INT), intent(out) :: shape_(3)

    ptr = C_NULL_PTR
    if (allocated(tallies)) then
      if (i >= 1 .and. i <= size(tallies)) then
        if (allocated(tallies(i) % results)) then
          ptr = C_LOC(tallies(i) % results(1,1,1))
          shape_(:) = shape(tallies(i) % results)
        end if
      end if
    end if
  end subroutine openmc_tally_results

  function to_f_string(c_string) result(f_string)
    character(kind=C_CHAR), intent(in) :: c_string(*)
    character(:), allocatable :: f_string

    integer :: i, n

    ! Determine length of original string
    n = 0
    do while (c_string(n + 1) /= C_NULL_CHAR)
      n = n + 1
    end do

    ! Copy C string character by character
    allocate(character(len=n) :: f_string)
    do i = 1, n
      f_string(i:i) = c_string(i)
    end do
  end function to_f_string

end module openmc_api
