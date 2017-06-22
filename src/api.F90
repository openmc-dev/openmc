module openmc_api

  use, intrinsic :: ISO_C_BINDING

  use constants,       only: K_BOLTZMANN
  use eigenvalue,      only: k_sum
  use finalize,        only: openmc_finalize
  use geometry,        only: find_cell
  use global
  use initialize,      only: openmc_init
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
  public :: openmc_material_get_densities
  public :: openmc_material_set_density
  public :: openmc_plot_geometry
  public :: openmc_reset
  public :: openmc_run
  public :: openmc_set_density
  public :: openmc_set_temperature
  public :: openmc_tally_results

contains

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

end module openmc_api
